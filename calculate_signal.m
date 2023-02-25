% calculateSignal Calculates the echo signal
%
%  [Signal,AuxSignals,Signals] = calculateSignal(System,Method,Nuclei,Clusters)
%
% Input:
%   System   system structure
%   Method   method structure
%   Nuclei   nuclei structure
%   Clusters information about clusters
%
% Output:
%   Signal
%   AuxSignals
%   Signals

% Clusters = Clusters(cluster index , 1:size ,order)
% Clusters(cluster index , size > order ,order) = 0.

function [total_signal,auxiliary_signals,order_n_signals] ... 
       = calculate_signal(System,Method,Nuclei,clusters,OutputData)

if Method.use_calculate_signal_ckpt
  [total_signal,auxiliary_signals,order_n_signals] ...
    = calculate_signal_ckpt(System,Method,Nuclei,clusters,OutputData);
else
  [total_signal,auxiliary_signals,order_n_signals] ...
    = calculate_signal_default(System,Method,Nuclei,clusters);
end

end

function [total_signal,auxiliary_signals,order_n_signals] ... 
       = calculate_signal_default(System,Method,Nuclei,clusters)

% Extract physical constants.     
ge = System.ge;
geff = System.gMatrix(3,3); 
muB = System.muB;
muN = System.muN;
mu0 = System.mu0;
hbar = System.hbar;

% Extract experimental parameters.
magneticField = System.magneticField;
dt1 = System.dt(1);
dt2 = System.dt(2);
N1 = System.nPoints(1);
N2 = System.nPoints(2);


B1y =  System.RF.B1y;

% ENUM
FID = 1; HAHN = 2; CPMG = 3; CPMG_CONST = 4; CPMG_2D = 5; %HAHN_TR = 6;
CP_N = 7; UHRIG_N = 8;

% Define the theory to use at the given cluster size.
Theory = System.Theory;
theory = Theory(Method.order,:);

useMeanFields = theory(10);
if useMeanFields
  Nuclear_Dipole_z_Z = zeroDiag(Nuclei.Statistics.Nuclear_Dipole);
  Nuclear_Dipole_x_iy_Z = zeroDiag(Nuclei.Statistics.Nuclear_Dipole_x_iy_Z);
else
  Nuclear_Dipole_z_Z = [];
  Nuclear_Dipole_x_iy_Z = [];
end

if Theory(Method.order,10)
  Method_extraOrder = Method.extraOrder;
  maxSuperclusterSize = Method_extraOrder;
else
  Method_extraOrder = Method.order;
  maxSuperclusterSize = Method.order;
end

% Convert variable to gpu compatible forms
dimensionality = 1;
doTR = false;
switch System.experiment
  case 'FID'
    EXPERIMENT = FID;
    total_time = System.Time(end);
  case 'Hahn'
    EXPERIMENT = HAHN;
    total_time = 2*System.Time(end);
  case 'CPMG'
    EXPERIMENT = CPMG;
    total_time = 4*System.Time(end);
  case 'CPMG-const'
    EXPERIMENT = CPMG_CONST;
    total_time = 4*System.Time(end);
  case 'CPMG-2D'
    EXPERIMENT = CPMG_2D;
    total_time = 4*System.Time(end);
    dimensionality = 2;
  case 'Hahn-TR'
    EXPERIMENT = HAHN;
    total_time = 2*System.Time(end);
    doTR = true;
    B1y =  System.RF.B1y(1);
    B1x2 =  System.RF.B1x(2);
    B1y2 =  System.RF.B1y(2);
    nuRF2 = System.RF.nuRF(2);
  case 'CP_N'
    EXPERIMENT = CP_N;
    total_time = 2*System.nPulses*System.Time(end);
  case 'Uhrig_N'
    EXPERIMENT = UHRIG_N;
    total_time = 2*System.nPulses*System.Time(end);
  otherwise
    error('The experiment ''%s'' is not supported.',System.experiment);
end

B1x =  System.RF.B1x;
nuRF = System.RF.nuRF;

betaT = 2*pi*System.hbar/System.kT; % 1/Hz.

% Unpack spin operators.
Op = Nuclei.SpinOperators;
SpinXiXjOps = Nuclei.SpinXiXjOperators;

numberClusters = Nuclei.numberClusters(1:maxSuperclusterSize);
maxNumberClusters = max(numberClusters(1:maxSuperclusterSize));

% CluserArray(iCluster,:,clusterSize) = nuclear indices.
ClusterArray = zeros(maxNumberClusters,Method.order,Method.order);
 
% Change data type of Clusters to a 3D array.
for isize = 1:Method_extraOrder 
  for ii = 1:Nuclei.numberClusters(isize)
    ClusterArray(ii,1:isize,isize) = clusters{isize}(ii,:);
  end
end

lockRotors = System.Methyl.lockRotors;
methylMethod1 = System.Methyl.method==0;
% Initialize coherences and subcluster indices.
[cluster_signals,  ...
  rotationalMatrix_c1_m1, rotationalMatrix_d1_m1, ...
  rotationalMatrix_c2_m1, rotationalMatrix_c2_m2, ...
  rotationalMatrix_d2_m1, rotationalMatrix_d2_m2] ...
  = initializeCoherences(Method, numberClusters, ...
  Nuclei.rotationalMatrix,(N1+N2)^dimensionality,methylMethod1,...
  lockRotors,System.Methyl.methylMethylCoupling);


% Methyl Groups
isMethylCarbon = strcmp(Nuclei.Type,'CH3');
isMethylHydron = strcmp(Nuclei.Type,'CH3_1H');

[Nuclei.ZeemanStates,Nuclei.ZeemanSpinStates] = setRandomZeemanState(Nuclei);

for clusterSize = 1:Method.order
 
  % Find coherences
  numClusters = numberClusters(clusterSize);
  
  for iCluster = numClusters:-1:1
    Cluster = ClusterArray(iCluster,1:clusterSize,clusterSize); 
   
    if System.Methyl.method ~= 2 && any(isMethylHydron(Cluster) )
      error(['Error in calculateSignal(): ',...
        'A methyl hydron is misplaced.']);
    end
    if System.Methyl.method == 2 && any(isMethylCarbon(Cluster) )
      error(['Error in calculateSignal(): ',...
        'A methyl carbon is misplaced.']);
    end
    if(  length(unique(Cluster)) <  length( Cluster )  )
      error(['Error in calculateSignal(): ',...
        'Cluster constains duplicate indices.']);
    end
    
    % Check if iCluster is valid.
    if Cluster(1,1) == 0, continue; end
    
    % Determin the number of spins in the cluster
    thisClusterSize = clusterSize + 2*sum(isMethylCarbon(Cluster));
    
    % Initialize cluster
    thisCluster = zeros(1,thisClusterSize);
    
    % Initialize counters.
    thisIndex = 0;
    methyl_number = 0;
    
    % Loop over the cluster to assign methyl rotors.
    for ii = 1:clusterSize
      
      if isMethylCarbon(Cluster(ii))
        % Count the number of metyls.
        methyl_number = methyl_number + 1;      
        
        % Unpack spins from methyl pseudo-particle.
        if Method.useCentralSpinSystem

          hydrons = ...
            nonzeros(Nuclei.associatedParticleMatrix(:,Cluster(ii)+1 ))' -1;
        else
           hydrons = Cluster(ii) + [1,2,3];
        end
        thisCluster(thisIndex+1:thisIndex+numel(hydrons)) = hydrons; 
        thisIndex = thisIndex + numel(hydrons);
      else
        thisIndex = thisIndex + 1;
        thisCluster(thisIndex) = Cluster(ii);
      end

    end
    if System.Methyl.method ==1
      thisClusterSize = numel(thisCluster);
    end
    if(  length(unique(thisCluster)) <  length( thisCluster )  )

      error(['Error in calculateSignal(): ',...
        'This cluster constains duplicate indices.']);  
    end
    if System.Methyl.method ~= 2 && mod(sum(isMethylHydron(thisCluster) ),3) ~= 0
      error(['Error in calculateSignal(): ',...
        'This cluster contains a partial methyl group.']);  
    end
    if any(isMethylCarbon(thisCluster) )
      error(['Error in calculateSignal(): ',...
        'A methyl carbon is misplaced.']);
    end
    
    % Generate methyl permutations
    cyclicPermutation = thisCluster;
    if lockRotors || ~methylMethod1
      nPerm = 1;
    else
      cyclicPermutation = getCyclicPermutations(isMethylCarbon(Cluster));
      nPerm = size(cyclicPermutation,1);
      cyclicPermutation = thisCluster(cyclicPermutation);
    end
    
    % Decide if the cluster should be skipped:
    % check if all spin have the same I, and that if I >= 1 is enabled. 
    skipCluster =  (System.spinHalfOnly && Nuclei.StateMultiplicity(thisCluster(1)) > 2);
    skipCluster = skipCluster || (~all(Nuclei.StateMultiplicity(thisCluster) ...
      == Nuclei.StateMultiplicity(thisCluster(1))));
       
    if skipCluster
      cluster_signals{clusterSize}(:,iCluster) = 1.0;
    end
    

    
    % Select the appropriate spin operator.
    spinMultiplicity = Nuclei.StateMultiplicity(thisCluster(1));
    SpinOp = Op{spinMultiplicity}{thisClusterSize};
    if spinMultiplicity==3
      SpinXiXjOp = SpinXiXjOps{thisClusterSize};
    else
      SpinXiXjOp = [];
    end
    H_alphaMF = 0;
    H_betaMF = 0;
     
    % density matrix size
    SpinSpaceDim = prod(Nuclei.StateMultiplicity(thisCluster));
    HilbertSpaceDim = SpinSpaceDim*nPerm;
    
    switch nPerm
      case 1 % 0 rotors
        Hb = zeros(HilbertSpaceDim);
      case 3 % 1 rotor
        switch SpinSpaceDim
          case 8 % CH3
            Hb = rotationalMatrix_c1_m1;
          case 27 % CD3
            Hb = rotationalMatrix_d1_m1;
          case 16 % H, CH3
            Hb = rotationalMatrix_c2_m1;
          case 81 % D, CD3
            Hb = rotationalMatrix_d2_m1;
        end
      case 9 % 2 rotors  
        switch SpinSpaceDim
          case 64 % CH3, CH3 
            Hb = rotationalMatrix_c2_m2;
          case 729 % CD3, CD3
            Hb = rotationalMatrix_d2_m2;
        end
      otherwise
        Hb = zeros(HilbertSpaceDim);
    end
    Ha = Hb;
    
    for iave = 1:System.nStates(clusterSize)
      if useMeanFields
        
        ZeemanStates = Nuclei.ZeemanStates(iave,:);
        spinState = Nuclei.ZeemanSpinStates(iave,:);
        spinState(thisCluster) = 0;
        if(size(spinState,2)>1)
          spinState = spinState';
        end
        
        mean_Dipole_z_Z = Nuclear_Dipole_z_Z*spinState;
        mean_Dipole_x_iy_Z = Nuclear_Dipole_x_iy_Z*spinState;
        
        %       mean_Dipole_z_Z = mean_Dipole_z_Z;
        %       mean_Dipole_x_iy_Z = mean_Dipole_x_iy_Z;
      else
        mean_Dipole_z_Z = [];
        mean_Dipole_x_iy_Z = [];
      end
    
      % Loop over cyclic permutations.
      for iPerm = 1:nPerm
        
        % Permute rotor spins.
        thisCluster = cyclicPermutation(iPerm,:);
        
        % Get interaction tensors.
        [tensors,~] = pairwisetensors(Nuclei.Nuclear_g, ...
          Nuclei.Coordinates,thisCluster,Nuclei.Atensor,magneticField,ge,geff,...
          muB,muN,mu0,hbar,theory,B1x,B1y,nuRF, ...
          mean_Dipole_z_Z, mean_Dipole_x_iy_Z);

        if ~System.limitToSpinHalf
          qtensors = Nuclei.Qtensor(:,:,thisCluster);
        else
          qtensors = 0;
        end
        
        if doTR
          [tensors_TR,~] = pairwisetensors(Nuclei.Nuclear_g, ...
            Nuclei.Coordinates,thisCluster,Nuclei.Atensors,magneticField,ge,geff,...
            muB, muN, mu0, hbar,theory,B1x2,B1y2,nuRF2,...
            mean_Dipole_z_Z, mean_Dipole_x_iy_Z);
          if ~System.limitToSpinHalf
            qtensors = Nuclei.Qtensors(:,:,thisCluster);
          else
            qtensors = 0;
          end
        else
          tensors_TR = [];
        end
        
        
        
        % Build Hamiltonians.
        % [H_alpha,H_beta] = assembleHamiltonian(state_multiplicity,...
        %   tensors,SpinOp,Qtensors,SpinXiXjOp,...
        %   theory,...
        %   isMethyl,...
        %   methylMethod, ...
        %   methylTunnelingSplitting, ...
        %   methylID)
        [H_alpha,H_beta] = ...
          assembleHamiltonian(Nuclei.StateMultiplicity(thisCluster),...
          tensors,SpinOp,qtensors,SpinXiXjOp, theory,...
          isMethylHydron(thisCluster),...
          System.Methyl.method,...
          Nuclei.methylTunnelingSplitting(thisCluster,thisCluster),...
          Nuclei.MethylID(thisCluster) ...
          );
        
        if ~isempty(tensors_TR)
        [H_alpha_TR,H_beta_TR] = ...
          assembleHamiltonian(Nuclei.StateMultiplicity(thisCluster),...
          tensors,SpinOp,qtensors,SpinXiXjOp, theory,...
          isMethylHydron(thisCluster), ...
          System.Methyl.method, ...
          Nuclei.methylTunnelingSplitting(thisCluster,thisCluster),...
          Nuclei.MethylID(thisCluster));
        else
        H_alpha_TR = [];
        H_beta_TR = [];  
        end
        
     
          
        if System.useThermalEnsemble
          densityMatrix = [];
        else
          densityMatrix = getDensityMatrix(ZeemanStates(thisCluster),...
            Nuclei.StateMultiplicity(thisCluster),thisCluster);
        end
            
        % Add mean field term to the Hamiltonian.
        Halpha = H_alpha + H_alphaMF;
        Hbeta = H_beta + H_betaMF;
       
        % Determine which block in the rotational Hamiltonian
        % this configuration is assigned to.
        hRange = (iPerm-1)*SpinSpaceDim+(1:SpinSpaceDim);
        
        % Add configuration Hamiltonian to the full Hamiltonian.
        Hb(hRange,hRange) = Hbeta;
        Ha(hRange,hRange) = Halpha;
        
        if ~isempty(H_alpha_TR)
          Hb_TR(hRange,hRange) = H_beta_TR;
          Ha_TR(hRange,hRange) = H_alpha_TR;
        else
          Hb_TR = [];
          Ha_TR = [];
        end
      end
 
      % Calculate cluster coherence.
      new_sig = propagate(total_time,N1+N2,dt1,dt2,N1,Hb,Ha,Hb_TR,Ha_TR, ...
        EXPERIMENT,densityMatrix, System.useThermalEnsemble, betaT, ...
        System.nPulses);

      cluster_signals{clusterSize}(:,iCluster) ...
        = cluster_signals{clusterSize}(:,iCluster)  +  new_sig;
      
      
    end
    if abs(cluster_signals{clusterSize}(1,iCluster)-1) >1e-12
    cluster_signals{clusterSize}(:,iCluster) = ...
      cluster_signals{clusterSize}(:,iCluster)...
      ./cluster_signals{clusterSize}(1,iCluster);
    end
  end
end

% Unused coherences need to be re-initialized to unity.
for clusterSize = 1:Method.order
  unusedClusters = find(  cluster_signals{clusterSize}(1,:)==0  )';
  if ~isempty(unusedClusters)
    cluster_signals{clusterSize}(:,unusedClusters) = 1;
  end
end

% Calculate signal
%-------------------------------------------------------------------------------
auxiliary_signals = get_auxiliary_signals(cluster_signals,clusters);

order_n_signals = do_CCE_product(auxiliary_signals);

total_signal = order_n_signals(Method.order,:);
 
%-------------------------------------------------------------------------------
end


% ==============================================================================
% Generate Density Matrix
% ==============================================================================
function densityMatrix = getDensityMatrix(states,multiplicities,Cluster)

clustersize = length(Cluster);

prod_state = states(clustersize);

offset_factor = multiplicities(clustersize);
for iSpin = clustersize-1:-1:1 
  prod_state = prod_state + offset_factor*(states(iSpin) - 1);
  offset_factor = offset_factor*multiplicities(iSpin);
end
densityMatrix = zeros(offset_factor);
densityMatrix(prod_state,prod_state) = 1;

end

% ========================================================================
% Propagate Function
% ========================================================================
function Signal = propagate(total_time,Npoints,dt,dt2,N1,...
  Ham_beta,Ham_alpha,...
  Ham_beta_TR,Ham_alpha_TR ,...
  EXPERIMENT, densityMatrix, useThermalEnsemble, betaT, nPulses)

% ENUM
FID = 1; HAHN = 2; CPMG = 3; CPMG_CONST = 4; CPMG_2D = 5; HAHN_TR = 6;
CP_N = 7; UHRIG_N = 8;

% Ensure sub-Hamiltonians are Hermitian.
Ham_beta = (Ham_beta+Ham_beta')/2; 
Ham_alpha = (Ham_alpha+Ham_alpha')/2;
 
% Get density matrix.
if useThermalEnsemble
  DensityMatrix = propagator_eig((Ham_alpha+Ham_beta)/2,-1i*betaT,0);
else
  DensityMatrix = densityMatrix;
end

% Vectorize density matrix.
vecDensityMatrixT = reshape(DensityMatrix.',1,[])/trace(DensityMatrix);

t_Uhrig = 0;

% Generate propagators for time steps dt and dt2
[dU_beta,dU_beta2] = propagator_eig(Ham_beta,dt,dt2);
[dU_alpha,dU_alpha2] = propagator_eig(Ham_alpha,dt,dt2);

doTR = ~isempty(Ham_beta_TR);
if doTR  
  % Ensure subHamiltonians are Hermitian
  Ham_beta_TR = (Ham_beta_TR+Ham_beta_TR')/2;
  Ham_alpha_TR = (Ham_alpha_TR+Ham_alpha_TR')/2;
  % Generate propagators for time steps dt and dt2
  [dU_beta_TR, dU_beta2_TR] = propagator_eig(Ham_beta_TR,dt,dt2);
  [dU_alpha_TR, dU_alpha2_TR] = propagator_eig(Ham_alpha_TR,dt, dt2);
end

% Initialize propagators
nStates = length(Ham_beta);
U_beta = eye(nStates);
U_alpha = eye(nStates);

% Initialize signal
if EXPERIMENT==CPMG_2D
  v = ones(Npoints,Npoints);
else
  v = ones(Npoints,1);
end

% Loop over time points
for iTime = 1:Npoints
  switch EXPERIMENT
    
    case FID
      U = U_beta'*U_alpha;
      v(iTime) = vecDensityMatrixT*U(:);
      
    case HAHN
      %       U = U_beta'*U_alpha'*U_beta*U_alpha;
      U = (U_alpha*U_beta)'*U_beta*U_alpha;
      v(iTime) = vecDensityMatrixT*U(:);
      
    case CPMG
      
      U = U_alpha'*U_beta'  * ...
        (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha*U_beta;
      v(iTime) = vecDensityMatrixT*U(:);
          
    case CPMG_CONST
      
      % THIS NEEDS TO BE UPDATED TO USE dt2 WHEN iTime > N1.
      [U_beta, U_beta_2] = propagator_eig(...
        Ham_beta,(iTime-1)*dt, total_time/4-(iTime-1)*dt);
      [U_alpha, U_alpha_2] = propagator_eig(...
        Ham_alpha,(iTime-1)*dt, total_time/4-(iTime-1)*dt);
      
      U = U_alpha_2'*U_beta_2'  *  ...
        (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
      
      v(iTime) = vecDensityMatrixT*U(:);
       
    case CPMG_2D
      
      U_beta_2 = eye(nStates);
      U_alpha_2 = eye(nStates);
      
      % Loop over second experimental dimension.
      for jTime = 1:iTime
        
        % Generate time dependent detection operator.
        U = U_alpha_2'*U_beta_2'  *  ...
          (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
        
        v(iTime,jTime) = vecDensityMatrixT*U(:);
        v(jTime,iTime) = v(iTime,jTime)';
        
        % Increment propagator.
        if jTime < N1
          U_beta_2 = dU_beta*U_beta_2;
          U_alpha_2 = dU_alpha*U_alpha_2;
        else
          U_beta_2 = dU_beta2*U_beta_2;
          U_alpha_2 = dU_alpha2*U_alpha_2;
        end
        
      end
    
    case HAHN_TR
      U = U_beta'*U_beta_TR'*U_alpha'*U_alpha_TR' ...
         * U_beta_TR*U_beta*U_alpha_TR*U_alpha;
      
      v(iTime) = vecDensityMatrixT*U(:);
    
    case CP_N
      U_aa = U_alpha*U_alpha;
      U_bb = U_beta*U_beta;

      exponent = fix((nPulses-1)/2);
      AABB = (U_aa*U_bb)^exponent;
      BBAA = (U_bb*U_aa)^exponent;
      if mod(nPulses,2)==0
        U = ( U_alpha*BBAA*U_bb*U_alpha )' ...
               *( U_beta*AABB*U_aa*U_beta  );
      else
        U = (  U_beta*AABB*U_alpha )' ...
            *( U_alpha*BBAA*U_beta  );
        
      end

      v(iTime) = vecDensityMatrixT*U(:); 

    case UHRIG_N
      % t_i = (2*tau)*sin^2[ (pi*i)/(2*(N+1)) ].

      if iTime< N1
        t_Uhrig = t_Uhrig + dt*2*nPulses;
      else 
        t_Uhrig = t_Uhrig + dt2*2*nPulses;
      end
        
      v(iTime) = getUhrigN(t_Uhrig, Ham_beta,Ham_alpha,...
      vecDensityMatrixT, nPulses);

  end
  
  % Increment propagators
  if iTime < N1
    U_beta = dU_beta*U_beta;
    U_alpha = dU_alpha*U_alpha;
  else
    U_beta = dU_beta2*U_beta;
    U_alpha = dU_alpha2*U_alpha;
  end


  if doTR
    if iTime < N1
      U_beta_TR = dU_beta_TR*U_beta_TR;
      U_alpha_TR = dU_alpha_TR*U_alpha_TR;
    else
      U_beta_TR = dU_beta2_TR*U_beta_TR;
      U_alpha_TR = dU_alpha2_TR*U_alpha_TR;
    end
  end

end

% Return signal as a vector
if EXPERIMENT == CPMG_2D
    Signal = reshape(v.',1,[]);
else
    Signal = v;
end

if any(abs(v) - 1 > 1e-9)
  error('Coherence error: coherence cannot be larger than 1.');
end

end


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function v_iTime = getUhrigN(total_time, Hamiltonian_beta,Hamiltonian_alpha,...
      vecDensityMatrixT, nPulses)
% t_i = T*sin^2[ (pi*i)/(2*(N+1)) ].

% delays = T*1/2*(cos(pi*([indices,N+1]-1)/( (N+1) ) ) ...
%        - cos(pi*[indices,N+1]/( (N+1) ) ));

% Since the Uhrig pulse sequence is symmetric about the  half-way point,
% only half the delays need to explicitly propogated.


if mod(nPulses,2)==0
  % Since the Uhrig pulse sequence is symmetric about the  half-way point,
  % only half the delays need to explicitly propogated.
  indices = 0:ceil( (nPulses+1)/2) ;
  delta_t = total_time * diff(sin(indices*pi/( 2*(nPulses+1) ) ).^2);

  % Split the last delat in two parts if needed.
  delta_t(end) = delta_t(end)/2;

  U_ket0 = 1;
  U_ket1 = 1;
  U_bra0 = 1;
  U_bra1 = 1;
  useHbeta = true;
  for t = delta_t
    if useHbeta
      H_ket = Hamiltonian_beta;
      H_bra = Hamiltonian_alpha;
    else
      H_ket = Hamiltonian_alpha;
      H_bra = Hamiltonian_beta;
    end
    useHbeta = ~useHbeta;

    U_prop = propagator_eig(H_ket,t);
    U_ket0 = U_prop*U_ket0;
    U_ket1 = U_ket1*U_prop;

    U_prop = propagator_eig(H_bra,t);
    U_bra0 = U_prop*U_bra0;
    U_bra1 = U_bra1*U_prop;
  end

else
  
%   indices = 0:ceil(nPulses+1);
%   delta_t = total_time * diff(sin(indices*pi/( 2*(nPulses+1) ) ).^2);
indices = 0:( (nPulses+1)/2) ;
delta_t = total_time * diff(sin(indices*pi/( 2*(nPulses+1) ) ).^2);

  U_ket0 = 1;
  U_ket1 = 1;
  U_bra0 = 1;
  U_bra1 = 1;

  useHbeta = true;
  for t = delta_t
    if useHbeta
      H_ket = Hamiltonian_beta;
      H_bra = Hamiltonian_alpha;
    else
      H_ket = Hamiltonian_alpha;
      H_bra = Hamiltonian_beta;
    end
    useHbeta = ~useHbeta;

    U_prop = propagator_eig(H_ket,t);
    U_ket0 = U_prop*U_ket0;
    U_bra1 = U_bra1*U_prop;

    U_prop = propagator_eig(H_bra,t);
    U_bra0 = U_prop*U_bra0;
    U_ket1 = U_ket1*U_prop;
  end

end
U_ = (U_bra1*U_bra0)' * (U_ket1*U_ket0);

v_iTime = vecDensityMatrixT*U_(:);

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


% ========================================================================
% Calculate propagator using diagonalization
% ========================================================================

function [U1,U2] = propagator_eig(Ham,dt1,dt2)
[EigenVectors,EigenValues] = eig(Ham);
Udiag1 = exp(-2i*pi*diag(EigenValues)*dt1);
U1 = EigenVectors*diag(Udiag1)*EigenVectors';

if nargin>2 && dt2>0
  Udiag2 = exp(-2i*pi*diag(EigenValues)*dt2);
  U2 = EigenVectors*diag(Udiag2)*EigenVectors';
else
  U2 = 1;
end
end



% ==============================================================================
% Generate cyclic permutation
% ==============================================================================
function cycPerm = getCyclicPermutations(boolCluster)

N = length(boolCluster);
numberRotors = sum(boolCluster);

% no methyls
if numberRotors == 0
  cycPerm = 1:N;
  return;
end

% 1-clusters
if N == 1
  cycPerm = [1,2,3;2,3,1;3,1,2];
  return;
end

% 2-clusters
if N == 2
  if numberRotors == 2
    cycPerm = [1,2,3,4,5,6;...
      2,3,1,4,5,6; ...
      3,1,2,4,5,6; ...
      1,2,3,5,6,4;...
      2,3,1,5,6,4; ...
      3,1,2,5,6,4; ...
      1,2,3,6,4,5;...
      2,3,1,6,4,5; ...
      3,1,2,6,4,5];
  elseif boolCluster(1)
    cycPerm = [1,2,3,4;2,3,1,4;3,1,2,4];
  elseif boolCluster(2)
    cycPerm = [1,2,3,4;1,3,4,2;1,4,2,3];
  end
  return;
end

% 3-clusters
if N==3
  if numberRotors == 3
    cycPerm = [1,2,3, 4,5,6, 7,8,9;...
      2,3,1, 4,5,6, 7,8,9; ...
      3,1,2, 4,5,6, 7,8,9; ...
      1,2,3, 5,6,4, 7,8,9;...
      2,3,1, 5,6,4, 7,8,9; ...
      3,1,2, 5,6,4, 7,8,9; ...
      1,2,3, 6,4,5, 7,8,9;...
      2,3,1, 6,4,5, 7,8,9; ...
      3,1,2, 6,4,5, 7,8,9; ...
      2,3,1, 4,5,6, 7,8,9; ...
      3,1,2, 4,5,6, 7,8,9; ...
      1,2,3, 5,6,4, 9,7,8;...
      2,3,1, 5,6,4, 9,7,8; ...
      3,1,2, 5,6,4, 9,7,8; ...
      1,2,3, 6,4,5, 9,7,8; ...
      2,3,1, 6,4,5, 9,7,8; ...
      3,1,2, 6,4,5, 9,7,8; ...
      2,3,1, 4,5,6, 9,7,8; ...
      3,1,2, 4,5,6, 9,7,8; ...
      1,2,3, 5,6,4, 8,9,7;...
      2,3,1, 5,6,4, 8,9,7; ...
      3,1,2, 5,6,4, 8,9,7; ...
      1,2,3, 6,4,5, 8,9,7;...
      2,3,1, 6,4,5, 8,9,7; ...
      3,1,2, 6,4,5, 8,9,7];
  elseif numberRotors == 2
    cycPerm = [1,2,3,4,5,6;...
      2,3,1,4,5,6; ...
      3,1,2,4,5,6; ...
      1,2,3,5,6,4;...
      2,3,1,5,6,4; ...
      3,1,2,5,6,4; ...
      1,2,3,6,4,5;...
      2,3,1,6,4,5; ...
      3,1,2,6,4,5];
    
    if ~boolCluster(3)
      cycPerm = [cycPerm,7*ones(9,1)];
    elseif ~boolCluster(2)
      cycPerm(:,4:6) = cycPerm(:,4:6) + 1;
      cycPerm = [cycPerm(:,1:3), 4*ones(9,1), cycPerm(:,4:6)];
    else
      cycPerm = cycPerm + 1;
      cycPerm = [1*ones(9,1),cycPerm];
    end
    
  elseif numberRotors == 1
    
  cycPerm = [1,2,3;2,3,1;3,1,2];
  if boolCluster(1)
    cycPerm = [cycPerm,[4,5].*ones(3,2)];
  elseif boolCluster(2)
    cycPerm = cycPerm + 1;
    cycPerm = [1*ones(3,1), cycPerm, 5*ones(3,1)];
  else
    cycPerm = cycPerm + 2;
    cycPerm = [[1,2].*ones(3,2),cycPerm];
  end
  end
  
  return; 
end

error('Cannot calculate cyclic permutations.');

end



% ========================================================================
% Initialize coherences
% ========================================================================

function [cluster_signals, ...
  rotationalMatrix_c1_m1, rotationalMatrix_d1_m1, ...
  rotationalMatrix_c2_m1, rotationalMatrix_c2_m2, ...
  rotationalMatrix_d2_m1, rotationalMatrix_d2_m2] ...
  = initializeCoherences(Method, numberClusters, ...
  Nuclei_rotationalMatrix, nt,  includeMethyls,lockRotor, ...
  Methyl_methylMethylCoupling)

cluster_signals = cell(Method.order,1);
for iclusterSize = 1:Method.order
  
  cluster_signals{iclusterSize} = zeros(nt,numberClusters(iclusterSize));
    

  % Initialize
  rotationalMatrix_c1_m1 = [];
  rotationalMatrix_d1_m1 = [];
  rotationalMatrix_c2_m1 = [];
  rotationalMatrix_c2_m2 = [];
  rotationalMatrix_d2_m1 = [];
  rotationalMatrix_d2_m2 = [];

  switch iclusterSize


    case 1
      if includeMethyls && ~lockRotor
        rotationalMatrix_c1_m1 = Nuclei_rotationalMatrix{1,1};
      else
        rotationalMatrix_c1_m1 = [];
        rotationalMatrix_d1_m1 = [];
      end

    case 2
      if includeMethyls && ~lockRotor
        rotationalMatrix_c2_m1 = Nuclei_rotationalMatrix{2,1};
        if Methyl_methylMethylCoupling
          rotationalMatrix_c2_m2 = Nuclei_rotationalMatrix{2,2};
        else
          rotationalMatrix_c2_m2 = [];
        end
      else
        rotationalMatrix_c2_m1 = [];
        rotationalMatrix_c2_m2 = [];
        rotationalMatrix_d2_m1 = [];
        rotationalMatrix_d2_m2 = [];
      end
  end

end

end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%-------------------------------------------------------------------------------
function subcluster_data = build_subcluster_data(cluster)

cluster_size = numel(cluster);

include_indices = include_vertices_in_which_subclusters(cluster_size);

buff_size = (cluster_size + 2)*2^(cluster_size-1) - cluster_size - 2;
subcluster_data = zeros(1,buff_size);
idx = 1;

for ii = 1:size(include_indices,1)
  subcluster = cluster(include_indices(ii,:));
  subcluster_size = numel(subcluster);

  subcluster_data(idx) = subcluster_size;
  idx = idx + 1;

  subcluster_data(idx:idx+subcluster_size-1) = subcluster;
  idx = idx + subcluster_size;

end


% Check resutls.
expected_size = (cluster_size + 2)*2^(cluster_size-1) - cluster_size - 2;
assert(numel(subcluster_data)==expected_size);
if cluster_size==3
  a = cluster(1); b = cluster(2); c = cluster(3);
  expected_data = [2,b,c, 2,a,c, 1,c, 2,a,b, 1,b, 1,a];
  assert(all(subcluster_data==expected_data))
end

end
%-------------------------------------------------------------------------------
function include_indices = include_vertices_in_which_subclusters(cluster_size)
  n_subclusters = 2^cluster_size;

  include_indices = false(n_subclusters-2,cluster_size);
  increments = 2.^((1:cluster_size)-1);
  keep = true(1,cluster_size);
  for ii = 2:n_subclusters-1
    sele = mod(ii-1,increments)==0;
    keep(sele) = ~keep(sele);
    include_indices(ii-1,keep) = true;
  end
end
%-------------------------------------------------------------------------------
% This function is not thread safe.  Do not use it inside parfor loops.
function auxiliary_signals = get_auxiliary_signals(cluster_signals,clusters)



cluster_map = build_cluster_map(clusters);

auxiliary_signals = cluster_signals;

for cluster_size = 1:numel(clusters)
  num_clusters = size(clusters{cluster_size},1);
  for icluster = 1:num_clusters

    % get_subclusters(cluster) returns all subcluster sizes and subclusters:
    % For the example cluster [a,b,c], build_subcluster_data([a,b,c]) 
    % returns [2,b,c, 2,a,c, 1,c, 2,a,b, 1,b, 1,a];
    % note that each subcluster is preceded by its size.
    cluster = clusters{cluster_size}(icluster,:);
    subcluster_data = build_subcluster_data(cluster); 

    idx = 1;
    while idx < numel(subcluster_data)
      subcluster_size = subcluster_data(idx);
      assert(subcluster_size<cluster_size);
      idx = idx + 1;
      subcluster = zeros(1,subcluster_size);
      for ivertex = 1:subcluster_size
        subcluster(ivertex) = subcluster_data(idx);
        idx = idx + 1;
      end

      %key = num2str(subcluster);
      key = sprintf('%d,',subcluster);
      %key = {subcluster}; % TODO: Faster, but incorect, fix.
      if ~cluster_map.isKey(key)
%         assert(sum(all(subcluster==clusters{subcluster_size},2))==0);
        continue; 
      end



%       assert(sum(all(subcluster==clusters{subcluster_size},2))==1);
     
      subcluster_index = cluster_map(key);
      
      auxiliary_signals{cluster_size}(:,icluster)...
        = auxiliary_signals{cluster_size}(:,icluster)...
        ./auxiliary_signals{subcluster_size}(:,subcluster_index);
    end

  end
end
end
%-------------------------------------------------------------------------------
function cluster_map = build_cluster_map(clusters)

cluster_map = dictionary(string([]),[]);
% cluster_map = dictionary({},[]);  % TODO: Faster, but incorect, fix.

for isize = 1:numel(clusters)
  for ii = 1:size(clusters{isize},1)
    cluster = clusters{isize}(ii,:);
    key = sprintf('%d,',cluster);
    cluster_map(key) = ii;
%     cluster_map({cluster}) = ii;
  end
end

end
%-------------------------------------------------------------------------------
function order_n_signals = do_CCE_product(auxiliary_signals)

max_order = numel(auxiliary_signals);
nt = size(auxiliary_signals{1},1);

order_n_signals = ones(nt,max_order);

for isize = 1:max_order
  order_n_signals(:,isize) = prod(auxiliary_signals{isize},2);
end

order_n_signals = cumprod(order_n_signals,2).';

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




