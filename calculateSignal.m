% Clusters = Clusters(cluster index , 1:size ,order)
% Clusters(cluster index , size > order ,order) = 0.

function [Signal, AuxiliarySignals, Signals] ... 
       = calculateSignal(System, Method, Nuclei,Clusters)

% Extract physical constants.     
ge=System.ge; geff=System.gMatrix(3,3); 
muB = System.muB;  muN = System.muN;  mu0 = System.mu0; hbar = System.hbar;

% Extract experimental parameters.
magneticField = System.magneticField;
timepoints = System.timepoints;
dt = System.dt;
t0 = System.t0;
dt2 = System.dt2;
Ndt = System.Ndt;
 
B1x =  System.RF.B1x;
B1y =  System.RF.B1y;
nuRF = System.RF.nuRF;
doTR = false;

% ENUM
FID = 1; HAHN = 2; CPMG = 3; CPMG_CONST = 4; CPMG_2D = 5; HAHN_TR = 6;

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

maxSpinClusterSize = Nuclei.maxClusterSize;

if Theory(Method.order,10)
  Method_extraOrder = Method.extraOrder;
  maxSuperclusterSize = Method_extraOrder;
else
  Method_extraOrder = Method.order;
  maxSuperclusterSize = Method.order;
end

if strcmp(Method.method,'HD-CCE')
  doHDCCE = true;
else
  doHDCCE = false;
end

if strcmp(Method.method,'OffsetCCE')
  doOffsetCCE = true;
else
  doOffsetCCE = false;
end

% Convert variable to gpu compatible forms
dimensionality = 1;
if strcmp(System.experiment,'FID')
  EXPERIMENT = FID;
  total_time = System.Time(end);
elseif strcmp(System.experiment,'Hahn')
  EXPERIMENT = HAHN;
  total_time = 2*System.Time(end);
elseif strcmp(System.experiment,'CPMG')
  EXPERIMENT = CPMG;
  total_time = 4*System.Time(end);
elseif strcmp(System.experiment,'CPMG-const')
  EXPERIMENT = CPMG_CONST;
  total_time = 4*System.Time(end);
elseif strcmp(System.experiment,'CPMG-2D')
  EXPERIMENT = CPMG_2D;
  total_time = 4*System.Time(end);
  dimensionality = 2;
elseif strcmp(System.experiment,'Hahn-TR')
  EXPERIMENT = HAHN;
  total_time = 2*System.Time(end);
  doTR = true;
  B1x =  System.RF.B1x(1);
  B1y =  System.RF.B1y(1);
  nuRF = System.RF.nuRF(1);
  B1x2 =  System.RF.B1x(2);
  B1y2 =  System.RF.B1y(2);
  nuRF2 = System.RF.nuRF(2);
else
  error('The experiment ''%s'' is not supported.',System.experiment);
end



betaT = 2*pi*System.hbar /System.kT; % 1/Hz.


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
    ClusterArray(ii,1:isize,isize) = Clusters{isize}(ii,:);
  end
end

lockRotors = System.Methyl.lockRotors;
methylMethod1 = System.Methyl.method==0;
% Initialize coherences and subcluster indices.
[Coherences, SubclusterIndices, ...
  rotationalMatrix_c1_m1, rotationalMatrix_d1_m1, ...
  rotationalMatrix_c2_m1, rotationalMatrix_c2_m2, ...
  rotationalMatrix_d2_m1, rotationalMatrix_d2_m2] ...
  = initializeCoherences(Method, numberClusters, ...
  Nuclei.rotationalMatrix,timepoints^dimensionality,methylMethod1,...
  lockRotors,System.Methyl.methylMethylCoupling);


% Methyl Groups
isMethylCarbon = strcmp(Nuclei.Type,'CH3');
isMethylHydron = strcmp(Nuclei.Type,'CH3_1H');


[Nuclei.ZeemanStates,Nuclei.ZeemanSpinStates] = ...
  setRandomZeemanState(Nuclei);
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
    
    % Determin the number of spins in the cluster.
    thisClusterSize = clusterSize + 2*sum(isMethylCarbon(Cluster));
    
    % Initialize cluster.
    thisCluster = zeros(1,thisClusterSize);
    
    % Initialize counters.
    thisIndex = 0;
    methyl_number = 0;
    
    % Loop over the cluster to assign methyl rotors.
    for ii = 1:clusterSize

      thisIndex = thisIndex + 1;
      
      if isMethylCarbon(Cluster(ii))
        % Count the number of metyls.
        methyl_number = methyl_number + 1;      
        
        % Unpack spins from methyl pseudo-particle.
        thisCluster(thisIndex:thisIndex+2) = Cluster(ii) + [1,2,3];
        thisIndex = thisIndex + 2;
      else
        thisCluster(thisIndex) = Cluster(ii);
      end

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
    
    
    % Generate methyl permutations.
    cyclicPermutation = thisCluster;
    if lockRotors || ~methylMethod1
      nPerm = 1;
    else
      cyclicPermutation = getCyclicPermutations(isMethylCarbon(Cluster));
      nPerm = size(cyclicPermutation,1);
      cyclicPermutation = thisCluster(cyclicPermutation);
    end
    
    % Decide if the cluster should be skipped:
    % check if all spin have the same I,
    % and that if I >= 1 is enabled. 
    skipCluster =  (System.spinHalfOnly && Nuclei.StateMultiplicity(thisCluster(1)) > 2);
    skipCluster = skipCluster || (~all(Nuclei.StateMultiplicity(thisCluster) ...
      == Nuclei.StateMultiplicity(thisCluster(1))));
   
    
    if skipCluster
      Coherences{clusterSize}(iCluster,:) = ones(size(Coherences{clusterSize}(iCluster,:) ));
    end
    
    % Set subcluster indices.
    SubclusterIndices{clusterSize}(:,:,iCluster) = findSubclusters_gpu(...
      ClusterArray,clusterSize,iCluster,clusterSize);
    
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
          case 729 % CH3, CH3
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
        [tensors,zeroIndex] = pairwisetensors_gpu(Nuclei.Nuclear_g, ...
          Nuclei.Coordinates,thisCluster,Nuclei.Atensor,magneticField, ge,geff,...
          muB, muN, mu0, hbar,theory,B1x,B1y,nuRF, ...
          mean_Dipole_z_Z, mean_Dipole_x_iy_Z);
        qtensors = Nuclei.Qtensor(:,:,thisCluster);
        
        if doTR        
          [tensors_TR,~] = pairwisetensors_gpu(Nuclei.Nuclear_g, ...
          Nuclei.Coordinates,thisCluster,Nuclei.Atensors,magneticField, ge,geff,...
          muB, muN, mu0, hbar,theory,B1x2,B1y2,nuRF2,...
          mean_Dipole_z_Z, mean_Dipole_x_iy_Z);
        qtensors = Nuclei.Qtensors(:,:,thisCluster);
        else
          tensors_TR = [];
        end
        
        
        
        % Build Hamiltonians.
        % [H_alpha,H_beta] = assembleHamiltonian_gpu(state_multiplicity,...
        %   tensors,SpinOp,Qtensors,SpinXiXjOp,...
        %   theory,...
        %   isMethyl,...
        %   methylMethod, ...
        %   methylTunnelingSplitting, ...
        %   methylID)
        [H_alpha,H_beta] = ...
          assembleHamiltonian_gpu(Nuclei.StateMultiplicity(thisCluster),...
          tensors,SpinOp,qtensors,SpinXiXjOp, theory,...
          isMethylHydron(thisCluster),...
          System.Methyl.method,...
          Nuclei.methylTunnelingSplitting(thisCluster,thisCluster),...
          Nuclei.MethylID(thisCluster) ...
          );
        
        if ~isempty(tensors_TR)
        [H_alpha_TR,H_beta_TR] = ...
          assembleHamiltonian_gpu(Nuclei.StateMultiplicity(thisCluster),...
          tensors,SpinOp,qtensors,SpinXiXjOp, theory,...
          isMethylHydron(thisCluster),System.Methyl.method, ...
          Nuclei.methylTunnelingSplitting(thisCluster,thisCluster),...
          Nuclei.MethylID(thisCluster));
        else
        H_alpha_TR = [];
        H_beta_TR = [];  
        end
        
     
          
        if System.useThermalEnsemble
          densityMatrix = [];
        else
          densityMatrix = getDensityMatrix(ZeemanStates(iave,thisCluster),...
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
      Coherences{clusterSize}(iCluster,:) = Coherences{clusterSize}(iCluster,:)  +  ...
        propagate(total_time,timepoints,dt,dt2,Ndt,Hb,Ha,Hb_TR,Ha_TR, ...
        EXPERIMENT,densityMatrix, System.useThermalEnsemble, betaT);
      
    end
    Coherences{clusterSize}(iCluster,:) = ...
      Coherences{clusterSize}(iCluster,:)/Coherences{clusterSize}(iCluster,1);
  end
end

% Unused coherences need to be re-initialized to unity.
for clusterSize = 1:Method.order
  unusedClusters = find(  Coherences{clusterSize}(:,1)==0  )';
  if ~isempty(unusedClusters)
    Coherences{clusterSize}(unusedClusters,:) = 1;
  end
end

% Calculate signal
%-------------------------------------------------------------------------------
[Signals, AuxiliarySignals] = doClusterCorrelationExpansion(...
  Coherences,ClusterArray,SubclusterIndices,...
  timepoints,dimensionality, Method.order,numberClusters);
 

Signal = Signals(Method.order,:);
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
function Signal = propagate(total_time,timepoints,dt,dt2,Ndt,...
  Hamiltonian_beta,Hamiltonian_alpha,...
  Hamiltonian_beta_TR,Hamiltonian_alpha_TR ,...
  EXPERIMENT, densityMatrix, useThermalEnsemble, betaT)

% ENUM
FID = 1; HAHN = 2; CPMG = 3; CPMG_CONST = 4; CPMG_2D = 5; HAHN_TR = 6;

% Ensure subHamiltinians are Hermitian.
Hamiltonian_beta =(Hamiltonian_beta+Hamiltonian_beta')/2; 
Hamiltonian_alpha =(Hamiltonian_alpha+Hamiltonian_alpha')/2;
 

% Get density matrix.
if useThermalEnsemble
  DensityMatrix = propagator_eig(...
    (Hamiltonian_alpha+Hamiltonian_beta)/2,-1i*betaT,0);
else
  DensityMatrix = densityMatrix;
end

% Vectorize density matrix.
vecDensityMatrixT = reshape(DensityMatrix.',1,[])/trace(DensityMatrix);

% Generate propagators for time element dt.
[dU_beta, dU_beta2] = propagator_eig(Hamiltonian_beta,dt,dt2);
[dU_alpha, dU_alpha2] = propagator_eig(Hamiltonian_alpha,dt, dt2);

% Find the dimensionality of the Hilbert space.
nStates = length(Hamiltonian_beta);

U_beta = eye(nStates);
U_alpha = eye(nStates);



doTR = ~isempty(Hamiltonian_beta_TR);
if doTR
  
  % Ensure subHamiltinians are Hermitian.
  Hamiltonian_beta_TR =(Hamiltonian_beta_TR+Hamiltonian_beta_TR')/2;
  Hamiltonian_alpha_TR =(Hamiltonian_alpha_TR+Hamiltonian_alpha_TR')/2;
  
  % Generate propagators for time element dt.
  [dU_beta_TR, dU_beta2_TR] = propagator_eig(Hamiltonian_beta_TR,dt,dt2);
  [dU_alpha_TR, dU_alpha2_TR] = propagator_eig(Hamiltonian_alpha_TR,dt, dt2);
end


% Get second experimental dimesion propagators.
if EXPERIMENT==CPMG
    U_beta_2 = eye(nStates);
    U_alpha_2 = eye(nStates);
end

% Get second experimental dimesion propagators.
if EXPERIMENT == CPMG_CONST
  U_beta_2 = propagator_eig(Hamiltonian_beta,total_time);
  U_alpha_2 = propagator_eig(Hamiltonian_alpha,total_time);
end

% Initialize signal.
v= ones(1 ,timepoints);

% Loop over time points.
for iTime = 1:timepoints
  
  % Find the correct experiment
  switch EXPERIMENT
    
    case FID
      % Generate time dependent detection operator.
      U_ = U_beta'*U_alpha;
      v(iTime) = vecDensityMatrixT*U_(:);
      
    case HAHN
      % Generate time dependent detection operator.
      U_ = U_beta'*U_alpha'*U_beta*U_alpha;
      v(iTime) = vecDensityMatrixT*U_(:);
      
    case CPMG
      % Generate time dependent detection operator.
      
      U_ = U_alpha'*U_beta'  * ...
        (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha*U_beta;
      v(iTime) = vecDensityMatrixT*U_(:);
      
      
    case CPMG_CONST
      
      % THIS NEEDS TO BE UPDATED TO USE dt2 WHEN iTime > Ndt.
      [U_beta, U_beta_2] = propagator_eig(...
        Hamiltonian_beta,(iTime-1)*dt, total_time/4-(iTime-1)*dt);
      [U_alpha, U_alpha_2] = propagator_eig(...
        Hamiltonian_alpha,(iTime-1)*dt, total_time/4-(iTime-1)*dt);
      
      % Generate time dependent detection operator.
      U_ = U_alpha_2'*U_beta_2'  *  ...
        (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
      
      v(iTime) = vecDensityMatrixT*U_(:);
      
      
    case CPMG_2D
      
      U_beta_2 = eye(nStates);
      U_alpha_2 = eye(nStates);
      
      % Loop over second experimental dimension.
      for jTime = 1:iTime %timepoints
        
        % Generate time dependent detection operator.
        U_ = U_alpha_2'*U_beta_2'  *  ...
          (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
        
        v(iTime,jTime) = vecDensityMatrixT*U_(:);
        v(jTime,iTime) = v(iTime,jTime)';
        
        % Increment propagator.
        if jTime < Ndt
          U_beta_2 = dU_beta*U_beta_2;
          U_alpha_2 = dU_alpha*U_alpha_2;
        else
          U_beta_2 = dU_beta2*U_beta_2;
          U_alpha_2 = dU_alpha2*U_alpha_2;
        end
        
      end
    
    case HAHN_TR
      % Generate time dependent detection operator.
      U_ = U_beta'*U_beta_TR'*U_alpha'*U_alpha_TR' ...
        * U_beta_TR*U_beta*U_alpha_TR*U_alpha;
      
      v(iTime) = vecDensityMatrixT*U_(:);
      
      
  end
  
  
  % Increment propagator.
  if iTime< Ndt
    U_beta = dU_beta*U_beta;
    U_alpha = dU_alpha*U_alpha;
  else
    U_beta = dU_beta2*U_beta;
    U_alpha = dU_alpha2*U_alpha;
  end
  
  if doTR
    % Increment propagator.
    if iTime< Ndt
      U_beta_TR = dU_beta_TR*U_beta_TR;
      U_alpha_TR = dU_alpha_TR*U_alpha_TR;
    else
      U_beta_TR = dU_beta2_TR*U_beta_TR;
      U_alpha_TR = dU_alpha2_TR*U_alpha_TR;
    end
    
  end
end

% Return signal as a vector.
if EXPERIMENT == CPMG_2D
    Signal = reshape(v.',1,[]);
else
    Signal = v;
end

if any(abs(v) - 1 > 1e-9)
  error('Coherence error: coherence cannot be larger than 1.');
end 
end


% ========================================================================
% Subcluster Function
% ========================================================================

function Indices = findSubclusters(Clusters,clusterSize,iCluster,...
  reference_clusterSize)
% Given the ith cluster of size clusterSize, 
% Indices(jCluster,subCluster_size) = jth cluster of sizesubCluster_size
% that is a subcluster of the ith cluster of size clusterSize
%
% SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
% the jth cluster of size subCluster_size that is a subcluster of
% the ith ccluster of size clusterSize.
% In a valid cluster, all indices must be natural numbers.


Indices = zeros(nchoosek(reference_clusterSize, ...
  ceil(reference_clusterSize/2)),reference_clusterSize);
jCluster = 0;

% Each cluster is a subset of itself.
Indices(1,clusterSize)=iCluster;

if clusterSize==1
  % All non-empty subsets have been found.
  return;
end

% Loop over all cluster sizes up to clusterSize.
for index = 1:clusterSize
  % CluserArray(iCluster,:,clusterSize) = nuclear indices.
  Cluster = Clusters(iCluster, 1:clusterSize ,clusterSize);
%   Cluster(Cluster==0) = [];
  
  % Remove one element labeled by index from 
  % Clusters{clusterSize}(iCluster:) to get a sub-cluster with one less element.  
  SubCluster = [Cluster(1:index-1), Cluster(index+1:end)];
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:clusterSize -1, clusterSize -1)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2) );
  
  % Skip if there is no subclustere.
  if isempty(subclusterIndex)
    continue;
  end
  
  jCluster = jCluster + 1;
  Indices(jCluster,clusterSize -1) = subclusterIndex;
  
  if clusterSize > 2
    % Given the ith cluster of size clusterSize,
    % Indices(jCluster,subCluster_size) = jth cluster of sizesubCluster_size
    % that is a subcluster of the ith cluster of size clusterSize
    
    % Recursively find smaller subclusters
    Subindices = findSubclusters(Clusters,clusterSize -1, ...
      subclusterIndex,reference_clusterSize);
    
    % Assign subcluster indices. 
    Indices(:,1:clusterSize-1) = Subindices(:,1:clusterSize-1);
    
  end
  
end

end

% ========================================================================
% Calculate propagator using diagonalization
% ========================================================================

function [U1,U2] = propagator_eig(Ham,t1,t2)
%Ham = (Ham+Ham')/2; % "hermitianize" Hamiltonian
[EigenVectors, EigenValues] = eig(Ham);
Udiag1 = exp(-2i*pi*diag(EigenValues)*t1);
Udiag2 = exp(-2i*pi*diag(EigenValues)*t2);
U1 = EigenVectors*diag(Udiag1)*EigenVectors';
U2 = EigenVectors*diag(Udiag2)*EigenVectors';
end


% ========================================================================
% New Function TEMPORARY
% ========================================================================

function Indices = findSubclusters0(Clusters,clusterSize,iCluster)
% Indices{size} = list of all jCluster such that Clusters{size}(jCluster,:) 
% is a subcluster of Clusters{clusterSize}(iCluster,:)

% In a valid cluster, all indices must be natural numbers.
if Clusters{clusterSize}(iCluster,1) == 0
  Indices = [];
  return;
end

% Each cluster is a subset of itself.
Indices{clusterSize}=iCluster;

if clusterSize==1
  % All non-empty subsets have been found.
  return;
end

% Loop over all cluster sizes up to clusterSize.
for index = 1:clusterSize
  
  % Remove one element labeled by index from Clusters{clusterSize}(iCluster:) 
  % to get a sub-cluster with one less element.  
  SubCluster = [Clusters{clusterSize}(iCluster,1:index-1), ...
    Clusters{clusterSize}(iCluster,index+1:end)];
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters{clusterSize -1}==SubCluster;
  for ii = 2:size(Search,2)
    Search(:,1) = Search(:,1).*Search(:,ii);
  end
  subclusterIndex = find(Search(:,1)==1);
  
  if isempty(subclusterIndex)
    continue;
  end
  Indices{clusterSize -1} = [Indices{clusterSize -1} , subclusterIndex];
  
  if clusterSize > 2
    Subindices = findSubclusters0(Clusters,clusterSize -1,subclusterIndex);
    for isize = (clusterSize-2):-1:1
      if ~isempty(Subindices{isize})
        Indices{isize} = [Indices{isize},Subindices{isize}];
      end
    end
  end
  
end

for isize = 1:clusterSize
  Indices{isize} = unique(Indices{isize});
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

% 1-Clusters
if N ==1
  cycPerm = [1,2,3;2,3,1;3,1,2];
  return;
end

% 2-Clusters
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

% 3-Clusters
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
      3,1,2, 6,4,5, 7,8,9 ...
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


% ==============================================================================
% Select the Spin Operator
% ==============================================================================
function [SpinOp,SpinXiXjOp] = getSpinOps(thisClusterSize,Nuclei_Spins,...
  Op, SpinXiXjOp, ...
  SpinOperators_HD, SpinOperators_DH,SpinXiXjOperators_HD,SpinXiXjOperators_DH)

% "ENUM"
HD_SPIN_OP = -1; DH_SPIN_OP = -2;


spinOp = Op{spinMultiplicity}{thisClusterSize};
SpinXiXjOp = SpinXiXjOps{thisClusterSize};

if all(Nuclei_Spins == Nuclei_Spins(1))
  Nuclei_Spin = Nuclei_Spins(1);
elseif length(Nuclei_Spins)==2
  if Nuclei_Spins(1)==1/2 && Nuclei_Spins(2)==1
    Nuclei_Spin = HD_SPIN_OP;
  elseif Nuclei_Spins(2)==1/2 && Nuclei_Spins(1)==1
    Nuclei_Spin = DH_SPIN_OP;
  end
end

switch thisClusterSize
  case 1
    switch Nuclei_Spin
      case 1/2
        SpinOp = Spin2Op1;
        SpinXiXjOp = [];
      case 1
        SpinOp = Spin3Op1;
        SpinXiXjOp = SpinXiXjOp_1;
      case 3/2
        SpinOp = Spin4Op1;
        SpinXiXjOp = [];
      otherwise
        error('Cluster contains unrecognized spins.');
    end
    
  case 2
    switch Nuclei_Spin
      case 1/2
        SpinOp = Spin2Op2;
        SpinXiXjOp = [];
      case 1
        SpinOp = Spin3Op2;
        SpinXiXjOp = SpinXiXjOp_2;
      case 3/2
        SpinOp = Spin4Op2;
        SpinXiXjOp = [];
      case HD_SPIN_OP
        SpinOp = SpinOperators_HD;
        SpinXiXjOp = SpinXiXjOperators_HD;
      case DH_SPIN_OP
        SpinOp = SpinOperators_DH;
        SpinXiXjOp = SpinXiXjOperators_DH;
      otherwise
        error('Cluster contains unrecognized spins.');
    end
    
  case 3
    switch Nuclei_Spin
      case 1/2
        SpinOp = Spin2Op3;
        SpinXiXjOp = [];
      case 1
        SpinOp = Spin3Op3;
        SpinXiXjOp = SpinXiXjOp_3;
      case 3/2
        SpinOp = Spin4Op3;
        SpinXiXjOp = [];
      otherwise
        error('Cluster contains unrecognized spins.');
    end
    
  case 4
    switch Nuclei_Spin
      case 1/2
        SpinOp = Spin2Op4;
        SpinXiXjOp = [];
      case 1
        SpinOp = Spin3Op4;
        SpinXiXjOp = SpinXiXjOp_4;
      case 3/2
        SpinOp = Spin4Op4;
        SpinXiXjOp = [];
      otherwise
        error('Cluster contains unrecognized spins.');
    end
    
  case 5
    switch Nuclei_Spin
      case 1/2
        SpinOp = Spin2Op5;
        SpinXiXjOp = [];
      case 1
        SpinOp = Spin3Op5;
        SpinXiXjOp = SpinXiXjOp_5;
      case 3/2
        SpinOp = Spin4Op5;
        SpinXiXjOp = [];
      otherwise
        error('Cluster contains unrecognized spins.');
    end
    
  case 6
    switch Nuclei_Spin
      case 1/2
        SpinOp = Spin2Op6;
        SpinXiXjOp = [];
      case 1
        SpinOp = Spin3Op6;
        SpinXiXjOp = SpinXiXjOp_6;
      case 3/2
        SpinOp = Spin4Op6;
        SpinXiXjOp = [];
      otherwise
        error('Cluster contains unrecognized spins.');
    end
    
  otherwise
    error('Cluster size is too large: spin operators cannot be assigned.');
    
end
end

% ==============================================================================
% Unpack Spin Operators
% ==============================================================================

function [Spin2Op1,Spin3Op1,Spin4Op1, ...
  Spin2Op2,Spin3Op2,Spin4Op2, ...
  Spin2Op3,Spin3Op3,Spin4Op3, ...
  Spin2Op4,Spin3Op4,Spin4Op4, ...
  Spin2Op5,Spin3Op5,Spin4Op5, ...
  Spin2Op6,Spin3Op6,Spin4Op6, ...
  SpinXiXjOp_1, SpinXiXjOp_2, SpinXiXjOp_3, ...
  SpinXiXjOp_4, SpinXiXjOp_5, SpinXiXjOp_6]...
  = initializeSpinOps(maxClusterSize, Op, SpinXiXjOps)

% Set 1-cluster operators.
if maxClusterSize > 0
  Spin2Op1 = Op{2}{1};
  Spin3Op1 = Op{3}{1};
  Spin4Op1 = Op{4}{1};
  SpinXiXjOp_1 = SpinXiXjOps{1};
else
  Spin2Op1 = [];
  Spin3Op1 = [];
  Spin4Op1 = [];
  SpinXiXjOp_1 = [];
end

% Set 2-cluster operators.
if maxClusterSize > 1
  Spin2Op2 = Op{2}{2};
  Spin3Op2 = Op{3}{2};
  Spin4Op2 = Op{4}{2};
  SpinXiXjOp_2 = SpinXiXjOps{2};
else
  Spin2Op2 = [];
  Spin3Op2 = [];
  Spin4Op2 = [];
  SpinXiXjOp_2 =[];
end

% Set 3-cluster operators.
if maxClusterSize > 2
  Spin2Op3 = Op{2}{3};
  Spin3Op3 = Op{3}{3};
  Spin4Op3 = Op{4}{3};
  SpinXiXjOp_3 = SpinXiXjOps{3};
else
  Spin2Op3 = [];
  Spin3Op3 = [];
  Spin4Op3 = [];
  SpinXiXjOp_3 =[];
end

% Set 4-cluster operators.
if maxClusterSize > 3
  Spin3Op4 = Op{3}{4};
  Spin2Op4 = Op{2}{4};
  Spin4Op4 = Op{4}{4};
  SpinXiXjOp_4 = SpinXiXjOps{4};
else
  Spin3Op4 = [];
  Spin2Op4 = [];
  Spin4Op4 = [];
  SpinXiXjOp_4 =[];
end

% Set 5-cluster operators.
if maxClusterSize > 4
  Spin2Op5 = Op{2}{5};
  Spin3Op5 = Op{3}{5};
  Spin4Op5 = Op{4}{5};
  SpinXiXjOp_5 = SpinXiXjOps{5};
else
  Spin2Op5 = [];
  Spin3Op5 = [];
  Spin4Op5 = [];
  SpinXiXjOp_5 =[];
end

% Set 6-cluster operators.
if maxClusterSize > 5
  Spin2Op6 = Op{2}{6};
  Spin3Op6 = Op{3}{6};
  Spin4Op6 = Op{4}{6};
  SpinXiXjOp_6 = SpinXiXjOps{6};
else
  Spin2Op6 = [];
  Spin3Op6 = [];
  Spin4Op6 = [];
  SpinXiXjOp_6 =[];
end

% Set 7-cluster operators.
if maxClusterSize > 6
  Spin2Op7 = Op{2}{7};
  Spin3Op7 = Op{3}{7};
  Spin4Op7 = Op{4}{7};
  SpinXiXjOp_7 = SpinXiXjOps{7};
else
  Spin2Op7 = [];
  Spin3Op7 = [];
  Spin4Op7 = [];
  SpinXiXjOp_7 =[];
end

% Set 8-cluster operators.
if maxClusterSize > 7
  Spin2Op8 = Op{2}{8};
  Spin3Op8 = Op{3}{8};
  Spin4Op8 = Op{4}{8};
  SpinXiXjOp_8 = SpinXiXjOps{8};
else
  Spin2Op8 = [];
  Spin3Op8 = [];
  Spin4Op8 = [];
  SpinXiXjOp_8 =[];
end

% Set 9-cluster operators.
if maxClusterSize > 8
  Spin2Op9 = Op{2}{9};
  Spin3Op9 = Op{3}{9};
  Spin4Op9 = Op{4}{9};
  SpinXiXjOp_9 = SpinXiXjOps{9};
else
  Spin2Op9 = [];
  Spin3Op9 = [];
  Spin4Op9 = [];
  SpinXiXjOp_9 =[];
end

% Set 9-cluster operators.
if maxClusterSize > 9
  Spin2Op10 = Op{2}{10};
  Spin3Op10 = Op{3}{10};
  Spin4Op10 = Op{4}{10};
  SpinXiXjOp_10 = SpinXiXjOps{10};
else
  Spin2Op10 = [];
  Spin3Op10 = [];
  Spin4Op10 = [];
  SpinXiXjOp_10 =[];
end

% Set 11-cluster operators.
if maxClusterSize > 11
  Spin2Op11 = Op{2}{11};
  Spin3Op11 = Op{3}{11};
  Spin4Op11 = Op{4}{11};
  SpinXiXjOp_11 = SpinXiXjOps{11};
else
  Spin2Op11 = [];
  Spin3Op11 = [];
  Spin4Op11 = [];
  SpinXiXjOp_11 =[];
end

% Set 12-cluster operators.
if maxClusterSize > 12
  Spin2Op12 = Op{2}{12};
  Spin3Op12 = Op{3}{12};
  Spin4Op12 = Op{4}{12};
  SpinXiXjOp_12 = SpinXiXjOps{12};
else
  Spin2Op12 = [];
  Spin3Op12 = [];
  Spin4Op12 = [];
  SpinXiXjOp_12 =[];
end
end


% ========================================================================
% Initialize Coherences
% ========================================================================

function [Coherences, SubclusterIndices , ...
  rotationalMatrix_c1_m1, rotationalMatrix_d1_m1, ...
  rotationalMatrix_c2_m1, rotationalMatrix_c2_m2, ...
  rotationalMatrix_d2_m1, rotationalMatrix_d2_m2] ...
  = initializeCoherences(Method, numberClusters, ...
  Nuclei_rotationalMatrix, nt,  includeMethyls,lockRotor, ...
  Methyl_methylMethylCoupling)

% nt = timepoints^dimensionality;
% includeMethyls = System.Methyl.include
SubclusterIndices = cell(1,Method.order);

for iclusterSize = 1:Method.order
  
  Coherences{iclusterSize} = zeros(numberClusters(iclusterSize),nt);
    
  % SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
  % the jth cluster of size subCluster_size that is a subcluster of
  % the ith ccluster of size clusterSize.
  
  subCluster_size_= ceil(iclusterSize/2);
  
  SubclusterIndices{iclusterSize} = zeros(nchoosek(iclusterSize,subCluster_size_), iclusterSize , ...
    numberClusters(iclusterSize));
  
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

% ========================================================================
% Set Subcluster Indices
% ========================================================================
function [SubclusterIndices_2, SubclusterIndices_3,SubclusterIndices_4,...
  SubclusterIndices_5,SubclusterIndices_6, ...
  SubclusterIndices_7,SubclusterIndices_8,...
  SubclusterIndices_9,SubclusterIndices_10,SubclusterIndices_11, ...
  SubclusterIndices_12]...
  = setSubclusterIndices(ClusterArray,clusterSize,iCluster,...
  SubclusterIndices_2, SubclusterIndices_3,SubclusterIndices_4,...
  SubclusterIndices_5,SubclusterIndices_6, ...
  SubclusterIndices_7,SubclusterIndices_8,...
  SubclusterIndices_9,SubclusterIndices_10,SubclusterIndices_11, ...
  SubclusterIndices_12)


switch clusterSize
  % SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
  % the jth cluster of size subCluster_size that is a subcluster of
  % the ith cluster of size clusterSize.
  case 1
    % This is just to catch this case.  There is nothing to set.
  case 2
    SubclusterIndices_2(:,:,iCluster) = findSubclusters_gpu(...
      ClusterArray,clusterSize,iCluster,clusterSize);
  case 3
    SubclusterIndices_3(:,:,iCluster) = findSubclusters_gpu(...
      ClusterArray,clusterSize,iCluster,clusterSize);
  case 4
    SubclusterIndices_4(:,:,iCluster) = findSubclusters_gpu(...
      ClusterArray,clusterSize,iCluster,clusterSize);
  case 5
    SubclusterIndices_5(:,:,iCluster) = findSubclusters_gpu(...
      ClusterArray,clusterSize,iCluster,clusterSize);
  case 6
    SubclusterIndices_6(:,:,iCluster) = findSubclusters_gpu(...
      ClusterArray,clusterSize,iCluster,clusterSize);
  case 7
    SubclusterIndices_7(:,:,iCluster) = findSubclusters_gpu(...
      ClusterArray,clusterSize,iCluster,clusterSize);
  case 8
    SubclusterIndices_8(:,:,iCluster) = findSubclusters_gpu(...
      ClusterArray,clusterSize,iCluster,clusterSize);
  case 9
    SubclusterIndices_9(:,:,iCluster) = findSubclusters_gpu(...
      ClusterArray,clusterSize,iCluster,clusterSize);
  case 10
    SubclusterIndices_10(:,:,iCluster) = findSubclusters_gpu(...
      ClusterArray,clusterSize,iCluster,clusterSize);
  case 11
    SubclusterIndices_11(:,:,iCluster) = findSubclusters_gpu(...
      ClusterArray,clusterSize,iCluster,clusterSize);
  case 12
    SubclusterIndices_12(:,:,iCluster) = findSubclusters_gpu(...
      ClusterArray,clusterSize,iCluster,clusterSize);
  otherwise
    fprintf('Cannot calculate clusters of size %d.\n', clusterSize);
    error('Cluster Error');
    
    
end
end
