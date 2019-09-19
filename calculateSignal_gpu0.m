% Clusters = Clusters(cluster index , 1:size ,order)
% Clusters(cluster index , size > order ,order) = 0.

function [Signal, AuxiliarySignal_1,AuxiliarySignal_2,AuxiliarySignal_3,AuxiliarySignal_4,Signals] ... 
       = calculateSignal_gpu0(System, Method, Nuclei,Clusters, timepoints,dt)
% ENUM
FID = 1; HAHN = 2; CPMG = 3; CPMG_CONST = 4; CPMG_2D = 5;

Method_order = Method.order;  
System_full_Sz_Hyperfine = System.full_Sz_Hyperfine;
maxClusterSize = min(4,Method_order);
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
else
  error('The experiment ''%s'' is not supported.',System.experiment);
end
Nuclei_Coordinates = Nuclei.Coordinates;
Nuclei_Abundance = Nuclei.Abundance;
Nuclei_Spin = Nuclei.Spin;
Nuclei_g = Nuclei.Nuclear_g;
NumberStates = Nuclei.NumberStates;
ZeemanStates = Nuclei.ZeemanStates;
max_basis = max(NumberStates);
States = zeros(max_basis,Nuclei.number);
for ii = 1:Nuclei.number
  States(1:NumberStates(ii),ii) = Nuclei.State{ii};
end
Spin2Op1 = Nuclei.SpinOperators{2}{1};
Spin2Op2 = Nuclei.SpinOperators{2}{2};
Spin2Op3 = Nuclei.SpinOperators{2}{3};
Spin2Op4 = Nuclei.SpinOperators{2}{4};

Spin3Op1 = Nuclei.SpinOperators{3}{1};
Spin3Op2 = Nuclei.SpinOperators{3}{2};
Spin3Op3 = Nuclei.SpinOperators{3}{3};
Spin3Op4 = Nuclei.SpinOperators{3}{4};

numberClusters = Nuclei.numberClusters(1:maxClusterSize);
maxNumberClusters = max(numberClusters(1:maxClusterSize));

% CluserArray(iCluster,:,clusterSize) = nuclear indices.
ClusterArray = zeros(maxNumberClusters,maxClusterSize,maxClusterSize);
 

for isize = 1:Method_order
  
  for ii = 1:Nuclei.numberClusters(isize)
    ClusterArray(ii,1:isize,isize) = Clusters{isize}(ii,:);
  end
  
      switch isize 
        
        % SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
        % the jth cluster of size subCluster_size that is a subcluster of
        % the ith ccluster of size clusterSize.
        
        case 1
          Coherences_1 = ones(numberClusters(isize),timepoints^dimensionality);  
        case 2
          SubclusterIndices_2 = zeros(nchoosek(isize,1), isize , Nuclei.numberClusters(isize)); 
          Coherences_2 = ones(numberClusters(isize),timepoints^dimensionality);
          
        case 3
          SubclusterIndices_3 = zeros(nchoosek(isize,1), isize , Nuclei.numberClusters(isize));
          Coherences_3 = ones(numberClusters(isize),timepoints^dimensionality);
        case 4
          SubclusterIndices_4 = zeros(nchoosek(isize,2), isize , Nuclei.numberClusters(isize));
          Coherences_4 = ones(numberClusters(isize),timepoints^dimensionality);
      end
      
end

% Define placeholder variables.
for isize = Method_order+1:4
  
      switch isize 
        
        case 1
          Coherences_1 = 1;  
        case 2
          SubclusterIndices_2 = []; 
          Coherences_2 = 1;
          
        case 3
          SubclusterIndices_3 = [];
          Coherences_3 = 1;
        case 4
          SubclusterIndices_4 = 1;
          Coherences_4 = 1;
      end
      
end
ge=System.gMatrix(3,3);
magneticField = System.magneticField;
muB = System.muB;
muN = System.muN;
mu0 = System.mu0;
hbar = System.hbar;

Nuclei_ValidPair = Nuclei.ValidPair;
% ENUM
CONNECTED = 0;  COMPLETE = 1;

graphCriterion = CONNECTED;
if strcmp(Nuclei.graphCriterion,'complete')
  graphCriterion = COMPLETE;
end
%--------------------------------------------------------------------------
% gpu code
%--------------------------------------------------------------------------


for clusterSize = 1:Method_order
   
  % Find coherences
  numClusters = numberClusters(clusterSize);
  
  for iCluster = numClusters:-1:1
    Cluster = ClusterArray(iCluster,1:clusterSize,clusterSize); 
   
    if Cluster(1,1) == 0 || ~validateCluster( Cluster  , Nuclei_ValidPair, graphCriterion, clusterSize)
      continue;
    end
    
    switch clusterSize 
      % SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
      % the jth cluster of size subCluster_size that is a subcluster of
      % the ith ccluster of size clusterSize.
      case 2
        SubclusterIndices_2(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 3
        SubclusterIndices_3(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 4
        SubclusterIndices_4(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
    end
    
    % Select the appropriate spin operator.
     switch clusterSize
      case 1
        switch Nuclei_Spin(Cluster(1))
          case 1/2
            SpinOp = Spin2Op1;
          case 1
            SpinOp = Spin3Op1;
        end
        
       case 2
        switch Nuclei_Spin(Cluster(1))
          case 1/2
            SpinOp = Spin2Op2;
          case 1
            SpinOp = Spin3Op2;
        end
        
       case 3
        switch Nuclei_Spin(Cluster(1))
          case 1/2
            SpinOp = Spin2Op3;
          case 1
            SpinOp = Spin3Op3;
        end
        
       case 4
        switch Nuclei_Spin(Cluster(1))
          case 1/2
            SpinOp = Spin2Op4;
          case 1
            SpinOp = Spin3Op4;
        end
         
     end
     
    [Hamiltonian,zeroIndex] = pairwiseHamiltonian_gpu(Nuclei_g, Nuclei_Coordinates,Cluster,magneticField, ge, muB, muN, mu0, hbar);
    
    ms = -1/2;
    Hb = assembleHamiltonian_gpu(Hamiltonian,SpinOp,Cluster,NumberStates,System_full_Sz_Hyperfine, ms,zeroIndex,clusterSize);
    ms= 1/2;
    Ha = assembleHamiltonian_gpu(Hamiltonian,SpinOp,Cluster,NumberStates,System_full_Sz_Hyperfine, ms,zeroIndex,clusterSize);
  
    % get density matrix
    DensityMatrix = getDensityMatrix(ZeemanStates,NumberStates,Cluster);
    
    % get cluster coherence
    switch clusterSize
      case 1
        Coherences_1(iCluster,:) = propagate(total_time, DensityMatrix,timepoints,dt,Hb,Ha,EXPERIMENT); 
      case 2
        Coherences_2(iCluster,:) = propagate(total_time, DensityMatrix,timepoints,dt,Hb,Ha,EXPERIMENT); 
      case 3
        Coherences_3(iCluster,:) = propagate(total_time, DensityMatrix,timepoints,dt,Hb,Ha,EXPERIMENT); 
        if clusterSize==3 && all(sort(Cluster)==[6,7,11])
          1+1;
        end
      case 4
        Coherences_4(iCluster,:) = propagate(total_time, DensityMatrix,timepoints,dt,Hb,Ha,EXPERIMENT); 

    end
    
    
  end
  
  
end


% Calculate signal
%-------------------------------------------------------------------------------
[Signals, AuxiliarySignal_1,AuxiliarySignal_2,AuxiliarySignal_3,AuxiliarySignal_4] ...
  = doClusterCorrelationExpansion_gpu(Coherences_1,Coherences_2,Coherences_3,Coherences_4,ClusterArray, ...
  SubclusterIndices_2,SubclusterIndices_3,SubclusterIndices_4,...
  timepoints,dimensionality, Method_order,numberClusters, Nuclei_Abundance);
% 
% if EXPERIMENT == CPMG_2D
%     Signal = reshape(Signals(Method_order,:).',timepoints,timepoints)';
% else
%     Signal = Signals(Method_order,:);
% end

    Signal = Signals(Method_order,:);
end

% ========================================================================
% Density Matrix Function
% ========================================================================

function DensityMatrix = getDensityMatrix(States,NumberStates, Cluster)

dimension = prod(NumberStates(Cluster));

DensityMatrix = eye(dimension);

switchIndex = dimension;
for inucleus = Cluster
  switchIndex = switchIndex/NumberStates(inucleus);
  counter = 0;
  stateNumber = 1;
  for jj = 1:dimension
    if counter == switchIndex
      stateNumber = mod(stateNumber,NumberStates(inucleus)) + 1;
    end
    counter = mod(counter,switchIndex) + 1;
    DensityMatrix(jj,jj) = DensityMatrix(jj,jj)*States(inucleus,stateNumber)^2;
  end
end
%DensityMatrix = DensityMatrix.^2;
if abs(trace(DensityMatrix)-1)>1e-9
  error('State is not normalized.');
end

end

% ========================================================================
% Propagate Function
% ========================================================================
function Signal = propagate(total_time,DensityMatrix,timepoints,dt,Hamiltonian_beta,Hamiltonian_alpha,EXPERIMENT)

% ENUM
FID = 1; HAHN = 2; CPMG = 3; CPMG_CONST = 4; CPMG_2D = 5;

vecDensityMatrixT = reshape(DensityMatrix.',1,[]);


dU_beta = propagator_eig(Hamiltonian_beta,dt);
dU_alpha = propagator_eig(Hamiltonian_alpha,dt);

nStates = length(Hamiltonian_beta);
U_beta = eye(nStates);
U_alpha = eye(nStates);

if EXPERIMENT==CPMG
  U_beta_2 = eye(nStates);
  U_alpha_2 = eye(nStates);
end

if EXPERIMENT == CPMG_CONST
  U_beta_2 = propagator_eig(Hamiltonian_beta,total_time);
  U_alpha_2 = propagator_eig(Hamiltonian_alpha,total_time);
end


v= ones(1,timepoints);
for iTime = 1:timepoints
  
  switch EXPERIMENT
    case FID
      U_ = U_beta'*U_alpha;
      v(iTime) = vecDensityMatrixT*U_(:);
      
    case HAHN
      U_ = U_beta'*U_alpha'*U_beta*U_alpha;
      v(iTime) = vecDensityMatrixT*U_(:);
      
    case CPMG
      U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
      
      v(iTime) = vecDensityMatrixT*U_(:);
      
      U_beta_2 = dU_beta*U_beta_2;
      U_alpha_2 = dU_alpha*U_alpha_2;
      
    case CPMG_CONST
      
      U_beta = propagator_eig(Hamiltonian_beta,(iTime-1)*dt);
      U_alpha = propagator_eig(Hamiltonian_alpha,(iTime-1)*dt);
      
      U_beta_2 = propagator_eig(Hamiltonian_beta, total_time/4-(iTime-1)*dt);
      U_alpha_2 = propagator_eig(Hamiltonian_alpha, total_time/4-(iTime-1)*dt);
      U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
      
      v(iTime) = vecDensityMatrixT*U_(:);
      
      
    case CPMG_2D
      
      U_beta_2 = eye(nStates);
      U_alpha_2 = eye(nStates);
      
      for jTime = 1:timepoints
        
        U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
        
        v(iTime,jTime) = vecDensityMatrixT*U_(:);
        
        U_beta_2 = dU_beta*U_beta_2;
        U_alpha_2 = dU_alpha*U_alpha_2;
      end
      
  end
  
  
  U_beta = dU_beta*U_beta;
  U_alpha = dU_alpha*U_alpha;
  
  
end
 
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

function Indices = findSubclusters(Clusters,clusterSize,iCluster,reference_clusterSize)
% from Indices{subCluster_size} = list of all jCluster such that Clusters{subCluster_size}(jCluster,:) is a subcluster of Clusters{clusterSize}(iCluster,:)
% to
% Given the ith cluster of size clusterSize, 
% Indices(jCluster,subCluster_size) = jth cluster of sizesubCluster_size
% that is a subcluster of the ith cluster of size clusterSize
%
% SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
% the jth cluster of size subCluster_size that is a subcluster of
% the ith ccluster of size clusterSize.
% In a valid cluster, all indices must be natural numbers.


Indices = zeros(nchoosek(reference_clusterSize, ceil(reference_clusterSize/2)),reference_clusterSize);
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
  
  % Remove one element labeled by index from Clusters{clusterSize}(iCluster:) to get a sub-cluster with one less element.  
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
    Subindices = findSubclusters(Clusters,clusterSize -1,subclusterIndex,reference_clusterSize);
    
    % Assign subcluster indices. 
    Indices(:,1:clusterSize-1) = Subindices(:,1:clusterSize-1);
    
  end
  
end

end


% ========================================================================
% Subcluster Function gpu
% ========================================================================

function Indices = findSubclusters_gpu(Clusters,clusterSize,iCluster,reference_clusterSize)
% from Indices{subCluster_size} = list of all jCluster such that Clusters{subCluster_size}(jCluster,:) is a subcluster of Clusters{clusterSize}(iCluster,:)
% to
% Given the ith cluster of size clusterSize, 
% Indices(jCluster,subCluster_size) = jth cluster of sizesubCluster_size
% that is a subcluster of the ith cluster of size clusterSize
%
% SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
% the jth cluster of size subCluster_size that is a subcluster of
% the ith ccluster of size clusterSize.
% In a valid cluster, all indices must be natural numbers.


Indices = zeros(nchoosek(reference_clusterSize, ceil(reference_clusterSize/2)),reference_clusterSize);
jCluster = 0;

% Each cluster is a subset of itself.
Indices(1,clusterSize)=iCluster;

switch clusterSize
  case 1
    % All non-empty subsets have been found.
    return;
    
  case 3
    
    numberSubClusters = [3,3,1];
   
    possibleSubClusters_2 = [0,1,1; 1, 0,1; 1,1,0];
    
  case 4
    
    numberSubClusters = [4,6,4,1];
    
    possibleSubClusters_2 = [0,0,1,1; 0,1,0,1; 1,0,0,1; 0,1,1,0; 1,0,1,0; 1,1,0,0];
    possibleSubClusters_3 = [0,1,1,1; 1,0,1,1; 1,1,0,1; 1,1,1,0];
end  

Cluster = Clusters(iCluster, 1:clusterSize ,clusterSize);
Indices(1:clusterSize,1) = Cluster;

% for jCluster = 1:numberSubClusters(1)
%   
%   SubCluster = Cluster(possibleSubClusters_1==1);
%   
%   % Search for a valid cluster that equals the subcluster.
%   Search = Clusters(:, 1:clusterSize -1, clusterSize -1)==SubCluster;
%  
%   % Locate where the subcluster is.
%   subclusterIndex = find( all(Search,2));
%   
%   Indices(jCluster,1) = subclusterIndex;
% end

if clusterSize == 2
  return;
end
  
for jCluster = 1:numberSubClusters(2)
  
  SubCluster = Cluster(possibleSubClusters_2(jCluster,:)==1);
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:2, 2)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2));
 
  if isempty(subclusterIndex)
    subclusterIndex = 0;
  end
  
  Indices(jCluster,2) = subclusterIndex;
end

if clusterSize == 3
  return;
end
  
for jCluster = 1:numberSubClusters(3)
  
  SubCluster = Cluster(possibleSubClusters_3(jCluster,:)==1);
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:3, 3)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2));
 
  if isempty(subclusterIndex)
    subclusterIndex = 0;
  end
  
  Indices(jCluster,3) = subclusterIndex;
end


end
% ========================================================================
% Calculate propagator using diagonalization
% ========================================================================

function U = propagator_eig(Ham,t)
%Ham = (Ham+Ham')/2; % "hermitianize" Hamiltonian
[EigenVectors, EigenValues] = eig(Ham);
Udiag = exp(-2i*pi*diag(EigenValues)*t);
U = EigenVectors*diag(Udiag)*EigenVectors';
end


% ========================================================================
% New Function TEMPORARY
% ========================================================================

function Indices = findSubclusters0(Clusters,clusterSize,iCluster)
% Indices{size} = list of all jCluster such that Clusters{size}(jCluster,:) is a subcluster of Clusters{clusterSize}(iCluster,:)

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
  
  % Remove one element labeled by index from Clusters{clusterSize}(iCluster:) to get a sub-cluster with one less element.  
  SubCluster = [Clusters{clusterSize}(iCluster,1:index-1),Clusters{clusterSize}(iCluster,index+1:end)];
  
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


% ========================================================================
% New Function
% ========================================================================
function Coherences = propagate0(System, Method, DensityMatrix,timepoints,dt,Hamiltonian_beta,Hamiltonian_alpha,isotopeProbability, linearTimeAxis,verbose)
vecDensityMatrixT = reshape(DensityMatrix.',1,[]);

if strcmp('time-domain',Method.propagationDomain)
  
if linearTimeAxis
  dU_beta = propagator_eig(Hamiltonian_beta,dt);
  dU_alpha = propagator_eig(Hamiltonian_alpha,dt);
    
  nStates = length(Hamiltonian_beta);
  U_beta = eye(nStates);
  U_alpha = eye(nStates);
  
  if strcmp(System.experiment,'CPMG')
    U_beta_2 = eye(nStates);
    U_alpha_2 = eye(nStates);
  end
  
  if strcmp(System.experiment,'CPMG-const')
    U_beta_2 = propagator_eig(Hamiltonian_beta,System.total_time);
    U_alpha_2 = propagator_eig(Hamiltonian_alpha,System.total_time);
  end
end

Signal_= ones(1,timepoints);
for iTime = 1:timepoints
  if ~linearTimeAxis
    t = System.Time(iTime);
    U_beta = propagator_eig(Hamiltonian_beta,t);
    U_alpha = propagator_eig(Hamiltonian_alpha,t);
  end
  
  if strcmp(System.experiment,'FID')
    U_ = U_beta'*U_alpha;
    Signal_(iTime) = vecDensityMatrixT*U_(:);
    
  elseif strcmp(System.experiment,'Hahn')
    U_ = U_beta'*U_alpha'*U_beta*U_alpha;
    Signal_(iTime) = vecDensityMatrixT*U_(:);
  elseif strcmp(System.experiment,'CPMG')
    
    U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
    
    Signal_(iTime) = vecDensityMatrixT*U_(:);
    
    U_beta_2 = dU_beta*U_beta_2;
    U_alpha_2 = dU_alpha*U_alpha_2;
    
  elseif strcmp(System.experiment,'CPMG-const')
    
    U_beta = propagator_eig(Hamiltonian_beta,(iTime-1)*dt);
    U_alpha = propagator_eig(Hamiltonian_alpha,(iTime-1)*dt);

    U_beta_2 = propagator_eig(Hamiltonian_beta,System.total_time/4-(iTime-1)*dt);
    U_alpha_2 = propagator_eig(Hamiltonian_alpha,System.total_time/4-(iTime-1)*dt);
    U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
      
    Signal_(iTime) = vecDensityMatrixT*U_(:);
    
%     U_beta_2 = dU_beta'*U_beta_2;
%     U_alpha_2 = dU_alpha'*U_alpha_2;
    
  elseif strcmp(System.experiment,'CPMG-2D')
    
    U_beta_2 = eye(nStates);
    U_alpha_2 = eye(nStates);
    
    for jTime = 1:timepoints
    
      U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;  
      
      Signal_(iTime,jTime) = vecDensityMatrixT*U_(:);
      
      U_beta_2 = dU_beta*U_beta_2;
      U_alpha_2 = dU_alpha*U_alpha_2;
    end
   
  else
    error('The experiment ''%s'' is not supported.',System.experiment);
  end
  
%   Signal_(iTime) = vecDensityMatrixT*U_(:);
  
  if linearTimeAxis
    U_beta = dU_beta*U_beta;
    U_alpha = dU_alpha*U_alpha;
  end
  
  
end
 
elseif strcmp('frequency-domain',Method.propagationDomain)
 
% Diagonalize sub-Hamiltonians for alpha and beta electron manifolds
[Vec_beta, E_beta] = eig(2*pi*Hamiltonian_beta);
[Vec_alpha, E_alpha] = eig(2*pi*Hamiltonian_alpha);
E_beta = diag(E_beta);
E_alpha = diag(E_alpha);
M = Vec_alpha'*Vec_beta; % Mims matrix = <a|b> overlap matrix
Madj = M';

prefactor = 1; % include orienation weights etc here

G = prefactor*M; % prepared density matrix before first varied evolution time
D = M; % detection matrix after last varied evolution time
T1left = Madj; % left transfer matrix due to pi pulse
T1right = Madj; % right transfer matrix due to pi pulse

x=2^20;
nPoints = x;

IncSchemeID = 2; % for 2p ESEEM
buffRe = zeros(1,nPoints); % real part of spectral histogram
buffIm = zeros(1,nPoints); % imag part of spectral histogram
dt = System.dt/2/pi; % time step, in microseconds

sf_peaks(IncSchemeID,buffRe,buffIm,dt,[1 2],[2 1],E_alpha,E_beta,G,D,T1left,T1right);

Signal_ =ifft(buffRe+1i*buffIm);

Signal_ = Signal_(1:timepoints);
Signal_=Signal_./Signal_(1);

elseif strcmp('frequency-domain --matlab',Method.propagationDomain)
  
  nStates = length(Hamiltonian_beta);
  
  [Vec_beta, E_beta] = eig(2*pi*Hamiltonian_beta);
  [Vec_alpha, E_alpha] = eig(2*pi*Hamiltonian_alpha);
  E_beta = diag(E_beta);
  E_alpha = diag(E_alpha);
  dE_beta = E_beta - E_beta.';
  dE_alpha = E_alpha - E_alpha.';
  M = Vec_alpha'*Vec_beta; % Mims matrix
  Madj = M';
  
  % frequency domain variables
  n_omega = 2^20 - 1;

  domega = 2*pi/dt/(n_omega-1);
  
  
  omega_amp = zeros(1,n_omega);
  
  for ii = 1:nStates
    for jj = 1:nStates
      for kk = 1:nStates
        for ll = 1:nStates

          amplitude = Madj(ii,jj) * M(jj,kk) * Madj(kk,ll) * M(ll,ii)/2;
          omega_ = -( dE_beta(ii,kk) + dE_alpha(jj,ll));
          idx = mod(round(omega_./domega),n_omega) + 1;
          
          omega_amp(idx) = omega_amp(idx) + amplitude;
          
        end
      end
    end
  end

  Signal_ = ifft(omega_amp,n_omega);
  Signal_ = Signal_(1:timepoints);
  Signal_ =Signal_ /Signal_(1);
end
if Method.gpu && gpuDeviceCount > 0
  Signal_ = gather(Signal_);
end

if any(abs(Signal_) - 1 > 1e-9)
  error('Coherence error: coherence cannot be larger than 1.');
end
% Coherences = 1 + isotopeProbability*(Signal_ - 1);
Coherences = Signal_;
end
