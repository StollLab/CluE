% Clusters = Clusters(cluster index , 1:size ,order)
% Clusters(cluster index , size > order ,order) = 0.

function [Signal, AuxiliarySignal_1,AuxiliarySignal_2,AuxiliarySignal_3,AuxiliarySignal_4,Signals] ... 
       = calculateSignal_gpu(System, Method, Nuclei,Clusters)

ge=System.ge;
geff=System.gMatrix(3,3);
magneticField = System.magneticField;
muB = System.muB;
muN = System.muN;
mu0 = System.mu0;
hbar = System.hbar;
timepoints = System.timepoints;
dt = System.dt;
t0 = System.t0;
dt2 = System.dt2;
Ndt = System.Ndt;
maxSize = 6;     
     
% ENUM
FID = 1; HAHN = 2; CPMG = 3; CPMG_CONST = 4; CPMG_2D = 5;

Method_order = Method.order;


% Define the theory to use at the given cluster size.
Theory = System.Theory;
Method_useMultipleBathStates = Method.useMultipleBathStates;
Method_useInterlacedClusters = Method.useInterlacedClusters;

theory = Theory(Method_order,:);
useMeanField = theory(10);
useMultipleBathStates = Method_useMultipleBathStates & useMeanField;
useInterlacedClusters = Method_useInterlacedClusters & useMeanField;

maxClusterSize = min(maxSize,Method_order);

if Theory(Method_order,10)
  Method_extraOrder = Method.extraOrder;
  maxSuperclusterSize = Method_extraOrder;
else
  Method_extraOrder = Method_order;
  maxSuperclusterSize = maxClusterSize;
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
else
  error('The experiment ''%s'' is not supported.',System.experiment);
end
Nuclei_Coordinates = Nuclei.Coordinates;
Nuclei_Abundance = Nuclei.Abundance;
Nuclei_Spin = Nuclei.Spin;
Nuclei_g = Nuclei.Nuclear_g;
NumberStates = Nuclei.NumberStates;
state_multiplicity = Nuclei.StateMultiplicity;
ZeemanStates = Nuclei.ZeemanStates;
MeanFieldCoefficients = Nuclei.MeanFieldCoefficients;
MeanFieldTotal = Nuclei.MeanFieldTotal;
max_basis = max(NumberStates);
useThermalEnsemble = System.useThermalEnsemble;
States = zeros(max_basis,Nuclei.number);
betaT = 2*pi*System.hbar /System.kT; % 1/Hz.
Qtensors = Nuclei.Qtensor;
Atensors = Nuclei.Atensor;


FermiContact = Nuclei.FermiContact;
nStates0 =System.nStates;
nStates =System.nStates;
% Unpackspin operators.
Op = Nuclei.SpinOperators;
SpinXiXjOps = Nuclei.SpinXiXjOperators;

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

numberClusters = Nuclei.numberClusters(1:maxSuperclusterSize);
maxNumberClusters = max(numberClusters(1:maxSuperclusterSize));


% CluserArray(iCluster,:,clusterSize) = nuclear indices.
ClusterArray = zeros(maxNumberClusters,maxClusterSize,maxClusterSize);
 

for isize = 1:Method_extraOrder
  
  for ii = 1:Nuclei.numberClusters(isize)
    ClusterArray(ii,1:isize,isize) = Clusters{isize}(ii,:);
  end
end

for isize = 1:Method_order
      switch isize 
        
        % SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
        % the jth cluster of size subCluster_size that is a subcluster of
        % the ith ccluster of size clusterSize.
        
        case 1
          Coherences_1 = zeros(numberClusters(isize),timepoints^dimensionality);  
        case 2
          SubclusterIndices_2 = zeros(nchoosek(isize,1), isize , Nuclei.numberClusters(isize)); 
          Coherences_2 = zeros(numberClusters(isize),timepoints^dimensionality);
          
        case 3
          SubclusterIndices_3 = zeros(nchoosek(isize,1), isize , Nuclei.numberClusters(isize));
          Coherences_3 = zeros(numberClusters(isize),timepoints^dimensionality);
        case 4
          SubclusterIndices_4 = zeros(nchoosek(isize,2), isize , Nuclei.numberClusters(isize));
          Coherences_4 = zeros(numberClusters(isize),timepoints^dimensionality);
        case 5
          SubclusterIndices_5 = zeros(nchoosek(isize,3), isize , Nuclei.numberClusters(isize));
          Coherences_5 = zeros(numberClusters(isize),timepoints^dimensionality);
        case 6
          SubclusterIndices_6 = zeros(nchoosek(isize,3), isize , Nuclei.numberClusters(isize));
          Coherences_6 = zeros(numberClusters(isize),timepoints^dimensionality);
        case 7
          SubclusterIndices_7 = zeros(nchoosek(isize,3), isize , Nuclei.numberClusters(isize));
          Coherences_7 = zeros(numberClusters(isize),timepoints^dimensionality);
        case 8
          SubclusterIndices_8 = zeros(nchoosek(isize,3), isize , Nuclei.numberClusters(isize));
          Coherences_8 = zeros(numberClusters(isize),timepoints^dimensionality);
        case 9
          SubclusterIndices_9 = zeros(nchoosek(isize,3), isize , Nuclei.numberClusters(isize));
          Coherences_9 = zeros(numberClusters(isize),timepoints^dimensionality);
        case 10
          SubclusterIndices_10 = zeros(nchoosek(isize,3), isize , Nuclei.numberClusters(isize));
          Coherences_10 = zeros(numberClusters(isize),timepoints^dimensionality);
        case 11
          SubclusterIndices_11 = zeros(nchoosek(isize,3), isize , Nuclei.numberClusters(isize));
          Coherences_11 = zeros(numberClusters(isize),timepoints^dimensionality);
        case 12
          SubclusterIndices_12 = zeros(nchoosek(isize,3), isize , Nuclei.numberClusters(isize));
          Coherences_12 = zeros(numberClusters(isize),timepoints^dimensionality);
      end
      
end

% Define placeholder variables for variables that need to be defined but do
% not contribute to the calculation for the selected order.
for isize = Method_order+1:maxSize
  
      switch isize 
        
        case 1
          Coherences_1 = 0;  
        case 2
          SubclusterIndices_2 = []; 
          Coherences_2 = 0;
          
        case 3
          SubclusterIndices_3 = [];
          Coherences_3 = 0;
        case 4
          SubclusterIndices_4 = 1;
          Coherences_4 = 0;
          
        case 5
          SubclusterIndices_5 = [];
          Coherences_5 = 0;
        case 6
          SubclusterIndices_6 = 1;
          Coherences_6 = 0;
        case 7
          SubclusterIndices_7 = 1;
          Coherences_7 = 0;
        case 8
          SubclusterIndices_8 = 1;
          Coherences_8 = 0;
        case 9
          SubclusterIndices_9 = 1;
          Coherences_9 = 0;
        case 10
          SubclusterIndices_10 = 1;
          Coherences_10 = 0;
        case 11
          SubclusterIndices_11 = 1;
          Coherences_11 = 0;
        case 12
          SubclusterIndices_12 = 1;
          Coherences_12 = 0;
      end
      
end

Adjacency = Nuclei.Adjacency;
% ENUM
CONNECTED = 0;  COMPLETE = 1;

graphCriterion = CONNECTED;
if strcmp(Nuclei.graphCriterion,'complete')
  graphCriterion = COMPLETE;
end

% Methyl Groups
IsMethyl = strcmp(Nuclei.Type,'CH3');
if ~isempty(Nuclei.Methyl_Data)
  Methyl_P3 = Nuclei.Methyl_Data.Projection_3;
  Methyl_P4 = Nuclei.Methyl_Data.Projection_4;
  Methyl_P5 = Nuclei.Methyl_Data.Projection_5;
  Methyl_P6 = Nuclei.Methyl_Data.Projection_6;
end
Methyl_gA = [0,0];
%--------------------------------------------------------------------------
% gpu code
%--------------------------------------------------------------------------


for clusterSize = 1:Method_order
 
  % Find coherences
  numClusters = numberClusters(clusterSize);
  
  for iCluster = numClusters:-1:1
    Cluster = ClusterArray(iCluster,1:clusterSize,clusterSize); 
   
    if Cluster(1,1) == 0 || ~validateCluster(Cluster,Adjacency,graphCriterion)
      continue;
    end
    
    thisClusterSize = clusterSize + 2*sum(IsMethyl(Cluster));
    thisCluster = zeros(1,thisClusterSize);
    thisIndex = 0;
    MethylID = zeros(1,thisClusterSize);
    methyl_number = 0;
    for ii = 1:clusterSize

      thisIndex = thisIndex + 1;
      
      if IsMethyl(Cluster(ii))
        methyl_number = methyl_number + 1;
        Methyl_gA(methyl_number) = Nuclei_Abundance(Cluster(ii));
        MethylID(thisIndex:thisIndex+2) = methyl_number;
        thisCluster(thisIndex) = Cluster(ii) + 1;
        thisIndex = thisIndex + 1;
        thisCluster(thisIndex) = Cluster(ii) + 2;
        thisIndex = thisIndex + 1;
        thisCluster(thisIndex) = Cluster(ii) + 3;
      else
        thisCluster(thisIndex) = Cluster(ii);
      end
      
    end
    if ~all(state_multiplicity(thisCluster)==state_multiplicity(thisCluster(1)))
      switch clusterSize
        case 1
          Coherences_1(iCluster,:) = ones(size(Coherences_1(iCluster,:) ));
        case 2
          Coherences_2(iCluster,:) = ones(size(Coherences_2(iCluster,:) ));
        case 3
          Coherences_3(iCluster,:) = ones(size(Coherences_3(iCluster,:) ));
        case 4
          Coherences_4(iCluster,:) = ones(size(Coherences_4(iCluster,:) ));
        case 5
          Coherences_5(iCluster,:) = ones(size(Coherences_5(iCluster,:) ));
        case 6
          Coherences_6(iCluster,:) = ones(size(Coherences_6(iCluster,:) ));          
      end
    
      continue;
    end
    switch clusterSize 
      % SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
      % the jth cluster of size subCluster_size that is a subcluster of
      % the ith ccluster of size clusterSize.
      case 1
        % This is just to catch this case.  There is nothing to set. 
      case 2
        SubclusterIndices_2(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);   
      case 3
        SubclusterIndices_3(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 4
        SubclusterIndices_4(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 5
        SubclusterIndices_5(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 6
        SubclusterIndices_6(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 7
        SubclusterIndices_7(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 8
        SubclusterIndices_8(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 9
        SubclusterIndices_9(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 10
        SubclusterIndices_10(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 11
        SubclusterIndices_11(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 12
        SubclusterIndices_12(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      otherwise
        fprintf('Cannot calculate clusters of size %d.\n', clusterSize);
        error('Cluster Error');
        continue;
        
    end
    
    % Select the appropriate spin operator.
     switch thisClusterSize
      case 1
        switch Nuclei_Spin(Cluster(1))
          case 1/2
            SpinOp = Spin2Op1;
            SpinXiXjOp = [];
          case 1
            SpinOp = Spin3Op1;
            SpinXiXjOp = SpinXiXjOp_1;
          case 3/2
            SpinOp = Spin4Op1;
            SpinXiXjOp = [];
        end
        
       case 2
         switch Nuclei_Spin(Cluster(1))
           case 1/2
             SpinOp = Spin2Op2;
             SpinXiXjOp = [];
           case 1
             SpinOp = Spin3Op2;
             SpinXiXjOp = SpinXiXjOp_2;
           case 3/2
             SpinOp = Spin4Op2;
             SpinXiXjOp = [];
         end
         
       case 3
         switch Nuclei_Spin(Cluster(1))
           case 1/2
             SpinOp = Spin2Op3;
             SpinXiXjOp = [];
           case 1
             SpinOp = Spin3Op3;
             SpinXiXjOp = SpinXiXjOp_3;
           case 3/2
             SpinOp = Spin4Op3;
             SpinXiXjOp = [];
         end
         
       case 4
         switch Nuclei_Spin(Cluster(1))
           case 1/2
             SpinOp = Spin2Op4;
             SpinXiXjOp = [];
           case 1
             SpinOp = Spin3Op4;
             SpinXiXjOp = SpinXiXjOp_4;
           case 3/2
             SpinOp = Spin4Op4;
             SpinXiXjOp = [];
         end
         
       case 5
         switch Nuclei_Spin(Cluster(1))
           case 1/2
             SpinOp = Spin2Op5;
             SpinXiXjOp = [];
           case 1
             SpinOp = Spin3Op5;
             SpinXiXjOp = SpinXiXjOp_5;
           case 3/2
             SpinOp = Spin4Op5;
             SpinXiXjOp = [];
         end
         
       case 6
         switch Nuclei_Spin(Cluster(1))
           case 1/2
             SpinOp = Spin2Op6;
             SpinXiXjOp = [];
           case 1
             SpinOp = Spin3Op6;
             SpinXiXjOp = SpinXiXjOp_6;
           case 3/2
             SpinOp = Spin4Op6;
             SpinXiXjOp = [];
         end
         
       otherwise
        continue;
         
     end
     
    [tensors,zeroIndex] = pairwisetensors_gpu(Nuclei_g, Nuclei_Coordinates,thisCluster,Atensors,magneticField, ge,geff, muB, muN, mu0, hbar,theory,MethylID);
    qtensors = Qtensors(:,:,thisCluster);
    


    
    
    if useInterlacedClusters(Method_order,clusterSize)
      
      % Get union of all clusters that contain the primary cluster.
      SuperclusterUnion = findSuperClusterUnion(ClusterArray,clusterSize,iCluster);
      
      % Get the set of nuclei not in the primary cluster that unions with the primary
      % cluster to give SuperclusterUnion.
      ClusterComplement = SuperclusterUnion(  ~any(SuperclusterUnion==thisCluster',1)  );
      
      % Determine the number of spin states required to describe SuperclusterUnion.
      tensors0 = tensors;
        if useThermalEnsemble
          nStates(clusterSize) = max( 1, min(nStates0(clusterSize),prod(state_multiplicity(ClusterComplement))));
        else
          nStates(clusterSize) = max( 1, min(nStates0(clusterSize),prod(state_multiplicity(SuperclusterUnion))));
        end
    else
      [H_alpha,H_beta] = ...
        assembleHamiltonian_gpu(state_multiplicity(thisCluster),tensors,SpinOp,qtensors,SpinXiXjOp,...
        theory,zeroIndex,methyl_number);
    end
    
    H_alphaMF = 0;
    H_betaMF = 0;
    
    for iave = 1:nStates(clusterSize)
      
      if useInterlacedClusters(Method_order,clusterSize)
        
        
        if ~useThermalEnsemble
          % density matrix size
          dm_size = prod(state_multiplicity(thisCluster));
          
          % Initialize density matrix.
          densityMatrix = zeros(dm_size);
          
          % Set cluster spin-state.
          dm_index = mod(iave,dm_size) + 1;
          densityMatrix(dm_index,dm_index) = 1;
          iState = floor(iave/dm_size)+1;
        else
          densityMatrix = [];
          
          iState = iave;
        end
        
        % Adjust tensors for external fields.
        if length(ClusterComplement) > length(thisCluster)
%           tensors = getInterlacedTensors(tensors0,zeroIndex, thisCluster,ClusterComplement,iState,...
%                     Nuclei_g,state_multiplicity,Nuclei_Coordinates, muN, mu0, hbar);
          tensors = tensors0; % TEMPORARY
        else
          tensors = tensors0;
          if iave >1
            continue;
          end
        end
        [H_alpha,H_beta] = ...
          assembleHamiltonian_gpu(state_multiplicity(thisCluster),tensors,SpinOp,qtensors,SpinXiXjOp,...
          theory,zeroIndex,methyl_number);
        
      else
        
        if useMultipleBathStates
          [H_alphaMF,H_betaMF] = assembleMeanFieldHamiltonian_gpu(state_multiplicity(thisCluster),tensors,SpinOp,qtensors,SpinXiXjOp,...
            theory,zeroIndex,methyl_number, MeanFieldCoefficients(:,:,:,iave), MeanFieldTotal(iave));
        end
        
        if useThermalEnsemble
          densityMatrix = [];
        else
          densityMatrix = getDensityMatrix(ZeemanStates(iave,thisCluster),state_multiplicity(thisCluster),thisCluster);
        end
        
      end
      
      Halpha = H_alpha + H_alphaMF;
      Hbeta = H_beta + H_betaMF;
      
      switch methyl_number
        case 0
          PA = eye(size(Halpha));
        case 1
          gA = Methyl_gA(1);
          gE = 1-gA;
          methyl_index = find(MethylID==1,1);
          switch thisClusterSize
            case 3
              PA = Methyl_P3(:,:,1);
              PEap = Methyl_P3(:,:,2);
              PEam = Methyl_P3(:,:,3);
              PEbp = Methyl_P3(:,:,4);
              PEbm = Methyl_P3(:,:,5);
              
            case 4
              switch methyl_index
                case 1
                  jmethyl = 6;
                case 2
                  jmethyl = 1;
              end
              PA = Methyl_P4(:,:,jmethyl);
              PEap = Methyl_P4(:,:,jmethyl+1);
              PEam = Methyl_P4(:,:,jmethyl+2);
              PEbp = Methyl_P4(:,:,jmethyl+3);
              PEbm = Methyl_P4(:,:,jmethyl+4);
              
            case 5
              switch methyl_index
                case 1
                  jmethyl = 11;
                case 2
                  jmethyl = 6;
                case 3
                  jmethyl = 1;
              end
              PA   = Methyl_P5(:,:,jmethyl);
              PEap = Methyl_P5(:,:,jmethyl+1);
              PEam = Methyl_P5(:,:,jmethyl+2);
              PEbp = Methyl_P5(:,:,jmethyl+3);
              PEbm = Methyl_P5(:,:,jmethyl+4);
              
            case 6
              switch methyl_index
                case 1
                  jmethyl = 16;
                case 2
                  jmethyl = 11;
                case 3
                  jmethyl = 6;
                case 4
                  jmethyl = 1;
              end
              PA   = Methyl_P6(:,:,jmethyl);
              PEap = Methyl_P6(:,:,jmethyl+1);
              PEam = Methyl_P6(:,:,jmethyl+2);
              PEbp = Methyl_P6(:,:,jmethyl+3);
              PEbm = Methyl_P6(:,:,jmethyl+4);
              
          end
        case 2
          gA = Methyl_gA(1)*Methyl_gA(2);
          gE = (1-Methyl_gA(1))*(1-Methyl_gA(2));
          
          jmethyl = 21;
          PA1   = Methyl_P6(:,:,jmethyl);
          PEap1 = Methyl_P6(:,:,jmethyl+1);
          PEam1 = Methyl_P6(:,:,jmethyl+2);
          PEbp1 = Methyl_P6(:,:,jmethyl+3);
          PEbm1 = Methyl_P6(:,:,jmethyl+4);
          
          PA2   = Methyl_P6(:,:,jmethyl+5);
          PEap2 = Methyl_P6(:,:,jmethyl+6);
          PEam2 = Methyl_P6(:,:,jmethyl+7);
          PEbp2 = Methyl_P6(:,:,jmethyl+8);
          PEbm2 = Methyl_P6(:,:,jmethyl+9);
          
          PA = PA1*PA2;
          PEap = PEap1*PEap2;
          PEam = PEam1*PEam2;
          PEbp = PEbp1*PEbp2;
          PEbm = PEbm1*PEbm2;
          
      end
      
      % Project Hamiltonian onto the methyls' C3 A state,
      % if there are methyls, otherwise multiply by the identity.
      Hb=PA*Hbeta*PA;
      Ha=PA*Halpha*PA;
   
      % get cluster coherence
      
      switch clusterSize
        case 1
          Coherences_1(iCluster,:) = Coherences_1(iCluster,:)  +  ...
            propagate(total_time,timepoints,dt,dt2,Ndt,Hb,Ha,EXPERIMENT,densityMatrix, useThermalEnsemble, betaT);
        case 2
          Coherences_2(iCluster,:) = Coherences_2(iCluster,:)  +...
            propagate(total_time,timepoints,dt,dt2,Ndt,Hb,Ha,EXPERIMENT,densityMatrix, useThermalEnsemble, betaT);
        case 3
          Coherences_3(iCluster,:) = Coherences_3(iCluster,:) + ...
            propagate(total_time,timepoints,dt,dt2,Ndt,Hb,Ha,EXPERIMENT,densityMatrix, useThermalEnsemble, betaT);
        case 4
          Coherences_4(iCluster,:) = Coherences_4(iCluster,:) + ...
            propagate(total_time,timepoints,dt,dt2,Ndt,Hb,Ha,EXPERIMENT,densityMatrix, useThermalEnsemble, betaT);
        case 5
          Coherences_5(iCluster,:) = Coherences_5(iCluster,:) + ...
            propagate(total_time,timepoints,dt,dt2,Ndt,Hb,Ha,EXPERIMENT,densityMatrix, useThermalEnsemble, betaT);
        case 6
          Coherences_6(iCluster,:) = Coherences_6(iCluster,:) + ...
            propagate(total_time,timepoints,dt,dt2,Ndt,Hb,Ha,EXPERIMENT,densityMatrix, useThermalEnsemble, betaT);
          
      end
    end
    
    switch clusterSize
        case 1
          Coherences_1(iCluster,:) = Coherences_1(iCluster,:)/Coherences_1(iCluster,1);
        case 2
          Coherences_2(iCluster,:) = Coherences_2(iCluster,:)/Coherences_2(iCluster,1);
        case 3
          Coherences_3(iCluster,:) = Coherences_3(iCluster,:)/Coherences_3(iCluster,1);
        case 4
          Coherences_4(iCluster,:) = Coherences_4(iCluster,:)/Coherences_4(iCluster,1);
        case 5
          Coherences_5(iCluster,:) = Coherences_5(iCluster,:)/Coherences_5(iCluster,1);
        case 6
          Coherences_6(iCluster,:) = Coherences_6(iCluster,:)/Coherences_6(iCluster,1);
          
      end
    
    if methyl_number==0
      continue;
    end
    
    % Project beta Hamiltonian onto the methyl's C3 E state.
    Hb_E =   PEap*Hbeta*PEap + PEam*Hbeta*PEam ...
           + PEbp*Hbeta*PEbp + PEbm*Hbeta*PEbm ...
           + PEap*Hbeta*PEbm + PEam*Hbeta*PEbp ...
           + PEbp*Hbeta*PEam + PEbm*Hbeta*PEap;
         
    % Project alpha Hamiltonian onto the methyl's C3 E state.
    Ha_E =   PEap*Halpha*PEap + PEam*Halpha*PEam ...
           + PEbp*Halpha*PEbp + PEbm*Halpha*PEbm ...
           + PEap*Halpha*PEbm + PEam*Halpha*PEbp ...
           + PEbp*Halpha*PEam + PEbm*Halpha*PEap;
            
         
    % Calculate the coherence.
    Coherences_E = propagate(total_time,timepoints,dt,dt2,Ndt,Hb_E,Ha_E,EXPERIMENT,densityMatrix, useThermalEnsemble,betaT);
    
    switch clusterSize
      case 1
        Coherences_1(iCluster,:) = gA*Coherences_1(iCluster,:) + gE*Coherences_E;
      case 2
        Coherences_2(iCluster,:) = gA*Coherences_2(iCluster,:) + gE*Coherences_E;
      case 3
        Coherences_3(iCluster,:) = gA*Coherences_3(iCluster,:) + gE*Coherences_E;
      case 4
        Coherences_4(iCluster,:) = gA*Coherences_4(iCluster,:) + gE*Coherences_E;
      case 5
        Coherences_5(iCluster,:) = gA*Coherences_5(iCluster,:) + gE*Coherences_E;
      case 6
        Coherences_6(iCluster,:) = gA*Coherences_6(iCluster,:) + gE*Coherences_E;
    end
    
     if methyl_number < 2
      continue;
     end
     
     % Project the Hamiltonian onto an A state for  
     % one methyl and an E state for the other.
     
     Hb_AE =   PEap2*Hbeta*PEap2 + PEam*Hbeta*PEam2 ...
       + PEbp2*Hbeta*PEbp2 + PEbm2*Hbeta*PEbm2 ...
       + PEap2*Hbeta*PEbm2 + PEam2*Hbeta*PEbp2 ...
       + PEbp2*Hbeta*PEam2 + PEbm2*Hbeta*PEap2;
     Hb_AE = PA1*Hb_AE*PA1;
     
     Ha_AE =   PEap2*Halpha*PEap + PEam2*Halpha*PEam2 ...
       + PEbp2*Halpha*PEbp + PEbm2*Halpha*PEbm2 ...
       + PEap2*Halpha*PEbm + PEam2*Halpha*PEbp2 ...
       + PEbp2*Halpha*PEam + PEbm2*Halpha*PEap2;
     Ha_AE = PA1*Ha_AE*PA1;

     
     % Calculate the coherence.
     Coherences_AE = propagate(total_time,timepoints,dt,Hb_AE,Ha_AE,EXPERIMENT,densityMatrix, useThermalEnsemble,betaT);
    
     % Project the Hamiltonian onto an A state for  
     % the other methyl and an E state for one.
     Hb_EA =   PEap1*Hbeta*PEap1 + PEam*Hbeta*PEam1 ...
       + PEbp1*Hbeta*PEbp1 + PEbm1*Hbeta*PEbm1 ...
       + PEap1*Hbeta*PEbm1 + PEam1*Hbeta*PEbp1 ...
       + PEbp1*Hbeta*PEam1 + PEbm1*Hbeta*PEap1;
     Hb_EA = PA2*Hb_EA*PA2;
     
     Ha_EA =   PEap1*Halpha*PEap + PEam1*Halpha*PEam1 ...
       + PEbp1*Halpha*PEbp + PEbm1*Halpha*PEbm1 ...
       + PEap1*Halpha*PEbm + PEam1*Halpha*PEbp1 ...
       + PEbp1*Halpha*PEam + PEbm1*Halpha*PEap1;
     Ha_EA = PA2*Ha_EA*PA2;
     
     % Calculate the coherence.
     Coherences_EA = propagate(total_time, timepoints,dt,dt2,Ndt,Hb_EA,Ha_EA,EXPERIMENT,densityMatrix, useThermalEnsemble, betaT);

     % Add the methyl coherences together, weighting the coherences by
     % a statistical factor.
     switch clusterSize
       case 2
         Coherences_2(iCluster,:) = Coherences_2(iCluster,:)  ...
           + Methyl_gA(1)*(1-Methyl_gA(2))*Coherences_AE ...
           +(1 - Methyl_gA(1))*Methyl_gA(2)*Coherences_EA;
       case 3
         Coherences_3(iCluster,:) = Coherences_3(iCluster,:) ...
           + Methyl_gA(1)*(1-Methyl_gA(2))*Coherences_AE ...
           +(1 - Methyl_gA(1))*Methyl_gA(2)*Coherences_EA;
       case 4
         Coherences_4(iCluster,:) = Coherences_4(iCluster,:) ...
           + Methyl_gA(1)*(1-Methyl_gA(2))*Coherences_AE ...
           +(1 - Methyl_gA(1))*Methyl_gA(2)*Coherences_EA;
       case 5
         Coherences_5(iCluster,:) = Coherences_5(iCluster,:) ...
           + Methyl_gA(1)*(1-Methyl_gA(2))*Coherences_AE ...
           +(1 - Methyl_gA(1))*Methyl_gA(2)*Coherences_EA;
       case 6
         Coherences_6(iCluster,:) = Coherences_6(iCluster,:)  ...
           + Methyl_gA(1)*(1-Methyl_gA(2))*Coherences_AE ...
           +(1 - Methyl_gA(1))*Methyl_gA(2)*Coherences_EA;
     end
  end
  
  
end


for clusterSize = 1:Method_order
  switch clusterSize
    case 1
      unusedClusters = find(Coherences_1(:,1)==0)';
      Coherences_1(unusedClusters,:) = 1;
    case 2
      unusedClusters = find(Coherences_2(:,1)==0)';
      Coherences_2(unusedClusters,:) = 1;
    case 3
      unusedClusters = find(Coherences_3(:,1)==0)';
      Coherences_3(unusedClusters,:) = 1;
    case 4
      unusedClusters = find(Coherences_4(:,1)==0)';
      Coherences_4(unusedClusters,:) = 1;
    case 5
      unusedClusters = find(Coherences_5(:,1)==0)';
      Coherences_5(unusedClusters,:) = 1;
    case 6
      unusedClusters = find(Coherences_6(:,1)==0)';
      Coherences_6(unusedClusters,:) = 1;
  end
end

% Calculate signal
%-------------------------------------------------------------------------------
[Signals, AuxiliarySignal_1,AuxiliarySignal_2,AuxiliarySignal_3,AuxiliarySignal_4] ...
  = doClusterCorrelationExpansion_gpu(Coherences_1,Coherences_2,Coherences_3,Coherences_4,Coherences_5,Coherences_6,ClusterArray, ...
  SubclusterIndices_2,SubclusterIndices_3,SubclusterIndices_4,SubclusterIndices_5,SubclusterIndices_6,...
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
% Generate Density Matrix
% ========================================================================
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
function Signal = propagate(total_time,timepoints,dt,dt2,Ndt,Hamiltonian_beta,Hamiltonian_alpha,EXPERIMENT, densityMatrix, useThermalEnsemble, betaT)

% ENUM
FID = 1; HAHN = 2; CPMG = 3; CPMG_CONST = 4; CPMG_2D = 5;

% Ensure subHamiltinians are Hermitian.
Hamiltonian_beta =(Hamiltonian_beta+Hamiltonian_beta')/2; 
Hamiltonian_alpha =(Hamiltonian_alpha+Hamiltonian_alpha')/2;

% Get density matrix.
if useThermalEnsemble
  DensityMatrix = propagator_eig((Hamiltonian_alpha+Hamiltonian_beta)/2,-1i*betaT,0);
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
%       U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
      
      U_ = U_alpha'*U_beta'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha*U_beta;
      v(iTime) = vecDensityMatrixT*U_(:);
      
      % Increment propagator.
%       if iTime<= Ndt
%         U_beta_2 = dU_beta*U_beta_2;
%         U_alpha_2 = dU_alpha*U_alpha_2;
%       else
%         U_beta_2 = dU_beta2*U_beta_2;
%         U_alpha_2 = dU_alpha2*U_alpha_2;
%       end
      
    case CPMG_CONST
      
      % THIS NEEDS TO BE UPDATED TO USE dt2 WHEN iTime > Ndt.
      [U_beta, U_beta_2] = propagator_eig(Hamiltonian_beta,(iTime-1)*dt, total_time/4-(iTime-1)*dt);
      [U_alpha, U_alpha_2] = propagator_eig(Hamiltonian_alpha,(iTime-1)*dt, total_time/4-(iTime-1)*dt);
      
      % Generate time dependent detection operator.
      U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
      
      v(iTime) = vecDensityMatrixT*U_(:);
      
      
    case CPMG_2D
      
      U_beta_2 = eye(nStates);
      U_alpha_2 = eye(nStates);
      
      % Loop over second experimental dimension.
      for jTime = 1:timepoints
        
        % Generate time dependent detection operator.
        U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
        
        v(iTime,jTime) = vecDensityMatrixT*U_(:);
        
        % Increment propagator.
        if jTime < Ndt
          U_beta_2 = dU_beta*U_beta_2;
          U_alpha_2 = dU_alpha*U_alpha_2;
        else
          U_beta_2 = dU_beta2*U_beta_2;
          U_alpha_2 = dU_alpha2*U_alpha_2;
        end
        
      end
      
  end
  % Increment propagator.
  if iTime< Ndt
    U_beta = dU_beta*U_beta;
    U_alpha = dU_alpha*U_alpha;
  else
    U_beta = dU_beta2*U_beta;
    U_alpha = dU_alpha2*U_alpha;
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

