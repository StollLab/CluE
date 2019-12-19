% ========================================================================
% New Function
% ========================================================================
function [Signal, Order_n_Signal,Statistics_ClusterCount, Statistics_NonClusterCount, Statistics_total_clusters, Statistics_total_nonclusters] ...
  = calculateSignal_conserveMemory_gpu(System, Method, Nuclei, timepoints,dt,OutputData,Progress)
%--------------------------------------------------------------------------
% general gpu compatibility
%--------------------------------------------------------------------------

% ENUM
FID = 1; HAHN = 2; CPMG = 3; CPMG_CONST = 4; CPMG_2D = 5;

Method_order = Method.order;  
Method_order_lower_bound = Method.order_lower_bound;
Method_parallelComputing = Method.parallelComputing;
Method_clear_partialSave = Method.clear_partialSave;
System_full_Sz_Hyperfine = System.full_Sz_Hyperfine;
% maxClusterSize = min(4,Method_order);
% Convert variable to gpu compatible forms
dimensionality = 1;
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
  otherwise
  error('The experiment ''%s'' is not supported.',System.experiment);
end
Nuclei_Coordinates = Nuclei.Coordinates;
Nuclei_Abundance = Nuclei.Abundance;
Nuclei_Spin = Nuclei.Spin;
Nuclei_g = Nuclei.Nuclear_g;

Nuclei_startSpin = Nuclei.startSpin;
Nuclei_endSpin = Nuclei.endSpin;
Nuclei_numberStartSpins = Nuclei.numberStartSpins;

NumberStates = Nuclei.NumberStates;
ZeemanStates = Nuclei.ZeemanStates;
max_basis = max(NumberStates);
States = zeros(max_basis,Nuclei.number);
for ii = 1:Nuclei.number
  States(1:NumberStates(ii),ii) = Nuclei.State{ii};
end


HNQ = Nuclei.Qtensor;

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
 
% maxPossibleNumberSubClusters = [0,2,3,6,10,20,35,70,126,252];
maxPossibleNumberSubClusters = zeros(Method_order);
for iorder = 1:Method_order
  for jsize = 1:iorder
    maxPossibleNumberSubClusters(iorder,jsize) = NchooseK(iorder,jsize);  
  end
end

% Define placeholder variables.
ge = System.gMatrix(3,3);
magneticField = System.magneticField;
muB = System.muB;
muN = System.muN;
mu0 = System.mu0;
hbar = System.hbar;

%--------------------------------------------------------------------------
% parfor gpu compatibility
%--------------------------------------------------------------------------

% Determine number of divisions to break calculation into.

Method_divisions_numSpins = false; 

if ~ischar(Method.divisions)
  numDivisions = Method.divisions;
elseif strcmp(Method.divisions,'numSpins')
  
  Method_divisions_numSpins = true;
  numDivisions = double(Nuclei_numberStartSpins);
  
elseif strcmp(Method.divisions,'numCores')
  
  % minimum that will use all nodes
  numDivisions = feature('numcores');
  
else
  numDivisions = round(sqrt(double(Nuclei_numberStartSpins)));
end

Method_shuffle = Method.shuffle;
Method_seed = Method.seed;
Nuclei_ValidPair = Nuclei.ValidPair;
Nuclei_number = Nuclei.number;
% ENUM
CONNECTED = 0;  COMPLETE = 1;

graphCriterion = CONNECTED;
if strcmp(Nuclei.graphCriterion,'complete')
  graphCriterion = COMPLETE;
end

useHamiltonian = System.useHamiltonian;

Method_record_clusters = Method.record_clusters;
Method_partialSave = Method.partialSave;
Method_MonteCarlo_Threshold = Method.MonteCarlo.Threshold;
Method_MonteCarlo_use = Method.MonteCarlo.use;

Possible_Clusters = zeros(Method_order,1);

Division_Cluster_Limit = ones(1,Method_order);
Division_Cluster_Increment = ones(1,Method_order);
Division_Cluster_Fraction = ones(1,Method_order);

for iorder = 1:Method_order  
  Possible_Clusters(iorder) = NchooseK(uint64(Nuclei_number),iorder);  
  if Method.MonteCarlo.use
    Division_Cluster_Limit(iorder) = ceil(Method.MonteCarlo.Cluster_Limit(iorder)/numDivisions);
    Division_Cluster_Increment(iorder) = ceil(Method.MonteCarlo.Increment(iorder)/numDivisions);
    Division_Cluster_Fraction(iorder) = Method.MonteCarlo.Fraction(iorder)/numDivisions;
  else
    Division_Cluster_Limit(iorder) = inf;
    Division_Cluster_Increment(iorder) = inf;
    Division_Cluster_Fraction(iorder) = 1;
  end
end

%--------------------------------------------------------------------------
% Transfer to gpu
%--------------------------------------------------------------------------

Statistics_ClusterCount = zeros(1,Method_order);
Statistics_NonClusterCount = zeros(1,Method_order);
Statistics_total_clusters = Nuclei_numberStartSpins;
Statistics_total_nonclusters = 0;
%--------------------------------------------------------------------------
% gpu code
%--------------------------------------------------------------------------

% initialize output
Signal = ones(1,timepoints^dimensionality);
Order_n_Signal = ones(Method_order,timepoints^dimensionality);

% Initialize internal records.
Cluster_Statistics = zeros(Method_order,numDivisions);
NonCluster_Statistics = zeros(Method_order,numDivisions);

% Get permutation method.
shuffle = 1:Nuclei_number;
% Shuffle = zeros(Nuclei_number);
if Method_shuffle
  
  % Set seed and rng method.
  rng(Method_seed, 'twister');
  
  % A constant seed ensures restarting will not give a new permutation.
  shuffle = randperm(Nuclei_number);
end
  
% Loop over order range.
for iorder = Method_order_lower_bound:Method_order
  
  % Initialize cluster limits.
  numberClusters  = zeros(1,iorder);
  for jorder = 1:iorder
    
    % Get the maximum possible number of subclusters of size iorder. 
    numberClusters(jorder) = nchoosek(iorder,jorder);
  end
  
  % Get cluster array and subcluster indices with elements labeled 1 to iorder. 
  [Reduced_ClusterArray, Reduced_SubclusterIndices_2,...
    Reduced_SubclusterIndices_3,Reduced_SubclusterIndices_4, ...
    Reduced_SubclusterIndices_5,Reduced_SubclusterIndices_6] = getReducedClusters_4(iorder);
  
  % Check if order is valid: n-CCE requires at least n nuclear spins.
  available_nuclei = Nuclei_number - Nuclei_startSpin + 1;
  if available_nuclei < iorder
    continue
  end
  
  % Check pigeonhole principle: each division must be unique.
  pigeonhole_required_nuclei = numDivisions + iorder - 1;
  if pigeonhole_required_nuclei > available_nuclei
    % Reduce the number of divisions to the maximum allowed value.
    numDivisions = available_nuclei - iorder + 1;
  end
  
  % Separate the clusters into bundles for each division.
  Bundle = uint32(ones(numDivisions,2));
  
  if Method_divisions_numSpins
    
    % Assign to the nth division all clusters where the index of the lowest
    % cluster is n.
    
    % Assign bundles while ensuring that there are at least
    % n nuclei with number of nuclei > indices >= n, for n-CCE.
    if Nuclei_endSpin + iorder - 1 > Nuclei_number
      Bundle(:,1) = Nuclei_startSpin:Nuclei_number -iorder + 1;
    else
      Bundle(:,1) = Nuclei_startSpin:Nuclei_endSpin;
    end
    Bundle(:,2) = Bundle(:,1);
    
  else
    
    % For N >> k, to approximately divide N choose k into c even parts,
    % the nth part should get all clusters whose lowest index in the ranges
    % N*[  ( (n-1)/c )^1/k, (n/c)^1/k ].
    
    % Assign the first bundle.
    Bundle(1,1) = 1;
    Bundle(1,2) =  Nuclei_startSpin ...
      + round( ...
      (Nuclei_number- Nuclei_startSpin + 1) ...
      *( 1 - ( 1 - 1/numDivisions)^(1/iorder) ) ...
      *(Nuclei_numberStartSpins/Nuclei_number) ...
      );
    
    % Assign the remaining bundles.
    for iCore = 2:numDivisions
      Bundle(iCore,1) = Bundle(iCore-1,2) + 1;
      Bundle(iCore,2) =  Nuclei_startSpin ...
        + round( ...
        (Nuclei_number- Nuclei_startSpin + 1) ...
        *(1 - (1 - iCore/numDivisions)^(1/iorder)  ) ...
        *(Nuclei_numberStartSpins/Nuclei_number) ...
        );
      
      Bundle(iCore,2) = max(Bundle(iCore,2) , Bundle(iCore,1));
    end
    
    Bundle(numDivisions,2) = Nuclei_endSpin;
  end
  
  % Generate unshuffled starting clusters.
  par_cluster = zeros(numDivisions,iorder);
  for iCore = 1:numDivisions
    par_cluster(iCore,:) = Bundle(iCore,1):(Bundle(iCore,1) + iorder -1);
  end
  
  % initialize division output
  partial_signal = ones(numDivisions,timepoints^dimensionality);
  
  if Method_parallelComputing
    parfor iCore = 1:numDivisions % parfor
    %  [partial_signal(iCore,:),Cluster_Statistics(iorder,iCore),NonCluster_Statistics(iorder,iCore)] = conserve_memory_for_loop(System,Method,iorder,OutputData,Nuclei,par_cluster{iCore},Bundle,iorder,iCore,Division_Cluster,shuffle, linearTimeAxis,numDivisions);
    
    [partial_signal(iCore,:),Cluster_Statistics(iorder,iCore),NonCluster_Statistics(iorder,iCore)] = ...
      conserve_memory_gpu_loop( ...
      OutputData, par_cluster(iCore,:), Bundle, iorder, iCore, ...
      Division_Cluster_Limit, Division_Cluster_Fraction, Division_Cluster_Increment, ...
      shuffle, timepoints,dt, ... % Method_order,
      EXPERIMENT, dimensionality, ...
      System_full_Sz_Hyperfine,total_time, ...
      Nuclei_Coordinates, Nuclei_ValidPair, Nuclei_Abundance, Nuclei_Spin, Nuclei_g, NumberStates, ZeemanStates, ...
      Spin2Op1, Spin2Op2, Spin2Op3, Spin2Op4, Spin2Op5, Spin2Op6,...
      Spin3Op1, Spin3Op2, Spin3Op3, Spin3Op4, Spin3Op5, Spin3Op6, ...
      numberClusters, ... % ClusterArray,SubclusterIndices_2,SubclusterIndices_3,SubclusterIndices_4, ...
      ge, magneticField, muB, muN, mu0, hbar, ...
      Method_seed, Nuclei_number,graphCriterion,Method_partialSave,Method_MonteCarlo_Threshold,maxPossibleNumberSubClusters, ...
      Reduced_ClusterArray,...
      Reduced_SubclusterIndices_2,Reduced_SubclusterIndices_3,Reduced_SubclusterIndices_4, ...
      Reduced_SubclusterIndices_5,Reduced_SubclusterIndices_6, ...
      Method_record_clusters,useHamiltonian);
    
    end
  else
    for iCore = 1:numDivisions % parfor
      % [partial_signal(iCore,:),Cluster_Statistics(iorder,iCore),NonCluster_Statistics(iorder,iCore)] = conserve_memory_for_loop(System,Method,iorder,OutputData,Nuclei,par_cluster{iCore},Bundle,iorder,iCore,Division_Cluster,shuffle, linearTimeAxis,numDivisions);
      
      [partial_signal(iCore,:),Cluster_Statistics(iorder,iCore),NonCluster_Statistics(iorder,iCore)] = ...
        conserve_memory_gpu_loop( ...
        OutputData, par_cluster(iCore,:), Bundle, iorder, iCore, ...
        Division_Cluster_Limit, Division_Cluster_Fraction, Division_Cluster_Increment, ...
        shuffle, timepoints,dt, ... % Method_order,
        EXPERIMENT, dimensionality, ...
        System_full_Sz_Hyperfine,total_time, ...
        Nuclei_Coordinates, Nuclei_ValidPair, Nuclei_Abundance, Nuclei_Spin, Nuclei_g, NumberStates, ZeemanStates, ...
        Spin2Op1, Spin2Op2, Spin2Op3, Spin2Op4, Spin2Op5, Spin2Op6,...
        Spin3Op1, Spin3Op2, Spin3Op3, Spin3Op4, Spin3Op5, Spin3Op6, ...
        numberClusters, ... % ClusterArray,SubclusterIndices_2,SubclusterIndices_3,SubclusterIndices_4, ...
        ge, magneticField, muB, muN, mu0, hbar, ...
        Method_seed, Nuclei_number,graphCriterion,Method_partialSave,Method_MonteCarlo_Threshold,maxPossibleNumberSubClusters, ...
        Reduced_ClusterArray,...
        Reduced_SubclusterIndices_2,Reduced_SubclusterIndices_3,Reduced_SubclusterIndices_4, ...
        Reduced_SubclusterIndices_5,Reduced_SubclusterIndices_6, ...
        Method_record_clusters,useHamiltonian);
      
    end
  end

  % combine node outputs to the final output
  
  
  Combined_AuxiliarySignal = ones(1,timepoints^dimensionality);
  for iCore = 1:numDivisions
    
    % main signal
    Combined_AuxiliarySignal = Combined_AuxiliarySignal.*partial_signal(iCore,:);
    % get signals for all n-CCE, with n<= Method_order
    %       for jclusterSize = Method_order:-1:iorder
    %         Order_n_Signal{jclusterSize} = Order_n_Signal{jclusterSize}.*partial_signal(iCore,:);
    %       end
    
  end
  
  
  Statistics_ClusterCount(iorder) = sum(Cluster_Statistics(iorder,:) );
  Statistics_NonClusterCount(iorder) = sum(NonCluster_Statistics(iorder,:) );
  if Method_MonteCarlo_use
    clusterFraction = (Statistics_NonClusterCount(iorder) + Statistics_ClusterCount(iorder))/Possible_Clusters(iorder);
    Combined_AuxiliarySignal = Combined_AuxiliarySignal.^(1/clusterFraction);
  end
  
  Signal = Signal.*Combined_AuxiliarySignal;
  Order_n_Signal(iorder,:) = Signal;
  
  % Save to file.
  Progress.order = [num2str(iorder) '/' num2str(Method_order) '-CCE'];
  save(OutputData,'Signal','Order_n_Signal','Progress','-append');
  
  if Method_clear_partialSave
    for iCore = 1:numDivisions
      partial_file = ['partial_', OutputData(1:end-4), '_',num2str(iorder),'cce_', num2str(iCore), '.mat'] ;
      if isfile(partial_file)
        delete(partial_file);
      end
    end
  end
  
end
  
for iorder = 2:Method_order
  %     Statistics_ClusterCount(iorder) = sum(Cluster_Statistics(iorder,:) );
  %     Statistics_NonClusterCount(iorder) = sum(NonCluster_Statistics(iorder,:) );
  
  Statistics_total_clusters = Statistics_total_clusters + Statistics_ClusterCount(iorder);
  Statistics_total_nonclusters = Statistics_total_nonclusters + Statistics_NonClusterCount(iorder);
  
end

end
