% CluE   Calculate spin dynamics using cluster-based expansions in Hilbert space
%
% Inputs:
%   System is a struct with fields: magneticField, Electron, temperature, Time.
%     Electron is a struct with fields: spin, g, Coordinates, wavefunction.
%   Method is a struct with fields: method, r0, order.
%   Data is a structure containing
%    .InputData  name of pdb file
%    .OutputData name of mat file to store results
%    .saveLevel  0 (standard), 1 (more), 2 (all)

function [SignalMean, experiment_time, ...
  TM_powder,Order_n_SignalMean,Nuclei,statistics] = CluE(System,Method,Data)

tic

% Set defaults base on specified parameters and for unspecified parameters
[System, Method, Data,statistics] = setDefaults(System,Method,Data);

[OutputData,doReturn,SignalMean,experiment_time,...
  TM_powder,Order_n_SignalMean,Nuclei,uncertainty] = setOutput(Data,Method);
if doReturn
  toc
  return
end

% Set verbosity
verbose = Method.verbose;
if verbose, fprintf('Setting up...\n');  end

if ~isempty(Method.seed)
  if verbose
    fprintf('Using rng seed to Method.seed = %d.\n',Method.seed);
  end
  rng(Method.seed);
end

% Time axis
switch System.experiment
  case 'FID'
    experiment_time = System.Time;
  case 'Hahn'
    experiment_time = 2*System.Time;
  case {'CPMG','CPMG-2D'}
    experiment_time = 4*System.Time;% + 2*System.Time_;
  case 'CPMG-const'
    experiment_time = 2*System.Time;
  case 'CP_N'
    experiment_time = 2*System.nPulses*System.Time;
  case 'Uhrig_N'
    experiment_time = 2*System.nPulses*System.Time;
end

% Save input data
if ~isfield(Data,'InputData') || ~ischar(Data.InputData)
  error('Data.InputData is not a valid file identifier.');
end
InputData = Data.InputData;
if ~isfile(InputData) 
  if System.RandomEnsemble.include
    InputData = 'System.RandomEnsemble.include';
  else
    error(['Could not find the specified input file ', InputData, '.']);
  end
end

% Initiate progress tracking
Progress.started = true;
Progress.complete = false;
Progress.Completed_Orders = zeros(1,Method.order);

% Save if filename is provided
if ~isempty(OutputData) && Data.saveLevel>=0
  Input.System = System;
  Input.Method = Method;
  Input.Data = Data;
  save(OutputData,'Input','experiment_time','Progress');
end

if strcmp(InputData(end-3:end),'.mat') % check to see if InputData is a saved file.
  newMethod = Method.method;
  load(InputData);
  Method.method = newMethod;
  clear newMethod
  Progress.LoadSavedData = true;
  
elseif min( (InputData(end-3:end)) == '.pdb') || strcmp(InputData,'user') ...
    || strcmp(InputData,'System.RandomEnsemble.include')
  
  if Method.useCentralSpinSystem
    System.pdb = parsePDBfile(Data.InputData, System.angstrom);


    if ~isempty(System.pdbTranslation)
      System.pdb = translatePDB(System.pdb,System);
    end

    if System.pdbRotate
      System.pdb = rotatePDB(System.pdb,System);
    end

  end
  [Nuclei, System] = parseNuclei(System, Method, Data, InputData);
  
  if Nuclei.number < 1
    Signals{1} = ones(size(System.Time));
    SignalMean = Signals{1};
    fprintf(2,'\n There are too few spins in the system for spin decoherence.\n')
    fprintf(2,' Try relaxing one or more cutoffs.\n')
    toc
    return
  end
  
else
  error('Input data in Data.InputData not recognized.')
end

Progress.DataLoaded = true;

if verbose, fprintf('Setup initialized system with %i nuclei.\n', Nuclei.number); end

if Nuclei.number < Method.order
  fprintf('Reducing the maximum cluster size to the system size of %d.\n',Nuclei.number)
  Method.order = double(Nuclei.number);
end

% Only possible for very small systems
if strcmp(Method.method,'full')
  Method.order = Nuclei.number;
end


% ========================================================================
% Compile list of connected clusters
% ========================================================================
Clusters = cell(Method.extraOrder,1);
for clusterSize = 1:Method.extraOrder
  Clusters{clusterSize} = [];
end

inClusters = {};
doFindClusters = ~Method.Ori_cutoffs;
if ~isempty(Data.ClusterData)  || isfield(Method,'Clusters')
  Clusters = {};
  if ~isempty(Data.ClusterData)
    try
      load(Data.ClusterData,'Clusters');      
    catch
      disp('Could not load clusters.')      
      if Data.exitOnFailedLoad
        error('Could not load clusters.')
      end
    end
  else
    Clusters = Method.Clusters;
  end
   
  inOrder = numel(Clusters); 
   
  for clusterSize = 1:inOrder
    
    Nuclei.numberClusters(clusterSize) = size(Clusters{clusterSize},1);
    if verbose
      fprintf('Loaded %d clusters of size %d.\n', Nuclei.numberClusters(clusterSize),clusterSize);
    end
   
  end
  if inOrder < Method.order
    inClusters = Clusters;
  else
    doFindClusters = isempty(Clusters);
  end
end


% Loop over cluster sizes, start at the largest (most time consuming) size
if doFindClusters
  
  if strcmp(Method.clusterization,'tree-search')
    if verbose
      fprintf('Finding clusters of up to size %d.\n',Method.extraOrder);
    end

    if Method.combineClusters
      Clusters = findClusters_treeSearch(Nuclei,Method.extraOrder,1,{});
      for clusterSize = 1:min(Method.order, numel(inClusters),Method)
        % Combine arrays.
        C = [Clusters{clusterSize}; inClusters{clusterSize}];
        
        % Sort clusters
        C = sortrows(C);
        
        % Remove duplicates
        keep = [true; any(C(1:end-1,:)~=C(2:end,:),2)];
        Clusters{clusterSize} = C(keep,:);
      end
    else
      if ~Method.cutoff.sizeDependent
      Clusters = findClusters_treeSearch(Nuclei,Method.extraOrder,1,...
        inClusters, Method);
      else
        Clusters = cell(1,Method.extraOrder);
        for adjacencyOrder = Method.extraOrder:-1:1
          Clusters_ = findClusters_treeSearch(Nuclei,Method.extraOrder,...
            adjacencyOrder,inClusters, Method);
          Clusters{adjacencyOrder} = Clusters_{adjacencyOrder};
        end
      end
    end
    Nuclei.numberClusters = zeros(1,Method.extraOrder);
    
    for clusterSize = 1:Method.extraOrder
      Nuclei.numberClusters(clusterSize) = size(Clusters{clusterSize},1);
      
      if verbose
        fprintf('  Found %d clusters of size %d.\n', Nuclei.numberClusters(clusterSize),clusterSize);
      end
    end
    
  else
    for clusterSize = Method.order:-1:1
      
      if strcmp(Method.method,'full')
        continue
      end
      
      if verbose, fprintf('Finding clusters of size %d.\n', clusterSize); end
      
      Clusters{clusterSize} = findClusters(Nuclei, clusterSize);
      
      % save the number of clusters of each size for export to the user
      Nuclei.numberClusters(clusterSize) = size(Clusters{clusterSize},1);
      
      if verbose
        fprintf('  Found %d clusters of size %d.\n', Nuclei.numberClusters(clusterSize),clusterSize);
      end
      
    end
  end
end

% Save Memory.
Nuclei.Adjacency = [];

% Save clusters and exit without calculation if requested
if strcmp(Method.method,'count clusters')
  SignalMean = Nuclei.numberClusters;
  experiment_time = 1:length(SignalMean);
  TM_powder = [];
  Order_n_SignalMean = [];
  if Method.exportClusters
    if ~isempty(OutputData)
      clusterSaveFile = [OutputData(1:end-4),'Clusters.mat'];
    else
      clusterSaveFile = 'Clusters.mat';
    end
    save(clusterSaveFile,'Clusters');
  end    
  toc
  return
end


% ========================================================================
% Starting main loop.
% ========================================================================
if verbose
  fprintf('  Computing cluster Hamiltonians and cluster signals.\n');
end
if Method.parallelComputing  
  % remove current pool if it exists.
  delete(gcp('nocreate'));
  
  % determine number of cores available
  numCores = min(feature('numcores'), Method.numberCores);
  
  % create parallel pool
  pc = parcluster('local');

  if isfield(Method,'JobStorageLocation')
    pc.JobStorageLocation = options.JobStorageLocation;
  elseif Method.slurm
    pc.JobStorageLocation = ...
      strcat(getenv('SCRATCH'),'/', getenv('SLURM_JOB_ID'));
  end 

  pool = parpool(numCores);
else
  pool = [];
end

%===============================================================================
% Set up powder averaging.
%===============================================================================
% Powder average settings
gridSize = System.gridSize;
GridInfo = [];
if gridSize==1
  System.averaging = 'none';
end

if strcmp(System.averaging,'powder')
  
  Grid = lebedev_grid(gridSize);
  % Use the inversion symmetry of the spin-Hamiltonian to skip all the
  
  removeGridPoints = 1:gridSize;
  removeGridPoints(Grid.z < 0) = -1;
  removeGridPoints( ((Grid.z==0) & (Grid.x<0)) ) = -1;
  removeGridPoints( ((Grid.z==0) & (Grid.x==0)) & (Grid.y<0) ) =-1;
  keep = removeGridPoints>0;
  
  Grid.w = Grid.w(keep);
  Grid.w = Grid.w/sum(Grid.w);
  
  Grid.x = Grid.x(keep);
  Grid.y = Grid.y(keep);
  Grid.z = Grid.z(keep);
  
  % Convert xyz coordinates to alpha/beta angles
  Alpha(numel(Grid.z)) = 0;
  Beta(numel(Grid.z)) = 0;
  for iOri = 1:numel(Grid.z)
    
    % Determine the beta Euler angle.
    Beta(iOri) = acos(Grid.z(iOri));
    
    % Determine the alpha Euler angle.
    if (Grid.x(iOri)^2+Grid.y(iOri)^2)>0 % check for the pole singularities.
      
      % Determine which quadrant of the xy-plane the grid point is in.
      quadZ = sign(Grid.y(iOri));
      if abs(quadZ) < 0.1
        quadZ = 1;
      end
      Alpha(iOri) = quadZ*acos(Grid.x(iOri)/sqrt(Grid.x(iOri)^2 + Grid.y(iOri)^2));
    else % set alpha to zero for beta = 0 or pi.
      Alpha(iOri) = 0;
    end
    
  end
  gridSize = length(Alpha);
  gridWeight = Grid.w;
  
  GridInfo.gridSize   = gridSize;
  GridInfo.gridWeight = gridWeight;
  GridInfo.Alpha      = Alpha;
  GridInfo.Beta       = Beta;
  
elseif strcmp(System.averaging,'none')
  
  % Use only the PDB file orientation.
  gridSize = 1;
  Alpha = 0;
  Beta = 0;
  gridWeight = 1;
  
elseif strcmp(System.averaging,'xy')
  
  % Average over rotations about the B0 direction.
  Alpha = linspace(0,pi,gridSize+1);
  Alpha(end) = [];
  Beta = ones(1,gridSize)*pi/2;
  gridWeight = ones(gridSize,1)/gridSize;
  
elseif strcmp(System.averaging,'custom')

  if isstruct(System.Grid)
    Alpha = System.Grid.Alpha;
    Beta = System.Grid.Beta;
    gridWeight = System.Grid.gridWeight;  
    gridSize = length(Alpha);
  end

elseif strcmp(System.averaging,'Nitroxide_Wband_Weights')

  GridSizes = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, ...
         266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,...
         3074, 3470, 3890, 4334, 4802, 5294, 5810];
  gridIndex = find(GridSizes==gridSize);
  
  load([Data.path2CluE, 'grids/Lebedev_weighted_Nitroxide_Wband_Weights.mat'],'Grids');
  Alpha = Grids{gridIndex}.Alpha;
  Beta = Grids{gridIndex}.Beta;
  gridWeight = Grids{gridIndex}.Weight;

  gridSize = length(Alpha);
  
  GridInfo.gridSize   = gridSize;
  GridInfo.gridWeight = gridWeight;
  GridInfo.Alpha      = Alpha;
  GridInfo.Beta       = Beta;
  
elseif strcmp(System.averaging,'random')
  
  % Generate random Euler angles.
  Alpha = rand(1,gridSize)*2*pi;
  Beta = acos(2*rand(1,gridSize)-1);
  gridWeight = ones(gridSize,1)/gridSize;

end
nOrientations = gridSize;
if verbose
  fprintf("Calculating signal for %d orientations.\n",nOrientations);
end

nTimePoints = sum(System.nPoints);

% initialize result variables
Signals{nOrientations} = [];
TM(nOrientations) = 0;
AuxiliarySignal{nOrientations} = [];
Calculate_Signal{nOrientations} = true;
Ori_Clusters = cell(nOrientations,1);

% update progress
Progress.Order_n_Mean = 'pending';
Order_n_Signals{nOrientations} = [];
if strcmp(Method.method,'HD-CCE')
  SignalMean = zeros(1,nTimePoints^System.dimensionality,length(System.deuteriumFraction));
else
SignalMean = zeros(1,nTimePoints^System.dimensionality);
end
% initialize
for iorder = Method.order:-1:1
  if strcmp(System.experiment,'CPMG-2D')
  Order_n_SignalMean{iorder} = zeros(1,nTimePoints^2);
  else
    Order_n_SignalMean{iorder} = zeros(1,nTimePoints);
  end
end

SignalsToCalculate = zeros(nOrientations,1);
counter_SignalsToCalculate = 0;
saveAll = Data.saveLevel==2;
Statistics = cell(nOrientations,1);
% Check to see if file already exists.
for iOri = 1:nOrientations
  
  % temporary file for partial saving
  temp_file = ['temp_', OutputData, '_sig_', num2str(iOri), '.mat'] ;
  
  Calculate_Signal{iOri} = true;
  
  % check if file alread exists
  if Method.partialSave && isfile(temp_file) && ~(saveAll || Method.getContributions)
    try
      % load partial save
      if Method.Ori_cutoffs
        load(temp_file,'signal','order_n','seed','statistics','iOri_Clusters');
      else
        load(temp_file,'signal','order_n','seed','statistics');
      end
      % check progress
      if seed ~= Method.seed %&& progress_powder
        error('RNG seeds do not match.  Loaded data may be inconsistant.')
      end  
      % use loaded data
      if verbose, fprintf(['Loading signal %d from ', temp_file, '.\n'],iOri); end
      
      % set values to simulation variables
      Calculate_Signal{iOri} = false;
      
      if Method.sparseMemory
        SignalMean = SignalMean + signal;
        for iorder = Method.order:-1:1
          Order_n_SignalMean{iorder} = Order_n_SignalMean{iorder} + order_n{iorder};
        end
      else
        Signals{iOri} = signal;
        Order_n_Signals{iOri} = order_n;
        Statistics{iOri} = statistics;
        if Method.Ori_cutoffs
          Ori_Clusters{iOri} = iOri_Clusters;
        end
      end
      

    catch     
      counter_SignalsToCalculate = counter_SignalsToCalculate + 1;
      SignalsToCalculate(counter_SignalsToCalculate) = iOri;
    end
  else
    counter_SignalsToCalculate = counter_SignalsToCalculate + 1;
    SignalsToCalculate(counter_SignalsToCalculate) = iOri;
  end
  
end
SignalsToCalculate(SignalsToCalculate==0) = [];
numberOfSignals = length(SignalsToCalculate);

TempSignals{numberOfSignals+1} = [];
Temp_Order_n_Signals{numberOfSignals+1} = [];


graphs = cell(numberOfSignals,1);

parallelComputing = Method.parallelComputing;


if parallelComputing
  
  parfor iOri = 1:numberOfSignals
    
    [TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics{iOri},graphs{iOri},Ori_Clusters{iOri}] ...
      = getOrientationSignals(System,Method,Nuclei,Clusters, Alpha,Beta, ...
      iOri,verbose,OutputData,Data,InputData,SignalsToCalculate,gridWeight,nOrientations);
    
    TempSignals{iOri} = TempSignals_;
    Temp_Order_n_Signals{iOri} = Temp_Order_n_Signals_;
    
    if saveAll || Method.getContributions
      AuxiliarySignal{iOri} = AuxiliarySignal_;
    end
  end
  
else
  
  for iOri = 1:numberOfSignals
    
    [TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics{iOri},graphs{iOri},Ori_Clusters{iOri}] ...
      = getOrientationSignals(System,Method,Nuclei,Clusters, Alpha,Beta, ...
      iOri,verbose,OutputData,Data,InputData,SignalsToCalculate,gridWeight,nOrientations);
    
    TempSignals{iOri} = TempSignals_;
    Temp_Order_n_Signals{iOri} = Temp_Order_n_Signals_;
    
    if Method.sparseMemory
    
      SignalMean = SignalMean + TempSignals_;
      for iorder = Method.order:-1:1
        Order_n_SignalMean{iorder} = Order_n_SignalMean{iorder} + Temp_Order_n_Signals_{iorder};
      end
      
    else
      
      TempSignals{iOri} = TempSignals_;
      Temp_Order_n_Signals{iOri} = Temp_Order_n_Signals_;
      
      if saveAll || Method.getContributions
        AuxiliarySignal{iOri} = AuxiliarySignal_;
      end
      
    end
    
  end
  
end

tempSignalIndex = 0;
if ~Method.sparseMemory 
  for iOri = 1:nOrientations
    if Calculate_Signal{iOri} 
      
      tempSignalIndex = tempSignalIndex  +1;
      Signals{iOri} = TempSignals{tempSignalIndex};
      Order_n_Signals{iOri} = Temp_Order_n_Signals{tempSignalIndex};
    end
    if ~strcmp(Method.method,'HD-CCE')
      TM(iOri) = getTM(experiment_time,Signals{iOri});
    end
  end
  
  clear('TempSignals');
end

if ~isempty(pool)
  delete(pool);
end

if System.doPruneNuclei && ~Method.reparseNuclei
  statistics = Nuclei.Isotopologue;
else
  statistics.Statistics = Statistics;
  statistics.uncertainty = uncertainty;
end
Nuclei.graphs = graphs;
% ========================================================================
% End of main loop
% ========================================================================

if ~Method.sparseMemory 
  if System.newIsotopologuePerOrientation
    Nuclei.PowderStatistics.Isotopologue.Mean_TypeNumber = [0,0,0];
    Nuclei.PowderStatistics.Isotopologue.Mean_Instance_2H_Number = [0,0,0];
    Nuclei.PowderStatistics.Isotopologue.Mean_Instance_2H_Fraction = [0,0,0];
  end
  
  for iOri = 1:nOrientations
    if ~Method.reparseNuclei && System.newIsotopologuePerOrientation
      Nuclei.PowderStatistics.Isotopologue.Mean_TypeNumber = Nuclei.PowderStatistics.Isotopologue.Mean_TypeNumber ...
        + Statistics{iOri}.Isotopologue.TypeNumber/nOrientations;
      
      Nuclei.PowderStatistics.Isotopologue.Mean_Instance_2H_Number = Nuclei.PowderStatistics.Isotopologue.Mean_Instance_2H_Number ...
        + Statistics{iOri}.Isotopologue.Instance_2H_Number/nOrientations;
      
      Nuclei.PowderStatistics.Isotopologue.Mean_Instance_2H_Fraction = Nuclei.PowderStatistics.Isotopologue.Mean_Instance_2H_Fraction ...
        + Statistics{iOri}.Isotopologue.InstanceFraction/nOrientations;
    end
    if Method.Ori_cutoffs
      Clusters =  combineClusters(Clusters,Ori_Clusters{iOri});
    end
    SignalMean = SignalMean + Signals{iOri};
    %Signals{isignal} = abs(Signals{isignal});
    if ~ischar(Order_n_SignalMean{iorder}) && ~isempty(Order_n_Signals{1})
      for iorder = 1:Method.order
        try
          Order_n_SignalMean{iorder} = Order_n_SignalMean{iorder} + Order_n_Signals{iOri}{iorder};
        catch
          fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
          
          fprintf('Error in Order %d mean at signal number %d.\n',iorder,iOri);
          disp('Could not evaluate');
          disp('Order_n_SignalMean{iorder} = Order_n_SignalMean{iorder} + Order_n_Signals{isignal}{iorder};');
          fprintf('Order_n_SignalMean{%d} = Order_n_SignalMean{%d} + Order_n_Signals{%d}{%d};\n', ...
            iorder,iorder,iOri,iorder);
          if exist('Order_n_Signals','var')
            if iscell(Order_n_Signals) && length(Order_n_Signals)>=iOri
              
              if iscell(Order_n_Signals{iOri}) && length(Order_n_Signals{iOri}) >= iorder
                fprintf('Length of Order_n_Signals{%d}{%d} is %d.\n',...
                  iOri,iorder, length(Order_n_Signals{iOri}{iorder}));
                
                fprintf('Length of Order_n_SignalMean{%d} is %d.\n',...
                  iorder, length(Order_n_SignalMean{iorder})  );
                
              else
                fprintf('Order_n_Signals{%d}{%d} does not exist.\n',iOri,iorder);
              end
              
            else
              fprintf('Order_n_Signals{%d} does not exist.\n',iOri);
            end
          else
            fprintf('Order_n_Signals does not exist.\n');
          end
          
          disp('Attempting to recover data...')
          
          
          
          if iorder == 1
            
            Order_n_Signals{iOri}{iorder} = gridWeight(iOri)*ones(size(experiment_time));
            
            disp('Recovered.');
            
          elseif iorder == Method.order
            Order_n_Signals{iOri}{iorder} = Signals{iOri};
            disp('Recovered.');
          else
            Progress.Order_n_Mean = false;
            disp('Failed.');
            Order_n_Signals{iOri}{iorder} = nan(size(experiment_time));
          end
          
          Order_n_SignalMean{iorder} = Order_n_SignalMean{iorder} + Order_n_Signals{iOri}{iorder};
          fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        end
      end
    end
    % Signals{isignal} = Signals{isignal}/max(Signals{isignal});
  end
  for iSize=1:Method.order
    Nuclei.numberClusters(iSize) = size(Clusters{iSize},1);
  end
end

% update progress
if Progress.Order_n_Mean, Progress.Order_n_Mean = true; end

if System.dimensionality==2 && min(size(SignalMean))==1
  SignalMean = reshape(SignalMean',nTimePoints,nTimePoints)';
  for iorder = 1:Method.order
    Order_n_SignalMean{iorder} = reshape(Order_n_SignalMean{iorder}',nTimePoints,nTimePoints)';
  end
end

% Save clusters if requested
if Method.exportClusters
  if ~isempty(OutputData)
    clusterSaveFile = [OutputData(1:end-4),'Clusters.mat'];
  else
    clusterSaveFile = 'Clusters.mat';
  end
  save(clusterSaveFile,'Clusters','Ori_Clusters');
end


if strcmp(Method.method,'count clusters')
  SignalMean = SignalMean(1:Method.order);
  experiment_time = 1:Method.order;
  SignalMean = SignalMean/SignalMean(1)*double(Nuclei.number);
end
TM_powder = getTM(experiment_time,abs(SignalMean));
Progress.complete = true;

if verbose
  fprintf('Processing\n');
end
uncertainty = [];
% decide what to save
if ~isempty(OutputData)
 
  T = array2table([experiment_time',SignalMean']);
  T.Properties.VariableNames(1:2) = {'time','signal'};
  writetable(T,[OutputData(1:end-4) , '.csv']);

  switch Data.saveLevel
    case 0
      if Method.sparseMemory
        save(OutputData,'SignalMean','Order_n_SignalMean','TM_powder',...
            'Progress','uncertainty','GridInfo','-append');
      else
        save(OutputData,'SignalMean','Signals','TM','TM_powder',...
            'Progress','uncertainty','GridInfo','-append');
      end
    case 1
      if Method.sparseMemory
        save(OutputData,'Nuclei','Order_n_Signals','GridInfo','-append');
      else
        save(OutputData,'SignalMean','Signals','TM','TM_powder','Progress',...
          'Nuclei','Order_n_SignalMean','Order_n_Signals','uncertainty',...
          'GridInfo','-append');
      end
    case 2
      save(OutputData,'-v7.3');
  end
end

% Delete temporary files
for iOri = 1:nOrientations
  temp_file = ['temp_', OutputData, '_sig_', num2str(iOri), '.mat'] ;
  if isfile(temp_file)
    delete(temp_file);
  end
end

% Calculate the contributions from each spin
% placed after the save since getSpinContributions() is still buggy
%NuclearContribution.findContributions = Method.getNuclearContributions;

if Method.getUncertainty
  uncertainty = getClusterError([], ...
    Nuclei, System, nOrientations, Clusters, AuxiliarySignal,Method, ...
    gridWeight, Order_n_SignalMean, Order_n_Signals); 
  save(OutputData,'uncertainty','-append');
end
if Method.getNuclearSpinContributions
  getNuclearSpinContributions([OutputData,'SpinContribution.mat'], ...
    Nuclei, System, nOrientations, Clusters, Signals, AuxiliarySignal, ...
    Method, experiment_time, gridWeight, TM_powder, Input, SignalMean, ...
    Order_n_SignalMean);
end
if Method.getClusterContributions
  if ~Method.Ori_cutoffs
  getClusterContributions([OutputData,'ClusterContribution.mat'], ...
    Nuclei, System, 1:nOrientations, Clusters, Signals, AuxiliarySignal, ...
    Method, experiment_time, gridWeight, TM_powder, Input, SignalMean, ...
    Order_n_SignalMean,Order_n_Signals);
  else
    for iOri = 1:nOrientations

      getClusterContributions(...
        [OutputData,'ClusterContribution_Ori_',num2str(iOri), '.mat'], ...
        Nuclei, System, iOri, Ori_Clusters{iOri}, Signals, AuxiliarySignal, ...
        Method, experiment_time, gridWeight, TM_powder, Input, SignalMean, ...
        Order_n_SignalMean,Order_n_Signals);
    end
  end
end

if isfield(Method, 'getNuclearContributions') && Method.getNuclearContributions
  if isnan(TM)
    fprintf('Cannot find nuclear contributions whe TM is Nan.');
  else
    NuclearContribution = getSpinContributions(System, Nuclei, Signals, ...
      AuxiliarySignal, Clusters,Method,OutputData,gridWeight, TM_powder);
    save(OutputData,'NuclearContribution','-append');
  end
end

if verbose
  fprintf('\nCompleted Nuclear Spin Diffusion\n');
  fprintf('\nNuclear spin decoherence time = %f µs. \n',TM_powder/1e-6);
end

if Method.conserveMemory
  if isfield(statistics,'Statistics')
    for jj = 1:numel(statistics.Statistics)
      if isfield(statistics.Statistics{jj},'Nuclear_Dipole')
        statistics.Statistics{jj} = ...
          rmfield(statistics.Statistics{jj},'Nuclear_Dipole');
      end
      if isfield(statistics.Statistics{jj},'Nuclear_Dipole_x_iy_Z')
        statistics.Statistics{jj} = ...
          rmfield(statistics.Statistics{jj},'Nuclear_Dipole_x_iy_Z');
      end

    end
  end
end


toc
end
% ========================================================================
% ========================================================================



% ========================================================================
% Calculates signal for a set of orientations
% ========================================================================
function [TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,...
  Statistics_isignal,graphs_isignal,iOri_Clusters] ...
    = getOrientationSignals(System,Method,Nuclei,inClusters, Alpha,Beta,...
    isignal,verbose,OutputData,Data,InputData,SignalsToCalculate,gridWeight,iSignal_max)
     
  
  % grid point indices
  if System.newIsotopologuePerOrientation && ~Method.reparseNuclei

    if Method.useCentralSpinSystem
      error(['Error in getOrientationSignals(): ',...
        'Method.useCentralSpinSystem is not yet compatible with',...
        ' System.newIsotopologuePerOrientation.']);
    end
    Nuclei = newHydronIsotopologue(Nuclei,System);
    if verbose
      fprintf('Generated a new hydron isotopologue for orientation %d.\n',isignal)
    end
    
    Clusters = inClusters;
    
  else
    Clusters = inClusters;
  end
  igrid = SignalsToCalculate(isignal);

if isempty(Clusters{1})
%   error(['Error in getOrientationSignals(): ',...
%     'empty cluster set.']);
end

% calculate coherence signal
[TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics_isignal,...
  graphs_isignal,iOri_Clusters] = ...
  beginCalculateSignal(System,Method,Nuclei,Clusters,...
  Alpha(igrid),Beta(igrid),verbose,OutputData,Data,InputData,isignal);

if System.newIsotopologuePerOrientation  && ~Method.reparseNuclei
  Statistics_isignal.Isotopologue = Nuclei.Isotopologue;
end

normalizing_factor = TempSignals_(1);
TempSignals_ = gridWeight(igrid)*TempSignals_/normalizing_factor;
if strcmp(Method.method,'full')
  return
end

% Order n signals
if ~ischar(Temp_Order_n_Signals_)
  for iorder = 1:Method.order
    % set the max amplitude to the weight
    Temp_Order_n_Signals_{iorder} = gridWeight(igrid)*Temp_Order_n_Signals_{iorder}/normalizing_factor;
  end
end

if verbose, fprintf('\nCompleted orientation %d/%d.\n',igrid,iSignal_max); end

% Save to file.
if Method.partialSave
  temp_file = ['temp_', OutputData, '_sig_', num2str(igrid), '.mat'] ;
  parsavefile = matfile(temp_file,'writable',true);
  parsavefile.signal = TempSignals_;
  parsavefile.order_n = Temp_Order_n_Signals_;
  parsavefile.progress_powder = true;
  parsavefile.seed = Method.seed;
  parsavefile.statistics = Statistics_isignal;
  if Method.Ori_cutoffs
    parsavefile.iOri_Clusters = iOri_Clusters;
  end
end

end

% ========================================================================
% New Function
% ========================================================================

function Indices = findSubclusters(Clusters,clusterSize,iCluster)
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
    Subindices = findSubclusters(Clusters,clusterSize -1,subclusterIndex);
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
%{
function DensityMatrix = getDensityMatrix(Nuclei, Cluster)

dimension = prod(Nuclei.NumberStates(Cluster));

DensityMatrix = eye(dimension);

switchIndex = dimension;
for inucleus = Cluster
  switchIndex = switchIndex/Nuclei.NumberStates(inucleus);
  counter = 0;
  stateNumber = 1;
  for jj = 1:dimension
    if counter == switchIndex
      stateNumber = mod(stateNumber,Nuclei.NumberStates(inucleus)) + 1;
    end
    counter = mod(counter,switchIndex) + 1;
    DensityMatrix(jj,jj) = DensityMatrix(jj,jj)*Nuclei.State{inucleus}(stateNumber)^2;
  end
end
%DensityMatrix = DensityMatrix.^2;
if abs(trace(DensityMatrix)-1)>1e-9
  error('State is not normalized.');
end

end

%}
% ========================================================================
% Calculate signal for one orientation
% ========================================================================
function [Signal, AuxiliarySignal,Order_n_Signal,Statistics,graphs,...
  Ori_Clusters] = ...
  beginCalculateSignal(System,Method,Nuclei,Clusters,Alpha,Beta,verbose,...
  OutputData,Data,InputData,isignal)

% Assign temporary value to AuxiliarySignal
AuxiliarySignal = 'pending';

% Get rotation matrix from PDB frame to lab frame, via Euler angles
R_pdb2lab = rotateZYZ(Alpha,Beta,0);

% Rotate nuclear coordinates.
Nuclei.Coordinates = Nuclei.Coordinates*R_pdb2lab';
Statistics = Nuclei.Statistics;
Ori_Clusters = [];
if Method.Ori_cutoffs

  isotopologueStatistics = [];
  if Method.reparseNuclei
    if isignal>1
      [Nuclei, System] = parseNuclei(System, Method, Data, InputData);
    end
    isotopologueStatistics.number_1H_exchangeable = Nuclei.number_1H_exchangeable;
    isotopologueStatistics.number_1H_nonExchangeable = Nuclei.number_1H_nonExchangeable;
    isotopologueStatistics.number_2H_exchangeable = Nuclei.number_2H_exchangeable;
    isotopologueStatistics.number_2H_nonExchangeable = Nuclei.number_2H_nonExchangeable;
  end
  
  StatisticsTemporary = getPairwiseStatistics(System, Method, Nuclei);
  Nuclei.Statistics = StatisticsTemporary;
  Adjacency = getAdjacencyMatrix(System, Nuclei,Method);
  Nuclei.Adjacency = Adjacency;
  
  Statistics.isotopologueStatistics = isotopologueStatistics;
  Nuclei.Statistics = Statistics;
  
  if ~Method.cutoff.sizeDependent
    Ori_Clusters = findClusters_treeSearch(Nuclei,Method.order,...
      Method.extraOrder,{}, Method);
  end

  for clusterSize = 1:Method.order
    if Method.cutoff.sizeDependent
      Ori_Clusters = findClusters_treeSearch(Nuclei,clusterSize,...
        clusterSize,{}, Method);
    end

    % Combine arrays.
    C = [Clusters{clusterSize}; Ori_Clusters{clusterSize}];
    
    % Sort clusters
    C = sortrows(C);
    
    % Remove duplicates
    keep = [true; any(C(1:end-1,:)~=C(2:end,:),2)];
    if ~isempty(C)
      Clusters{clusterSize} = C(keep,:);
    elseif Method.emptyClusterSetsOkay
      Clusters{clusterSize} = C;
    else
      error(['Error in beginCalculateSignal(): ', ...
        'no clusters found.']);

    end
    Nuclei.numberClusters(clusterSize) = size(Clusters{clusterSize},1); 
  
  
    if Method.verbose
      fprintf('Found %d orientation clusters of size %d.\n', ...
        size(Clusters{clusterSize},1),clusterSize);
    end
  end
end

% Rotate bath spin tensors.
for inucleus = 1:Nuclei.number
  
  Nuclei.Atensor(inucleus,:) = reshape( ...
    R_pdb2lab*reshape( full(Nuclei.Atensor(inucleus,:))',3,3)*R_pdb2lab' , ...
    9,1);
  if ~System.limitToSpinHalf
    Nuclei.Qtensor(:,:,inucleus) = R_pdb2lab*Nuclei.Qtensor(:,:,inucleus)*R_pdb2lab';
    % Elementwise Qtensor manipulation used for testing.  The default filter is ones(3);
    Nuclei.Qtensor(:,:,inucleus) = Nuclei.Qtensor(:,:,inucleus).*System.nuclear_quadrupole_filter;
  end
end

% Rotate the g-matrix.
if isfield(System,'gFrame')
  % Get rotation matrix.
  g2MolRotation = rotateZYZ(-System.gFrame(1),-System.gFrame(2),-System.gFrame(3));
  
  % Rotate to molecular frame.
  System.gMatrix = g2MolRotation*System.gMatrix_gFrame*g2MolRotation';
end

% Rotate the g-matrix to the lab frame.
System.gMatrix = R_pdb2lab*System.gMatrix_gFrame*R_pdb2lab';

% Get the g-value along the magnetic field direction.
System.Electron.g = System.gMatrix(3,3);

graphs = Nuclei.Adjacency; 

if Method.useMultipleBathStates
  [Nuclei.MeanFieldCoefficients, Nuclei.MeanFieldTotal]= getMeanFieldCoefficients(Nuclei,System);
else
  Nuclei.MeanFieldCoefficients = [];
  Nuclei.MeanFieldTotal = [];
end


if isfield(Method,'exportHamiltonian') && Method.exportHamiltonian
  
  % Generate Cartesian spin-spin coupling Hamiltonian.
  
  geff=System.gMatrix(3,3);
  
  [Tensors,zeroIndex] = pairwisetensors_gpu(Nuclei.Nuclear_g, Nuclei.Coordinates,...
    1:Nuclei.number,Nuclei.Atensor, System.magneticField, System.ge, geff, System.muB, System.muN, System.mu0, System.hbar,System.theory,[]);
  
  % set file name
  H_file = [OutputData(1:end-4), '_Hamiltonian.mat'];
  
  % save Hamiltonian
  if ~isempty(Hamiltonian)
    save(H_file,'Tensors','zeroIndex');
    clear Hamiltonian;
  end
  
end

% calculate restricted signal directly
if strcmp(Method.method,'rCE')
  if verbose
    fprintf('\nCalculating restricted Cluster Expansion.\n')
  end
  [Signal,AuxiliarySignal,Order_n_Signal] = ...
    doRestrictedCE(System,Method, Nuclei, Method.r0, verbose);
  if verbose
    fprintf('\nComplete.\n')
  end
  return
elseif strcmp(Method.method,'rCCE')
  if verbose
    fprintf('\nCalculating restricted Cluster Correlation Expansion.\n')
  end
  [Signal,AuxiliarySignal,Order_n_Signal] = ...
    doRestrictedCCE(System,Method, Nuclei, verbose);
  if verbose
    fprintf('\nComplete.\n')
  end
  return
  
elseif strcmp(Method.method,'LCE')
  if verbose
    fprintf('\nCalculating Linked Cluster Expansion.\n')
  end
  [Signal,AuxiliarySignal,Order_n_Signal] = ...
    doLCE(System,Method, Nuclei,Clusters, verbose);
  if verbose
    fprintf('\nComplete.\n')
  end
  return
end

timepoints = numel(System.Time);
dt = System.dt;
linearTimeAxis = true; % Code should be changed to enforce this.

% Calculate signal
if Method.mixed_eState
  [Signal,Order_n_Signal,~] = calculateSignal_pulse(System,Method, ...
    Nuclei,Clusters,timepoints,dt,linearTimeAxis,verbose);
else
  % Calculate signal and save extra parameters (RAM intensive)
  
  % Check if the theroy is the same for every cluster size.
  if all(all(System.Theory)==any(System.Theory)) && ...
      ~strcmp(Method.method,'HD-CCE')
    
    if Method.gpu
      [Signal, AuxiliarySignal_1,AuxiliarySignal_2,...
        AuxiliarySignal_3,AuxiliarySignal_4,Signals] ...
        = calculateSignal_gpu(System, Method, Nuclei,Clusters);
      
      AuxiliarySignal = {AuxiliarySignal_1,AuxiliarySignal_2,...
        AuxiliarySignal_3,AuxiliarySignal_4};
      
      Order_n_Signal = {Signals(1,:),Signals(2,:),Signals(3,:),...
        Signals(4,:),Signals(5,:),Signals(6,:)};
    else
      
      [Signal, AuxiliarySignal, Signals] ...
        = calculateSignal(System, Method, Nuclei,Clusters);
      
      Order_n_Signal = cell(1,Method.order);
      
      for ii = 1:Method.order
        Order_n_Signal{ii} = Signals(1,ii);
      end
    end
  else
    Order_n_Signal = cell(1,Method.order);
    AuxiliarySignal = cell(1,Method.order);
    
    % Remember Method.order.
    MethodOrder = Method.order;
    
    new_order = 1;
    
    % Loop over orders.
    for iorder = 1:MethodOrder
      if iorder < new_order
        continue;
      end
      
      % Check for orders with the same theory.
      for jorder = iorder:MethodOrder
        if ~all(   all(  System.Theory(iorder:jorder,:),1  )  ...
            == any(System.Theory(iorder:jorder,:),1 ) )
          break;
        end
        new_order = jorder;
      end
      
      % Set new order to the highest order with the same theory as iorder.
      Method.order = new_order;
      
      % Calculate signals.
      if strcmp(Method.method,'HD-CCE')
        
        fractions = System.deuteriumFraction;
        
        % Change all hydrons to protons.
        System.deuteriumFraction = 0;
        Nuclei = newHydronIsotopologue(Nuclei,System);
        
        % Get proton auxiliary signals.
        [ClusterArray, Coherences_1H,Coherences_2H,...
          SubclusterIndices_2H,dimensionality,~] ...
          = calculateSignal_gpu(System, Method, Nuclei,Clusters);
        
        % Change all hydrons to deuteron.
        System.deuteriumFraction = 1;
        Nuclei = newHydronIsotopologue(Nuclei,System);

        % Get deuteron auxiliary signals.

        [~, Coherences_1D,Coherences_2D,SubclusterIndices_2D,~,~] ...
          = calculateSignal_gpu(System, Method, Nuclei,Clusters);
        
        
        [~, ...
          AuxiliarySignal_1,AuxiliarySignal_2] = ...
          doHydrogenIsotopologueCCE(...
          Coherences_1H,Coherences_2H,Coherences_1D,Coherences_2D, ...
          fractions, ClusterArray, ...
          SubclusterIndices_2H,SubclusterIndices_2D,...
          timepoints,dimensionality, Method.order,...
          Nuclei.numberClusters,Nuclei.Exchangeable,Nuclei.MoleculeID);
        
        System.deuteriumFraction = fractions;

      else
        
        if Method.gpu
          [~, AuxiliarySignal_1,AuxiliarySignal_2,...
            AuxiliarySignal_3,AuxiliarySignal_4,~] ...
            = calculateSignal_gpu(System, Method, Nuclei,Clusters);
        else
          
          [~, AuxiliarySignal_ofOrder, ~] ...
            = calculateSignal(System, Method, Nuclei,Clusters);
          
        end
      end
      % Record the appropraite signals.
      for record_order = iorder:new_order
        
        if Method.gpu
        switch record_order
          case 1
            AuxiliarySignal{record_order} = AuxiliarySignal_1;
          case 2
            AuxiliarySignal{record_order} = AuxiliarySignal_2;
          case 3
            AuxiliarySignal{record_order} = AuxiliarySignal_3;
          case 4
            AuxiliarySignal{record_order} = AuxiliarySignal_4;
        end
        else
          AuxiliarySignal{record_order} = AuxiliarySignal_ofOrder{record_order};
        end
       
          Order_n_Signal{record_order} = prod(AuxiliarySignal{record_order},1);
          if record_order > 1
            Order_n_Signal{record_order} = ...
              Order_n_Signal{record_order}.*Order_n_Signal{record_order-1};
          end
        
      
      end
      
    end
    Signal = Order_n_Signal{MethodOrder};
    
    
  end

end

end

%==========================================================================
% Set Output Data
%==========================================================================
function [OutputData,doReturn,SignalMean,experiment_time,...
  TM_powder,Order_n_SignalMean,Nuclei,statistics] = setOutput(Data,Method)
doReturn = false;
SignalMean  = [];
experiment_time = [];
TM_powder = [];
Order_n_SignalMean =[];
Nuclei = [];
% uncertainty = [];
statistics = [];
if ~isfield(Data,'overwriteLevel')
  Data.overwriteLevel = 0;
end

OutputData = Data.OutputData;
if ~isempty(OutputData)
  
  if ~strcmp(OutputData(end-3:end),'.mat')
    OutputData = [OutputData, '.mat'];
  end
  
  if isfile(OutputData)
    disp('Pre-existing save file found.')
    
    switch Data.overwriteLevel
      
      % no overwrite, return save
      case 0
        
        disp('Checking completion status.')
        try
          
          sim_ = load(OutputData);
          
          if sim_.Progress.complete && ~(isempty(sim_.uncertainty) &&  Method.getUncertainty)
            
            
            disp('Complete: returning saved data.')
            SignalMean  = sim_.SignalMean;
            experiment_time = sim_.experiment_time;
            TM_powder = sim_.TM_powder;
            % uncertainty = sim_.uncertainty;
            if isfield(sim_,'statistics')
              statistics = sim_.statistics;
            end
            
            if isfield(sim_,'Order_n_SignalMean')
              Order_n_SignalMean = sim_.Order_n_SignalMean;
            end
            if isfield(sim_,'Nuclei')
              Nuclei = sim_.Nuclei;
            end
            doReturn = true;
            return;
          else
            clear('sim_')
          end
          
        catch
        end
        
        disp('Incomplete: data will be overwritten.')
        
        
        % no overwrite, backup
      case 1
        
        newOutputFile = ['%',OutputData];
        if isfile(newOutputFile)
          disp(['Replacing ',newOutputFile, ' as ', newOutputFile, '.'])
        else
          disp(['Backing up ',OutputData, ' as ', newOutputFile, '.'])
        end
        if isunix
          command = ['mv ', OutputData, ' ', newOutputFile];
          system(command);
        else
          movefile(OutputData, newOutputFile);
        end
        % overwrite
      case 2
        disp('Preparing to overwrite save file.');
    end
  end
end
end
