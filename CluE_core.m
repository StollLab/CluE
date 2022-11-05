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
  TM_powder,Order_n_SignalMean,Nuclei,statistics] = CluE_core(System,Method,Data)

tic

% Set defaults base on specified parameters and for unspecified parameters
[System, Method, Data,statistics] = setDefaults(System,Method,Data);

[OutputData,doReturn,SignalMean,experiment_time,...
  TM_powder,Order_n_SignalMean,Nuclei,uncertainty] = setOutput(Data,Method);
if doReturn
  toc
  return
end

[System,Method,Data,Nuclei,Progress,experiment_time] =...
  parseInput(System,Method,Data,OutputData);

% ========================================================================
% Compile list of connected clusters
% ========================================================================
[Clusters,Nuclei] = getClusters(Nuclei,Method,Data);

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



%===============================================================================
% Set up powder averaging.
%===============================================================================
% Powder average settings
[Alpha, Beta, gridWeight,GridInfo] = getOrientationGrid(System);


%===============================================================================
% Calculate average signal.
%===============================================================================
[SignalMean, experiment_time, ...
  TM_powder,Order_n_SignalMean,Nuclei,statistics] ...
  = averageSignals( ...
  System,Method,Data,OutputData,...
  experiment_time, Nuclei,Clusters,...
  Alpha,Beta,gridWeight,GridInfo,...
  uncertainty,Progress);

toc
end
% ========================================================================
% ========================================================================



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


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [System,Method,Data,Nuclei,Progress,experiment_time] = ...
  parseInput(System,Method,Data,OutputData)

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
Data.InputData = Data.InputData;
if ~isfile(Data.InputData) 
  if System.RandomEnsemble.include
    Data.InputData = 'System.RandomEnsemble.include';
  else
    error(['Could not find the specified input file ', Data.InputData, '.']);
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

if strcmp(Data.InputData(end-3:end),'.mat') % check to see if InputData is a saved file.
  newMethod = Method.method;
  load(Data.InputData);
  Method.method = newMethod;
  clear newMethod
  Progress.LoadSavedData = true;
  
elseif min( (Data.InputData(end-3:end)) == '.pdb') || strcmp(Data.InputData,'user') ...
    || strcmp(Data.InputData,'System.RandomEnsemble.include')
  
  if Method.useCentralSpinSystem
    System.pdb = parsePDBfile(Data.InputData, System.angstrom);


    if ~isempty(System.pdbTranslation)
      System.pdb = translatePDB(System.pdb,System);
    end

    if System.pdbRotate
      System.pdb = rotatePDB(System.pdb,System);
    end

  end
  [Nuclei, System] = parseNuclei(System, Method, Data, Data.InputData);
  
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

if verbose
  fprintf('Setup initialized system with %i nuclei.\n', Nuclei.number); 
end

if Nuclei.number < Method.order
  fprintf('Reducing the maximum cluster size to the system size of %d.\n',...
    Nuclei.number)
  Method.order = double(Nuclei.number);
end

% Only possible for very small systems
if strcmp(Method.method,'full')
  Method.order = Nuclei.number;
end

end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>