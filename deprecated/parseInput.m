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
