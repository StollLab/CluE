% nuclear_spin_diffusion   Calculate signals for nuclear spin diffusion
%
% Inputs:
%   System is a struct with fields: magneticField, Electron, temperature, Time.
%     Electron is a struct with fields: spin, g, Coordinates, wavefunction.
%   Method is a struct with fields: method, r0, order.
%   Data is a structure containing
%    .InputData  name of pdb file
%    .OutputData name of mat file to store results
%    .saveLevel  0 (standard), 1 (more), 2 (all)

function [SignalMean, experiment_time, TM_powder,Order_n_SignalMean,Nuclei] = nuclear_spin_diffusion(System,Method,Data)

tic

% set defaults base on specified parameters and for unspecified parameters
[System, Method, Data] = setSystemDefaults(System,Method,Data);

if ~isfield(Method,'r_min')
  Method.r_min = 0.1*System.meter*1e-10; % m.
end

% time axis

if strcmp(System.experiment,'FID')
  experiment_time = System.Time;
elseif strcmp(System.experiment,'Hahn')
  experiment_time = 2*System.Time;
elseif strcmp(System.experiment,'CPMG') || strcmp(System.experiment,'CPMG-2D')
  experiment_time = 4*System.Time;% + 2*System.Time_;
elseif strcmp(System.experiment,'CPMG-const')
    experiment_time = 2*System.Time;
end

% save input data

if ~isfield(Data,'InputData') || ~ischar(Data.InputData)
  error(['Data.InputData is not a valid file identifier.']);
end
InputData = Data.InputData;
if ~isfile(InputData)
  error(['Could not find the specified input file ', InputData, '.']);
end

% Determine the output file name.
if isfield(Data,'OutputData') && ~isempty(Data.OutputData)
  OutputData = Data.OutputData;
  if ~strcmp(OutputData(end-3:end),'.mat')
    OutputData = [OutputData, '.mat'];
  end
else
  OutputData = '';
end

% Initiate progress tracking.
Progress.started = true;
Progress.complete = false;
Progress.Completed_Orders = zeros(1,Method.order);

% Save if filename is provided.
if ~isempty(OutputData)
  Input.System = System;
  Input.Method = Method;
  Input.Data = Data;
  save(OutputData,'Input','experiment_time','Progress');
end

% Set verbosity.
verbose = Method.verbose;

if verbose, fprintf('Setting up...\n');  end

if strcmp( InputData(end-3:end), '.mat') % check to see if InputData is a saved file.
  newMethod = Method.method;
  load(InputData);
  Method.method = newMethod;
  clear newMethod
  Progress.LoadSavedData = true;
  
elseif min( (InputData(end-3:end)) == '.pdb') || strcmp(InputData,'user')
  Nuclei = parseNuclei(System, Method, InputData);
  System.Electron.Coordinates = [0,0,0];
  if Nuclei.number < 1
    Signals{1} = ones(size(System.Time));
    SignalMean = Signals{1};
    fprintf(2,'\n There are too few magnetic nuclei in the system for nuclear spin diffusion.\n')
    toc
    return;
  end
  
  
else
  error('Input data not recognized')
end

Progress.DataLoaded = true;

if verbose, fprintf('Setup initialized %i nuclei.\n', Nuclei.number); end

if Nuclei.number < Method.order
  fprintf('Reducing the maximum cluster size to the system size of %d.\n',Nuclei.number)
  Method.order = double(Nuclei.number);
end

% only possible for very small systems
if strcmp(Method.method,'full')
  Method.order = Nuclei.number;
end


% ========================================================================
% Compile list of connected clusters
% ========================================================================
Clusters = [];
if ~Method.conserveMemory
  % Loop over cluster sizes, start at the largest (most time consuming) size
  for clusterSize = Method.order:-1:1
    
    if strcmp(Method.method,'full')
      if clusterSize<Nuclei.number
        Clusters{clusterSize} = [];
      else
        Clusters{Nuclei.number} = Nuclei.Index;
      end
      continue
    end
    
    if verbose, fprintf('Finding clusters of size %d.\n', clusterSize); end
    
    Clusters{clusterSize} = findClusters(Nuclei, clusterSize);
    
    % save the number of clusters of each size for export to the user
    Nuclei.numberClusters(clusterSize) = size(Clusters{clusterSize},1);
    
    if verbose
      fprintf('  Found %d clusters of size %d.\n', size(Clusters{clusterSize},1),clusterSize);
    end
    
  end
  
  if strcmp(Method.method,'count clusters')
    SignalMean = Nuclei.numberClusters;
    experiment_time = 1:length(SignalMean);
    TM_powder = [];
    Order_n_SignalMean = [];
    toc
    return
  end
  
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
  numCores = feature('numcores');
  
  % create parallel pool
  pool = parpool(numCores);
else
  pool = [];
end

%===============================================================================
% Set up powder averaging.
%===============================================================================
% Powder average settings
gridSize = System.gridSize;
if gridSize==1
  System.averaging = 'none';
end
% useless feature (to be removed)
if isfield(System,'gammaGridSize')
  gammaGridSize = System.gammaGridSize; % number of point for the gamma Euler angle.
else
  gammaGridSize = 1; % number of point for the gamma Euler angle.
end

if strcmp(System.averaging,'powder')
  
  Grid = lebedev_grid(gridSize);
  % Use the inversion symmetry of the spin-Hamiltonian to skip all the
  % orientations with z < 0.
  keep = Grid.z>=0;
  
  % Double the weight of all points whose mirror is to be skipped.
  doubleWeight = Grid.z>0; 
  Grid.w(doubleWeight) = 2*Grid.w(doubleWeight); 
  Grid.w = Grid.w(keep);
  Grid.w = Grid.w/sum(Grid.w);
  
  Grid.x = Grid.x(keep);
  Grid.y = Grid.y(keep);
  Grid.z = Grid.z(keep);
  

  Gamma = linspace(0,2*pi,gammaGridSize+1); % grid of gamma values.
  
  % Convert xyz coordinates to alpha/beta angles
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
    
elseif strcmp(System.averaging,'none')
  
  % Use only the PDB file orientation.
  gridSize = 1;
  gammaGridSize = 1;
  Alpha = 0;
  Beta = 0;
  Gamma = 0;
  gridWeight = 1;
  
elseif strcmp(System.averaging,'xy')
  
  % Average over rotations about the B0 direction.
  gridSize = 15;
  gammaGridSize = 1;
  Alpha = linspace(0,pi,gridSize+1);
  Alpha(end) = [];
  Beta = ones(1,gridSize)*pi/2;
  Gamma = 0;
  gridWeight = ones(gridSize,1)/gridSize;
  
end
gridWeight = repmat(gridWeight(:),1,gammaGridSize);

nOrientations = gridSize*gammaGridSize;

% initialize result variables
Signals{nOrientations} = [];
TM(nOrientations) = 0;
AuxiliarySignal{nOrientations} = [];
Calculate_Signal{nOrientations} = true;

% update progress
Progress.Order_n_Mean = 'pending';
Order_n_Signals{nOrientations} = [];

SignalMean = zeros(1,System.timepoints^System.dimensionality);

% initialize
for iorder = Method.order:-1:1
  if strcmp(System.experiment,'CPMG-2D')
  Order_n_SignalMean{iorder} = zeros(1,System.timepoints^2);
  else
    Order_n_SignalMean{iorder} = zeros(1,System.timepoints);
  end
end

SignalsToCalculate = [];

% Check to see if file already exists.
for isignal = 1:nOrientations
  
  % temporary file for partial saving
  temp_file = ['temp_', OutputData, '_sig_', num2str(isignal), '.mat'] ;
  
  Calculate_Signal{isignal} = true;
  SignalsToCalculate(end+1) = isignal;
  
  % check if file alread exists
  if isfile(temp_file)
    try
      % load partial save
      load(temp_file,'signal','order_n','seed');
      
      % check progress
      if seed == Method.seed && progress_powder
        
        % use loaded data
        if verbose, fprintf(['Loading signal %d from ', temp_file, '.'],isignal); end
        
        % set values to simulation variables
        Calculate_Signal{isignal} = false;
        
        if Method.sparseMemory
          SignalMean = SignalMean + signal;
          for iorder = Method.order:-1:1
            Order_n_SignalMean{iorder} = Order_n_SignalMean{iorder} + order_n{iorder};
          end
        else
          Signals{isignal} = signal;
          Order_n_Signals{isignal} = order_n;
        end
        
      end
    end
  end
  
end

numberOfSignals = length(SignalsToCalculate);

TempSignals{numberOfSignals+1} = [];
Temp_Order_n_Signals{numberOfSignals+1} = [];

Statistics = cell(numberOfSignals,1);
graphs = cell(numberOfSignals,1);

parallelComputing = Method.parallelComputing && ~Method.conserveMemory;
saveAll = Data.saveLevel==2;

if parallelComputing
  
  parfor isignal = 1:numberOfSignals
    
    [TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics{isignal},graphs{isignal}] ...
      = getOrientationSignals(System,Method,Nuclei,Clusters, Alpha,Beta, ...
      Gamma,isignal,verbose,OutputData,Progress,SignalsToCalculate,gammaGridSize,gridWeight,nOrientations);
    
    TempSignals{isignal} = TempSignals_;
    Temp_Order_n_Signals{isignal} = Temp_Order_n_Signals_;
    
    if saveAll || Method.getNuclearContributions
      AuxiliarySignal{isignal} = AuxiliarySignal_;
    end
  end
  
else
  
  for isignal = 1:numberOfSignals
    
    [TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics{isignal},graphs{isignal}] ...
      = getOrientationSignals(System,Method,Nuclei,Clusters, Alpha,Beta, ...
      Gamma,isignal,verbose,OutputData,Progress,SignalsToCalculate,gammaGridSize,gridWeight,nOrientations);
    
    TempSignals{isignal} = TempSignals_;
    Temp_Order_n_Signals{isignal} = Temp_Order_n_Signals_;
    
    if Method.sparseMemory
    
      SignalMean = SignalMean + TempSignals_;
      for iorder = Method.order:-1:1
        Order_n_SignalMean{iorder} = Order_n_SignalMean{iorder} + Temp_Order_n_Signals_{iorder};
      end
      
    else
      
      TempSignals{isignal} = TempSignals_;
      Temp_Order_n_Signals{isignal} = Temp_Order_n_Signals_;
      
      if saveAll || Method.getNuclearContributions || Method.getNuclearSpinContributions
        AuxiliarySignal{isignal} = AuxiliarySignal_;
      end
      
    end
    
  end
  
end

tempSignalIndex = 0;
if ~Method.sparseMemory 
  for isignal = 1:nOrientations
    if Calculate_Signal{isignal} 
      
      tempSignalIndex = tempSignalIndex  +1;
      Signals{isignal} = TempSignals{tempSignalIndex};
      Order_n_Signals{isignal} = Temp_Order_n_Signals{tempSignalIndex};
    end
    TM(isignal) = getTM(experiment_time,Signals{isignal});
    
  end
  
  clear('TempSignals');
end
if ~isempty(pool)
  delete(pool);
end

Nuclei.Statistics = Statistics;
Nuclei.graphs = graphs;
% ========================================================================
% End of main loop
% ========================================================================

if ~Method.sparseMemory 
  for isignal = 1:nOrientations
    SignalMean = SignalMean + Signals{isignal};
    %Signals{isignal} = abs(Signals{isignal});
    if ~ischar(Order_n_SignalMean{iorder}) && ~isempty(Order_n_Signals{1})
      for iorder = 1:Method.order
        try
          Order_n_SignalMean{iorder} = Order_n_SignalMean{iorder} + Order_n_Signals{isignal}{iorder};
        catch
          fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
          
          fprintf('Error in Order %d mean at signal number %d.\n',iorder,isignal);
          disp('Could not evaluate');
          disp('Order_n_SignalMean{iorder} = Order_n_SignalMean{iorder} + Order_n_Signals{isignal}{iorder};');
          fprintf('Order_n_SignalMean{%d} = Order_n_SignalMean{%d} + Order_n_Signals{%d}{%d};\n', ...
            iorder,iorder,isignal,iorder);
          if exist('Order_n_Signals')
            if iscell(Order_n_Signals) && length(Order_n_Signals)>=isignal
              
              if iscell(Order_n_Signals{isignal}) && length(Order_n_Signals{isignal}) >= iorder
                fprintf('Length of Order_n_Signals{%d}{%d} is %d.\n',...
                  isignal,iorder, length(Order_n_Signals{isignal}{iorder}));
                
                fprintf('Length of Order_n_SignalMean{%d} is %d.\n',...
                  iorder, length(Order_n_SignalMean{iorder})  );
                
              else
                fprintf('Order_n_Signals{%d}{%d} does not exist.\n',isignal,iorder);
              end
              
            else
              fprintf('Order_n_Signals{%d} does not exist.\n',isignal);
            end
          else
            fprintf('Order_n_Signals does not exist.\n');
          end
          
          disp('Attempting to recover data...')
          
          
          
          if iorder == 1
            index_gamma = mod(isignal-1,gammaGridSize) + 1;
            iOri = 1 + (isignal - index_gamma)/gammaGridSize;
            
            Order_n_Signals{isignal}{iorder} = gridWeight(iOri,index_gamma)*ones(size(experiment_time));
            
            disp('Recovered.');
            
          elseif iorder == Method.order
            Order_n_Signals{isignal}{iorder} = Signals{isignal};
            disp('Recovered.');
          else
            Progress.Order_n_Mean = false;
            disp('Failed.');
            Order_n_Signals{isignal}{iorder} = nan(size(experiment_time));
          end
          
          Order_n_SignalMean{iorder} = Order_n_SignalMean{iorder} + Order_n_Signals{isignal}{iorder};
          fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        end
      end
    end
    % Signals{isignal} = Signals{isignal}/max(Signals{isignal});
  end
end

% update progress
if Progress.Order_n_Mean, Progress.Order_n_Mean = true; end

if System.dimensionality ==2 && min(size(SignalMean))==1
  SignalMean = reshape(SignalMean',System.timepoints,System.timepoints)';
  
  for iorder = 1:Method.order
    Order_n_SignalMean{iorder} = reshape(Order_n_SignalMean{iorder}',System.timepoints,System.timepoints)';
  end
end

if strcmp(Method.method,'count clusters')
  SignalMean = SignalMean(1:Method.order);
  experiment_time = 1:Method.order;
  SignalMean = SignalMean/SignalMean(1)*double(Nuclei.number);
end
%SignalMean = SignalMean/max(SignalMean);
TM_powder = getTM(experiment_time,abs(SignalMean));
Progress.complete = true;
if verbose
  fprintf('Processing\n');
end

% decide what to save
if ~isempty(OutputData)
  switch Data.saveLevel
    case 0
      if Method.sparseMemory
        save(OutputData,'SignalMean','Order_n_SignalMean','TM_powder','Progress','-append');
      else
        save(OutputData,'SignalMean','Signals','TM','TM_powder','Progress','-append');
      end
    case 1
      if Method.sparseMemory
        save(OutputData,'Nuclei','Order_n_Signals','-append');
      else
        save(OutputData,'SignalMean','Signals','TM','TM_powder','Progress',...
          'Nuclei','Order_n_SignalMean','Order_n_Signals','-append');
      end
    case 2
      save(OutputData,'-v7.3');
  end
end

% delete temporary files
for ii = 1:nOrientations
  temp_file = ['temp_', OutputData, '_sig_', num2str(ii), '.mat'] ;
  if isfile(temp_file)
    delete(temp_file);
  end
end

% calculate the contributions from each spin
% placed after the save since getSpinContributions() is still buggy
%NuclearContribution.findContributions = Method.getNuclearContributions;
if Method.getNuclearSpinContributions
  getNuclearSpinContributions([OutputData,'SpinContribution.mat'], ...
    Nuclei, System, nOrientations, Clusters, Signals, AuxiliarySignal,Method, experiment_time, gridWeight, TM_powder, Input, SignalMean, Order_n_SignalMean)
end

if isfield(Method, 'getNuclearContributions') && Method.getNuclearContributions
  if isnan(TM)
    NuclearContribution.TM=TM_powder;
    fprintf('Cannot find nuclear contributions whe TM is Nan.');
  else
    NuclearContribution = getSpinContributions(System, Nuclei, Signals, AuxiliarySignal, Clusters,Method,OutputData,gridWeight,gammaGridSize, TM_powder);
    save(OutputData,'NuclearContribution','-append');
  end
else
  NuclearContribution = 'not calculated';
end
if verbose
  fprintf('\nCompleted Nuclear Spin Diffusion\n');
  fprintf('\nNuclear spin decoherence time = %d s. \n',TM_powder);
end

if false
  v_ee = eeDecoherence(System, Method, experiment_time, Nuclei);
  save(OutputData,'v_ee','-append');
  SignalMean = abs(v_ee).*SignalMean;
  disp(v_ee(end));
end

if false
  v_ee = eeDecoherence(System, Method, experiment_time, Nuclei);
  save(OutputData,'v_ee','-append');
  %   SignalMean = abs(v_ee).*SignalMean;
  SignalMean = abs(v_ee);
  
  disp(v_ee(end));
end
toc
end
% ========================================================================
% ========================================================================



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
% Calculates signal for a set of orientations
% ========================================================================
function [TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics_isignal,graphs_isignal] ...
    = getOrientationSignals(System,Method,Nuclei,Clusters, Alpha,Beta,Gamma,isignal,verbose,OutputData,Progress,SignalsToCalculate,gammaGridSize,gridWeight,iSignal_max)
     
% grid point indices
adjusted_isignal = SignalsToCalculate(isignal);
index_gamma = mod(adjusted_isignal-1,gammaGridSize) + 1;
igrid = 1 + (adjusted_isignal - index_gamma)/gammaGridSize;

% calculate coherence signal
[TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics_isignal,graphs_isignal] = ...
  calculateSignal(System,Method,Nuclei,Clusters,...
  Alpha(igrid),Beta(igrid),Gamma(index_gamma),verbose,OutputData,Progress);

normalizing_factor = TempSignals_(1);
TempSignals_ = gridWeight(igrid,index_gamma)*TempSignals_/normalizing_factor;
if strcmp(Method.method,'full')
  return
end

% Order n signals
if ~ischar(Temp_Order_n_Signals_)
  for iorder = 1:Method.order
    % set the max amplitude to the weight
    Temp_Order_n_Signals_{iorder} = gridWeight(igrid,index_gamma)*Temp_Order_n_Signals_{iorder}/normalizing_factor;
  end
end

if verbose, fprintf('\nCompleted orientation %d/%d.\n',isignal,iSignal_max); end

% Save to file.
if Method.partialSave
  temp_file = ['temp_', OutputData, '_sig_', num2str(adjusted_isignal), '.mat'] ;
  parsavefile = matfile(temp_file,'writable',true);
  parsavefile.signal = TempSignals_;
  parsavefile.order_n = Temp_Order_n_Signals_;
  parsavefile.progress_powder = true;
  parsavefile.seed = Method.seed;
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


% ========================================================================
% Calculate signal for one orientation
% ========================================================================
function [Signal, AuxiliarySignal,Order_n_Signal,Statistics,graphs] = ...
  calculateSignal(System,Method,Nuclei,Clusters,Alpha,Beta,Gamma,verbose,OutputData,Progress)

% Assign temporary value to AuxiliarySignal
AuxiliarySignal = 'pending';

% Get rotation matrix from PDB frame to lab frame, via Euler angles
R_pdb2lab = rotateZYZ(Alpha,Beta,Gamma);

% Rotate nuclear coordinates.
Nuclei.Coordinates = Nuclei.Coordinates*R_pdb2lab';

% Rotate bath spin tensors.
for inucleus = 1:Nuclei.number
  Nuclei.Atensor(:,:,inucleus) = R_pdb2lab*Nuclei.Atensor(:,:,inucleus)*R_pdb2lab';
  Nuclei.Qtensor(:,:,inucleus) = R_pdb2lab*Nuclei.Qtensor(:,:,inucleus)*R_pdb2lab';
  % Elementwise Qtensor manipulation used for testing.  The default filer is ones(3); 
  Nuclei.Qtensor(:,:,inucleus) = Nuclei.Qtensor(:,:,inucleus).*System.nuclear_quadrupole_filter;
  
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

graphs = []; % getBathGraphs(Nuclei);

if System.useMeanField
  [Nuclei.MeanFieldCoefficients, Nuclei.MeanFieldTotal]= getMeanFieldCoefficients(Nuclei,System);
else
  Nuclei.MeanFieldCoefficients = [];
  Nuclei.MeanFieldTotal = [];
end

if Method.precalculateHamiltonian
  if strcmp(Method.HamiltonianType,'pairwise')
    
    % Generate Cartesian spin-spin coupling Hamiltonian.
    Hamiltonian_pairwise = pairwiseHamiltonian(System,Nuclei,1:Nuclei.number);
    
  else
    
    % calculate compressed Hailtonian
    [Hamiltonian_diagonal,Hamiltonian_offDiagonal,zeroIndex] = constructHamiltonian(System,Nuclei,[1:Nuclei.number]);
    
  end
end

if isfield(Method,'exportHamiltonian') && Method.exportHamiltonian
  
  % Generate Cartesian spin-spin coupling Hamiltonian.
  [Hamiltonian,zeroIndex] = pairwisetensors(System,Nuclei,[1:Nuclei.number]);
  
  % set file name
  H_file = [OutputData(1:end-4), '_Hamiltonian.mat'];
  
  % save Hamiltonian
  if ~isempty(Hamiltonian)
    save(H_file,'Hamiltonian','zeroIndex');
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
    doRestrictedCCE(System,Method, Nuclei, Method.r0, verbose);
  if verbose
    fprintf('\nComplete.\n')
  end
  return
end

timepoints = size(System.Time,2);
dt = System.dt;
linearTimeAxis = true; % Code should be changed to enforce this.

% Calculate signal
if Method.conserveMemory
  if Method.vectorized
    
    [Signal, Signals,Statistics.ClusterCount, Statistics.NonClusterCount, Statistics.total_clusters, Statistics.total_nonclusters] ...
      = calculateSignal_conserveMemory_gpu(System, Method, Nuclei, timepoints,dt,OutputData,Progress);
    
    %     Order_n_Signal = {Signals(1,:),Signals(2,:),Signals(3,:),Signals(4,:)};
    Order_n_Signal = cell(1,Method.order);
    for ii=1:Method.order
      Order_n_Signal{ii} = Signals(ii,:);
    end
    
  else
    [Signal, Order_n_Signal,Statistics] = calculateSignal_conserveMemory(System, Method, Nuclei, timepoints,dt, linearTimeAxis,verbose,OutputData,Progress);
  end
  
elseif Method.mixed_eState
  [Signal, Order_n_Signal,Statistics] = calculateSignal_pulse(System, Method, Nuclei,Clusters, timepoints,dt, linearTimeAxis,verbose);
else
  % calculate signal and save extra parameters (RAM intensive)
  
  if Method.vectorized
    [Signal, AuxiliarySignal_1,AuxiliarySignal_2,AuxiliarySignal_3,AuxiliarySignal_4,Signals] ... 
       = calculateSignal_gpu0(System, Method, Nuclei,Clusters);
    AuxiliarySignal = {AuxiliarySignal_1,AuxiliarySignal_2,AuxiliarySignal_3,AuxiliarySignal_4};
    
    Statistics = [];
    Order_n_Signal = {Signals(1,:),Signals(2,:),Signals(3,:),Signals(4,:),Signals(5,:),Signals(6,:)};
  else
    [Signal,AuxiliarySignal,Order_n_Signal,Statistics] = calculateSignal_standard(System, Method, Nuclei,Clusters, timepoints,dt, linearTimeAxis,verbose);
  end
  
end

end


% ========================================================================
% New Function
% ========================================================================
function Coherences = propagate(System, Method, DensityMatrix,timepoints,dt,Hamiltonian_beta,Hamiltonian_alpha,isotopeProbability, linearTimeAxis,verbose)
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

% ========================================================================
% Calculates signal without GPU or parallel capabilities
% ========================================================================

function [Signal, AuxiliarySignal, Order_n_Signal,Statistics] = calculateSignal_standard(System, Method, Nuclei,Clusters, timepoints,dt, linearTimeAxis,verbose)
verboseCounter = 0;
verboseThreshold = 0.05*size(Clusters{Method.order_lower_bound},1);

for clusterSize = Method.order:-1:1
  
  % Find coherences
  numClusters = size(Clusters{clusterSize},1);
  for iCluster = numClusters:-1:1
    if strcmp(Method.method,'full') && clusterSize~=Method.order
      Coherences{clusterSize,iCluster} = 1;
      continue
    end
    
    if Clusters{clusterSize}(iCluster,1)==0
      Coherences{clusterSize,iCluster} = 1;
      continue
    end
    
    % The isotope probability is used to determine the likelihood a
    % particular cluster will the assumed spins.
    isotopeProbability = prod(Nuclei.Abundance(Clusters{clusterSize}(iCluster,:)));
    
    if clusterSize > 1 && ~strcmp(Method.method,'full')
      SubclusterIndices{clusterSize,iCluster} = findSubclusters(Clusters,clusterSize,iCluster);
    elseif clusterSize == 1
      SubclusterIndices{clusterSize,iCluster} = [];
    end
        
    % Get cluster Hamiltonian.
    if ~Method.precalculateHamiltonian && strcmp(Method.HamiltonianType,'pairwise')
      
      % Generate Cartesian spin-spin coupling Hamiltonian.
      [Hamiltonian_pairwise,zeroIndex] = pairwiseHamiltonian(System,Nuclei,Clusters{clusterSize}(iCluster,:));
      % Generate spin Hamiltonian.
      [ClusterHamiltonian{1},ClusterHamiltonian{2}] = ...
        assembleHamiltonian(Hamiltonian_pairwise,Clusters{clusterSize}(iCluster,:),System,Nuclei,zeroIndex,clusterSize);
      
    elseif strcmp(Method.HamiltonianType,'pairwise')
      
      % Generate Cartesian spin-spin coupling Hamiltonian.
      %ClusterHamiltonian{eState} = assemblePairwiseHamiltonian(Hamiltonian_pairwise,Clusters{clusterSize}(iCluster,:),System, eState,Nuclei,Method.HamiltonianType);
      zeroIndex =  Clusters{clusterSize}(iCluster,1) - 1;
        [ClusterHamiltonian{1},ClusterHamiltonian{2}] = ...
          assembleHamiltonian(Hamiltonian_pairwise,Clusters{clusterSize}(iCluster,:),System,Nuclei,zeroIndex,clusterSize);
      
    elseif ~Method.precalculateHamiltonian
      
      % Generate a diassembled spin Hamiltonian.
      [Hamiltonian_diagonal,Hamiltonian_offDiagonal,zeroIndex] = constructHamiltonian(System,Nuclei,Clusters{clusterSize}(iCluster,:));
      % Generate spin Hamiltonian.
      for eState = (2*System.Electron.spin+1):-1:1
        ClusterHamiltonian{eState} = getClusterHamiltonian(Hamiltonian_diagonal,Hamiltonian_offDiagonal, eState, Clusters{clusterSize}(iCluster,:), Nuclei,zeroIndex);
      end
      
    else
      
      % Generate spin Hamiltonian.
      for eState = (2*System.Electron.spin+1):-1:1
        ClusterHamiltonian{eState} = getClusterHamiltonian(Hamiltonian_diagonal,Hamiltonian_offDiagonal, eState, Clusters{clusterSize}(iCluster,:), Nuclei,zeroIndex);
      end
      
    end
    
    % Get density matrix
    DensityMatrix = getDensityMatrix(Nuclei,Clusters{clusterSize}(iCluster,:));
    
    % Get cluster coherences
    Coherences{clusterSize,iCluster} = propagate(System, Method, DensityMatrix,timepoints,dt,ClusterHamiltonian{1},ClusterHamiltonian{2},isotopeProbability,linearTimeAxis,verbose);
    
    if verbose
      verboseCounter = verboseCounter + 1;
      if verboseCounter>=verboseThreshold
        verboseCounter = 0;
        fprintf('cluster size %d: %d/%d, (%s).\n', clusterSize,(1+size(Clusters{clusterSize},1)-iCluster),size(Clusters{clusterSize},1),datetime);
      end
    end
    
  end
  
  Statistics.total_clusters = Nuclei.numberStartSpins;
  Statistics.ClusterCount = zeros(1,Method.order);
  Statistics.ClusterCount(1) = Nuclei.numberStartSpins;
  for iorder = 2:Method.order
    Statistics.ClusterCount(iorder) = size(Clusters{iorder},1 );
    Statistics.total_clusters = Statistics.total_clusters + Statistics.ClusterCount(iorder);
    if verbose
      fprintf('Cluster size %d: %d found.\n', iorder,Statistics.ClusterCount(iorder));
    end
  end
  
  if verbose
    disp('-------------------------------------------------');
    for iorder = 1:Method.order
      fprintf('Cluster size %d: %d found.\n', iorder,Statistics.ClusterCount(iorder));
    end
    fprintf('Total clusters: %d found.\n',Statistics.total_clusters );
    disp('-------------------------------------------------');
  end
  
end


% Calculate signal
%-------------------------------------------------------------------------------
if verbose, disp('Calculating signal...'); end

switch Method.method
  case 'CE'
    % Use the cluster expansion method.
    if verbose, fprintf('\nCalculating Cluster Expansion.\n'); end
    [Signal, AuxiliarySignal, Order_n_Signal] = doClusterExpansion(Coherences,Clusters, SubclusterIndices,Method.order);
    
  case 'CCE'
    % Use the cluster correlation expansion method.
    if verbose, fprintf('\nCalculating Cluster Correlation Expansion.\n'); end
    [Signal, AuxiliarySignal, Order_n_Signal] = doClusterCorrelationExpansion(Coherences,Clusters,SubclusterIndices, Method.order,Nuclei);
    
    if Method.order_lower_bound  > 1
      Signal = Order_n_Signal{Method.order}./Order_n_Signal{Method.order_lower_bound-1};
    end
    
  case 'full'
    % Calculate the full answer.
    Signal = Coherences{Method.order,1};
    AuxiliarySignal = [];
    Order_n_Signal = [];
    % place holder only since nth order signal is not meaningful
    Order_n_Signal{1} = ones(size(Signal));
end

end

% ========================================================================
% New Function
% ========================================================================
function [Signal, Order_n_Signal,Statistics] = calculateSignal_conserveMemory(System, Method, Nuclei, timepoints,dt, linearTimeAxis,verbose,OutputData,Progress)

% initialize output
Signal = ones(size(System.Time));
for initialize_clusterSize = Method.order:-1:1
  Order_n_Signal{initialize_clusterSize} = ones(size(System.Time));
end

% Determine number of divisions to break calculation into.
if ~ischar(Method.divisions)
  numDivisions = Method.divisions;
elseif strcmp(Method.divisions,'numSpins')
  
  % maximum allowed by the current algorithm
  numDivisions = double(Nuclei.numberStartSpins);
  
elseif strcmp(Method.divisions,'numCores')
  
  % minimum that will use all nodes
  numDivisions = feature('numcores');
  
else
  numDivisions = round(sqrt(double(Nuclei.numberStartSpins)));
end


% Initialize internal records.
Cluster_Statistics = zeros(Method.order,numDivisions);
NonCluster_Statistics = zeros(Method.order,numDivisions);
Possible_Clusters = zeros(Method.order,1);
for ksize=1:Method.order
      Possible_Clusters(ksize) = nchoosek(uint64(Nuclei.number),ksize); 
end
Method_ = Method;

% Get permutation method.
shuffle = 1:Nuclei.number;
% Shuffle = zeros(Nuclei.number);
if Method.shuffle
  
  % Set seed and rng method.
  rng(Method.seed, 'twister');
  
  % A constant seed ensures restarting will not give a new permutation.
  shuffle = randperm(Nuclei.number);
end

  
% Loop over order range.
for iorder = Method.order_lower_bound:Method.order%:-1:2
  
  Method_.order = iorder;
  
  % Check if order is valid: n-CCE requires at least n nuclear spins.
  if Nuclei.number - Nuclei.startSpin + 1 < iorder
    continue;
  end
  
  % Check pigeonhole principle: each division must be unique.
  if numDivisions + iorder -1 > Nuclei.number -Nuclei.startSpin +1
    numDivisions = Nuclei.number -Nuclei.startSpin +1 - iorder + 1;
  end
  
  
  % Separate the clusters into bundles for each division.
  Bundle = uint32(ones(numDivisions,2));
  Division_Cluster.Possible_Clusters = Possible_Clusters(iorder); 
  if Method.MonteCarlo.use
  Division_Cluster.Limit = ceil(Method.MonteCarlo.Cluster_Limit(iorder)/numDivisions);
  Division_Cluster.Increment= ceil(Method.MonteCarlo.Increment(iorder)/numDivisions);
  Division_Cluster.Fraction = Method.MonteCarlo.Fraction(iorder)/numDivisions;
  else
    Division_Cluster.Limit = inf;
    Division_Cluster.Increment = inf;
    Division_Cluster.Fraction = 1;
  end
  if strcmp(Method.divisions,'numSpins')
    
    % Assign to the nth division all clusters where the index of the lowest
    % cluster is n.
    
    % Assign bundles based while ensuring that there are at least
    % n nuclei with number of nuclei > indices >= n, for n-CCE.
    if Nuclei.endSpin + iorder - 1 > Nuclei.number
      Bundle(:,1) = Nuclei.startSpin:Nuclei.number -iorder + 1;
    else
      Bundle(:,1) = Nuclei.startSpin:Nuclei.endSpin;
    end
    Bundle(:,2) = Bundle(:,1);
    
  else
    
    % For N >> k, to approximately divide N choose k into c even parts,
    % the nth part should get all clusters whose lowest index in the ranges
    % N*[  ( (n-1)/c )^1/k, (n/c)^1/k ].
    
    % Assign the first bundle.
    Bundle(1,1) = 1;
    Bundle(1,2) =  Nuclei.startSpin ...
      + round( ...
      (Nuclei.number- Nuclei.startSpin + 1) ...
      *( 1 - ( 1 - 1/numDivisions)^(1/iorder) ) ...
      *(Nuclei.numberStartSpins/Nuclei.number) ...
      );
    
    % Assign the remaining bundles.
    for iCore = 2:numDivisions
      Bundle(iCore,1) = Bundle(iCore-1,2) + 1;
      Bundle(iCore,2) =  Nuclei.startSpin ...
        + round( ...
        (Nuclei.number- Nuclei.startSpin + 1) ...
        *(1 - (1 - iCore/numDivisions)^(1/iorder)  ) ...
        *(Nuclei.numberStartSpins/Nuclei.number) ...
        );
      
      Bundle(iCore,2) = max(Bundle(iCore,2) , Bundle(iCore,1));
    end
    
    Bundle(numDivisions,2) = Nuclei.endSpin;
  end
  
  % Generate unshuffled starting clusters.
  par_cluster = cell(1,numDivisions);
  for iCore = 1:numDivisions
    par_cluster{iCore} = Bundle(iCore,1):(Bundle(iCore,1) + iorder -1);
  end
  
  % initialize division output
  partial_signal = cell(1,numDivisions);
   
  if Method.parallelComputing
    parfor iCore = 1:numDivisions % parfor
      [partial_signal{iCore},Cluster_Statistics(iorder,iCore),NonCluster_Statistics(iorder,iCore)] = conserve_memory_for_loop(System,Method,Method_,OutputData,Nuclei,par_cluster{iCore},Bundle,iorder,iCore,Division_Cluster,shuffle, linearTimeAxis,numDivisions);
    end
  else
    for iCore = 1:numDivisions % parfor
      [partial_signal{iCore},Cluster_Statistics(iorder,iCore),NonCluster_Statistics(iorder,iCore)] = conserve_memory_for_loop(System,Method,Method_,OutputData,Nuclei,par_cluster{iCore},Bundle,iorder,iCore,Division_Cluster,shuffle, linearTimeAxis,numDivisions);
    end
  end

  % Combine node outputs to the final output.
  
  Combined_AuxiliarySignal = ones(size(System.Time));
  for iCore = numDivisions:-1:1
    
    % main signal
    Combined_AuxiliarySignal = Combined_AuxiliarySignal.*partial_signal{iCore};
    
  end
  
  
  Statistics.ClusterCount(iorder) = sum(Cluster_Statistics(iorder,:) );
  Statistics.NonClusterCount(iorder) = sum(NonCluster_Statistics(iorder,:) );
  if Method.MonteCarlo.use
    clusterFraction = (Statistics.NonClusterCount(iorder) + Statistics.ClusterCount(iorder))/Possible_Clusters(iorder);
    Combined_AuxiliarySignal = Combined_AuxiliarySignal.^(1/clusterFraction);
  end
  
  Signal = Signal.*Combined_AuxiliarySignal;
  Order_n_Signal{iorder} = Signal;
  
  % Save to file.
  Progress.order = [num2str(iorder) '/' num2str(Method.order) '-CCE'];
  Progress.Completed_Orders(iorder) = true;
  save(OutputData,'Signal','Order_n_Signal','Progress','-append');
  
  if Method.clear_partialSave
    for iCore = 1:numDivisions
      partial_file = ['partial_', OutputData(1:end-4), '_',num2str(iorder),'cce_', num2str(iCore), '.mat'] ;
      if isfile(partial_file)
        delete(partial_file);
      end
    end
  end
  
end


if strcmp(Method.method,'count clusters')
  Signal = zeros(1,System.timepoints);
  Signal(1) = Nuclei.number;
  for iorder = 2:Method.order
    Signal(iorder) = sum(Cluster_Statistics(iorder,:) );
  end
end

Statistics.total_clusters = Nuclei.numberStartSpins;
Statistics.total_nonclusters = 0;
  
for iorder = 2:Method.order
  %     Statistics.ClusterCount(iorder) = sum(Cluster_Statistics(iorder,:) );
  %     Statistics.NonClusterCount(iorder) = sum(NonCluster_Statistics(iorder,:) );
  
  Statistics.total_clusters = Statistics.total_clusters + Statistics.ClusterCount(iorder);
  Statistics.total_nonclusters = Statistics.total_nonclusters + Statistics.NonClusterCount(iorder);
  
  fprintf('Cluster size %d: %d found.\n', iorder,Statistics.ClusterCount(iorder));
  fprintf('Cluster size %d: %d skipped.\n', iorder,Statistics.NonClusterCount(iorder));
end
  
if verbose
  disp('-------------------------------------------------');
  for iorder = 1:Method.order
    fprintf('Cluster size %d: %d found.\n', iorder,Statistics.ClusterCount(iorder));
  end
  fprintf('Total clusters: %d found.\n',Statistics.total_clusters );
  disp('-------------------------------------------------');
end

end

% ========================================================================
% New Function
% ========================================================================

function [partial_signal,Cluster_Statistics,NonCluster_Statistics] = conserve_memory_for_loop(System,Method,Method_,OutputData,Nuclei,par_cluster,Bundle,iorder,iCore,Division_Cluster,shuffle, linearTimeAxis,numDivisions)

% set signal to unperturbed value
partial_signal = ones(size(System.Time));
partial_signal0 = partial_signal;

Cluster_Statistics = 0;
cluster0 = 0;
NonCluster_Statistics = 0;

% Save to file.
if Method.partialSave
  partial_file = ['partial_', OutputData(1:end-4), '_',num2str(iorder),'cce_', num2str(iCore), '.mat'] ;
  if isfile(partial_file)
    try
      parload = load(partial_file,'parsignal','seed');
      if parload.seed==Method.seed
        partial_signal = parload.parsignal;
        return;
      end
    catch
    end
  end
end

% loop over all the clusters in the bundle
iBundle = Bundle(iCore,1);
bundle_counter = 0;

while iBundle <= Bundle(iCore,2)
  
  % Do not calculate more than Division_Cluster.Limit valid clusters.
  if Cluster_Statistics>Division_Cluster.Limit
    
    clusterFraction = (NonCluster_Statistics + Cluster_Statistics)/Division_Cluster.Possible_Clusters;
    partial_signal1 = partial_signal.^(1/clusterFraction);
    
    delta_sig = (partial_signal1-partial_signal0);
    rmsd = sqrt(  delta_sig*delta_sig'/System.timepoints  );
    
    
    if Method.verbose
      fprintf('cluster %d to cluster %d rmsd  = %d, (%s).\n', ...
        Division_Cluster.Limit,  cluster0, rmsd,datetime);
    end
    
    if rmsd < Method.MonteCarlo.Threshold(iorder) || clusterFraction > Division_Cluster.Fraction  
      break;
    end
    
    partial_signal0 = partial_signal1;
    Division_Cluster.Limit = Division_Cluster.Limit + Division_Cluster.Increment;
    cluster0 = Cluster_Statistics;

  end
  
  % Update location index.
  iBundle = par_cluster(1);
  
  % Shuffle nuclei.
  shuffled_cluster = shuffle(par_cluster);
    
  isvalid_1 = validateCluster(shuffled_cluster,Nuclei.ValidPair,Nuclei.graphCriterion);
  
  if ~isvalid_1
    
    % update internal records.
    NonCluster_Statistics = NonCluster_Statistics + 1;
    
    % get next cluster
    par_cluster = getNextCluster(par_cluster,Nuclei.number);
    
    % check if there are no more clusters
    if isempty(par_cluster)
      break;
    end
    
    % update location index
    iBundle = par_cluster(1);
    
    % start over
    continue;
  end
  
  % update internal records.
  Cluster_Statistics = Cluster_Statistics + 1;
  
  % index location check
  if iBundle > bundle_counter
    
    % update bundle counter
    bundle_counter = iBundle;
    
    if Method.verbose
      fprintf('division %d/%d, cluster %d, (%s).\n',iCore, numDivisions, Cluster_Statistics,datetime);
      %           disp(par_cluster);
      disp(shuffled_cluster);
    end
    
  end
  
  if strcmp(Method.method,'count clusters')
    % get next cluster
    par_cluster = getNextCluster(par_cluster,Nuclei.number);
    
    % checkk if next cluster exists
    if isempty(par_cluster)
      break;
    end
    
    % update location index
    iBundle = par_cluster(1);
    continue;
  end
  
  % re-index
  new_labels = zeros(1,max( shuffled_cluster  ));
  old_labels = zeros(1,max( shuffled_cluster  ));
  for ii=1:iorder
    new_labels( shuffled_cluster (ii) ) = ii;
    old_labels(new_labels( shuffled_cluster (ii) )) = shuffled_cluster (ii);
  end
  
  % get reduced system
  Cluster_ = new_labels(shuffled_cluster );
  
  
  % Reduced_Clusters{iorder} = new system with indices 1:iorder.
  Reduced_Clusters = getClusterSet(Cluster_);
  
  % Loop through subcluster sizes.
  for jSize = 2:iorder-1
    
    jSizeLimit = size(Reduced_Clusters{jSize},1);
    
    % Loop through subclusters.
    for jCluster = 1:jSizeLimit
      
%       cluster_ = Reduced_Clusters(jCluster,1:jsize,jsize);
      
      if ~validateCluster(shuffled_cluster(Reduced_Clusters{jSize}(jCluster,:)),Nuclei.ValidPair,Nuclei.graphCriterion)
        Reduced_Clusters{jSize}(jCluster,:) = zeros(size(Reduced_Clusters{jSize}(jCluster,:) ) );
      end
    end
  end
  
  % reduce system to only include cluster nuclei
  Reduced_SubclusterIndices = cell(iorder,1);
  for jSize = iorder:-1:1
    for jCluster = 1:size(Reduced_Clusters{jSize},1)
      % SubclusterIndices{clusterSize,iCluster} = findSubclusters(Clusters,clusterSize,iCluster);
      % findSubclusters gives Indices{size} = list of all jCluster such that Clusters{size}(jCluster,:) is a subcluster of Clusters{clusterSize}(iCluster,:)
      %
      % subC_index = SubclusterIndices{hostClusterSize,hostClusterIndex}{subclusterSize}(index)
      % --> Clusters{subclusterSize}(subC_index,:) is a subset of
      % Clusters{hostClusterSize}(hostClusterIndex,:)
      
      % find subclusters
      if Reduced_Clusters{jSize}(jCluster,1)~=0
        Reduced_SubclusterIndices{jSize,jCluster} = findSubclusters(Reduced_Clusters,jSize,jCluster);
      else
        Reduced_SubclusterIndices{jSize,jCluster} = [];
      end
      
    end
  end
  
  
  
  
  Clusters = Reduced_Clusters;
  for jSize = 1:iorder
    for jCluster = 1:size(Reduced_Clusters{jSize},1)
      if Reduced_Clusters{jSize}(jCluster,1) ~= 0
        Clusters{jSize}(jCluster,:) = shuffled_cluster(Reduced_Clusters{jSize}(jCluster,:));
      end
    end
  end
  [~, Cluster_AuxiliarySignal, ~] = calculateSignal_standard(System, Method_, Nuclei,Clusters, System.timepoints,System.dt, linearTimeAxis,false);
  
  % collect cluster contribution to a the node output
  isotopeProbability = prod(Nuclei.Abundance(Clusters{iorder}(1,:) ));
  v_ = 1 + isotopeProbability*(Cluster_AuxiliarySignal{iorder,1}-1);

  % partial_signal  = partial_signal.*Cluster_AuxiliarySignal{iorder,1};
  partial_signal  = partial_signal.*v_;
  
  % get next cluster
  par_cluster = getNextCluster(par_cluster,Nuclei.number);
  
  % checkk if next cluster exists
  if isempty(par_cluster)
    break;
  end
  
  % update location index
  iBundle = par_cluster(1);
end


% Save to file.
if Method.partialSave
  partial_file = ['partial_', OutputData(1:end-4), '_',num2str(iorder),'cce_', num2str(iCore), '.mat'] ;
  parsavefile = matfile(partial_file,'writable',true);
  parsavefile.parsignal = partial_signal;
  parsavefile.seed = Method.seed;
end


end




% ========================================================================
% New Function
% ========================================================================

% Find all subsets of a sequential set, {1,2,...n}.
function Clusters = getClusterSet(Cluster)

% Determine cluster size.
numberSpins = length(Cluster);

% Loop over all proper sub-cluster sizes.
for icluster = numberSpins:-1:1
  
  % Determine number of sub-clusters
  n_clusters = nchoosek(numberSpins,icluster);
  
  % Initialize Clusters
  Clusters{icluster} = zeros(n_clusters,icluster);
  
  % Set the first subset.
  Clusters{icluster}(1,:) = Cluster(1:icluster);
  
  % Set the remaining subsets.
  for ii = 2:n_clusters
    Clusters{icluster}(ii,:) = getNextCluster(Clusters{icluster}(ii-1,:),numberSpins);
  end
  
end
end
