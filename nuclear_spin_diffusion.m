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

function [SignalMean, experiment_time, TM_powder,Order_n_SignalMean] = nuclear_spin_diffusion(System,Method,Data)

tic

% set defaults base on specified parameters and for unspecified parameters
[System, Method, Data] = setDefaults(System,Method,Data);

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
  save(OutputData,'Input','experiment_time');
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
  if Nuclei.number < 2
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
  %  warning('Parallel Computing is disabled.');
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

parallelComputing = Method.parallelComputing && ~Method.conserveMemory;
saveAll = Data.saveLevel==2;

if parallelComputing
  
  parfor isignal = 1:numberOfSignals
    
    [TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics{isignal}] ...
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
    
    [TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics{isignal}] ...
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
      
      if saveAll || Method.getNuclearContributions
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
      save(OutputData);
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
function [TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics_isignal] ...
    = getOrientationSignals(System,Method,Nuclei,Clusters, Alpha,Beta,Gamma,isignal,verbose,OutputData,Progress,SignalsToCalculate,gammaGridSize,gridWeight,iSignal_max)
     
% grid point indices
adjusted_isignal = SignalsToCalculate(isignal);
index_gamma = mod(adjusted_isignal-1,gammaGridSize) + 1;
igrid = 1 + (adjusted_isignal - index_gamma)/gammaGridSize;

% calculate coherence signal
[TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics_isignal] = ...
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
function [Signal, AuxiliarySignal,Order_n_Signal,Statistics] = ...
  calculateSignal(System,Method,Nuclei,Clusters,Alpha,Beta,Gamma,verbose,OutputData,Progress)

% Assign temporary value to AuxiliarySignal
AuxiliarySignal = 'pending';

% Get rotation matrix from PDB frame to lab frame, via Euler angles
R_pdb2lab = rotateZYZ(Alpha,Beta,Gamma);

% Rotate nuclear coordinates.
Nuclei.Coordinates = Nuclei.Coordinates*R_pdb2lab';

% Rotate nuclear quadrupole tensors.
for inucleus = 1:Nuclei.number  
  Nuclei.Qtensor(:,:,inucleus) = R_pdb2lab*Nuclei.Qtensor(:,:,inucleus)*R_pdb2lab';
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
dt = System.Time(2)-System.Time(1);
t0 = System.t0;
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
       = calculateSignal_gpu0(System, Method, Nuclei,Clusters, timepoints,dt,t0);
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
% Find decay time (where signal drops to 1/e of initial amplitude
% ========================================================================
function TM = getTM(t,Signal)

% normalize signal maximum to 1
Signal = abs(Signal/Signal(1));

if any(isnan(Signal))
  TM = inf;
  return
end

e_inv = exp(-1);
if min(Signal) > e_inv
  TM = inf;
  return
end

% translate signal to put TM at the x-intercept
searchVec = abs(Signal) - e_inv;
negatives = find(searchVec<=0);
if isempty(negatives)
  TM = inf;
  return
end

% initial search index
iTM = negatives(1);
if length(iTM)>1
  iTM = iTM(1);
end
% number of time points
N = length(t);

% initial TM guess
TM = t(iTM);

% search
searching = true;
while searching
  
  % get TM guess
  v1 = Signal(iTM);
  i_mismatch = 0;
  
  % examine guess
  if v1 == e_inv
    
    % TM found
    searching = false;
    break;
    
  elseif v1 > e_inv
    
    % TM guess too soon
    % make second guess
    jTM = iTM + 1;
    i_mismatch = +1;
    
  elseif v1 < e_inv
    
    % TM guess too late
    % make second guess
    jTM = iTM - 1;
    i_mismatch = -1;
    
  end
  
  % checkguess validity
  if (jTM<1)
    
    % guess outside range of data
    TM = nan;
    return;
    
  end
  if (jTM>N)
    
    % guess outside range of data
    TM = inf;
    return;
    
  end
  
  % new guess
  v2 = Signal(jTM);
  jmismatch = 0;
  
  % examine new guess
  if v2 == e_inv
    
    % TM found
    TM = t(jTM);
    iTM = jTM;
    v1 = v2;
    searching = false;
    break;
  elseif v1 > e_inv
    
    % TM guess too soon
    jmismatch = +1;
    
  elseif v2 < e_inv
    
    % TM guess too late
    jmismatch = -1;
    
  end
  
  % error parameter
  mismatch = i_mismatch*jmismatch;
  
  %decide how to make another guess
  if mismatch > 0
    
    % shift over one index
    TM = t(jTM);
    iTM = jTM;
    v1=v2;
    
  else
    
    % one guess is too late, the other too soon
    searching = false;
    
    % interpolate TM between the two guesses
    d2tau = t(jTM) - t(iTM);
    dv   =  v2-v1;
    TM = (d2tau/dv)*(e_inv - v1 + t(iTM)*dv/d2tau);
    
  end
  
end
end

% ========================================================================
% New Function
% ========================================================================

function [System,Method, Data] = setDefaults(System,Method, Data)

if ~isfield(Method, 'sparseMemory')
  Method.sparseMemory = false;
end
if ~isfield(Method, 'conserveMemory')
  Method.conserveMemory = false;
end

if Method.conserveMemory
  % set the Method to Memory conserving mode
  
  % Auxiliary signals are not saved so getNuclearContributions will fail.
  Method.getNuclearContributions = false;
  
  % The following behaviors are are used regardless of settings,
  % so the settings are modified to reflect what is done.
  
  Method.precalculateHamiltonian = false;
  Method.HamiltonianType = 'pairwise';
  
end
if ~isfield(Method,'gpu')
  Method.gpu = false;
end
if ~isfield(Method,'vectorized')
  Method.vectorized = true;
end
if ~isfield(Method,'exportHamiltonian')
  Method.exportHamiltonian = false;
end
if ~isfield(Method,'propagationDomain')
  Method.propagationDomain='time-domain'; % fastest method
end
% Toggle for saving each orientation,
if ~isfield(Method,'partialSave')
  Method.partialSave = true;
end
if ~isfield(Method,'clear_partialSave')
  Method.clear_partialSave = true;
end

if ~isfield(Method,'parallelComputing')
  Method.parallelComputing = false;
end

% Set verbosity.
if ~isfield(Method,'verbose')
  Method.verbose = false;
end

% cutoff criteria
if ~isfield(Method,'Criteria') || isempty(Method.Criteria)
  Method.Criteria = {'neighbor'};
end
if ~isfield(Method,'cutoff') 
  Method.cutoff.modulation = 0;
  Method.cutoff.dipole = 0;
  Method.cutoff.max_distance = inf;
  Method.cutoff.min_distance = 0;
  Method.cutoff.hyperfine_sup = inf;
  Method.cutoff.hyperfine_inf = 0;
end

% if ~isfield(Method,'dipole_cutoff')
%   Method.dipole_cutoff = 0;
% end

% cluster order max
if ~isfield(Method,'order')
  Method.order = 2;
end

% cluster order min
if ~isfield(Method,'order_lower_bound')
  Method.order_lower_bound = 1;
end
if ~isfield(Method,'shuffle')
  Method.shuffle = true;
end

Method.order_lower_bound = max(1,floor(Method.order_lower_bound));
Method.order_lower_bound= min(Method.order,Method.order_lower_bound);

if ~isfield(Method,'graphCriterion')
  Method.graphCriterion = 'connected';
end

% Monte Carlo option
if ~isfield(Method,'MonteCarlo')
    Method.MonteCarlo.use = false;
end

if ~isfield(Method.MonteCarlo,'Cluster_Limit')
  Method.MonteCarlo.Cluster_Limit = inf*(1:Method.order); 
elseif length(Method.MonteCarlo.Cluster_Limit) < Method.order
  Cluster_Limit = inf*(1:Method.order);
  Cluster_Limit(1:length(Method.MonteCarlo.Cluster_Limit)) = Method.MonteCarlo.Cluster_Limit;
  for ii = length(Method.MonteCarlo.Cluster_Limit):Method.order
    Cluster_Limit(ii) = Method.MonteCarlo.Cluster_Limit(length(Method.MonteCarlo.Cluster_Limit));
  end
  Method.MonteCarlo.Cluster_Limit = Cluster_Limit;
end


if ~isfield(Method.MonteCarlo,'Increment')
  Method.MonteCarlo.Increment = 1000*(1:Method.order);
end
if length(Method.MonteCarlo.Increment)==1
  Method.MonteCarlo.Increment = Method.MonteCarlo.Increment*ones(1,Method.order);
end
if length(Method.MonteCarlo.Increment) < Method.order
  Increment = 1000*(1:Method.order);
  Increment(1:length(Method.MonteCarlo.Increment)) = length(Method.MonteCarlo.Increment);
  for ii = length(Method.MonteCarlo.Increment):Method.order
    Increment(ii) = Method.MonteCarlo.Increment(length(Method.MonteCarlo.Increment));
  end
  Method.MonteCarlo.Increment = Increment;
end

if ~isfield(Method.MonteCarlo,'Fraction')
  Method.MonteCarlo.Fraction = ones(1,Method.order);
 end
if length(Method.MonteCarlo.Fraction)==1
  Method.MonteCarlo.Fraction = Method.MonteCarlo.Fraction*ones(1,Method.order);
end
if length(Method.MonteCarlo.Fraction) < Method.order
    Method.MonteCarlo.Fraction = ones(1,Method.order);
  for ii = Method.order_lower_bound:Method.order
    Fraction(ii) = Method.MonteCarlo.Fraction(1+ii-Method.order_lower_bound);
  end
  Method.MonteCarlo.Fraction = Fraction;
end


if ~isfield(Method.MonteCarlo,'Threshold')
  Method.MonteCarlo.Threshold = ones(1,Method.order)*1e-3;
end
if length(Method.MonteCarlo.Threshold)==1
  Method.MonteCarlo.Threshold = Method.MonteCarlo.Threshold*ones(1,Method.order);
end
if length(Method.MonteCarlo.Threshold) < Method.order
  Threshold = 1000*(1:Method.order);
  for ii = Method.order_lower_bound:Method.order
    Threshold(ii) = Method.MonteCarlo.Threshold(1+ii-Method.order_lower_bound);
  end
  Method.MonteCarlo.Threshold= Threshold;
end

if ~isfield(Method,'seed')
  Method.seed = 42;
end

if ~isfield(Method,'divisions')
  Method.divisions = 'numSpins';
end

if ~isfield(Method,'startSpin')
  Method.startSpin = 0;
end
if ~isfield(Method,'endSpin')
  Method.endSpin = inf;
end

if ~isfield(Method,'record_clusters')
  Method.record_clusters = false;
end

% Toggle between calculating the spin Hamiltonian or pairwise couplings.
if ~isfield(Method,'HamiltonianType')
  Method.HamiltonianType = 'pairwise';
end

if ~isfield(Method,'precalculateHamiltonian')
  Method.precalculateHamiltonian = false;
  % If there is not enough memory to save the entire Hamiltonian at once
  % setting precalculateHamiltonian to false may allow the calculation to
  % proceed.
end

% allowing for alternate inputs
if strcmp(Method.method,'restrictedCE'),  Method.method = 'rCE';  end
if strcmp(Method.method,'restrictedCCE'),  Method.method = 'rCCE';  end

% The methods rCE and rCE do not use precomputed Hamiltonians.
if strcmp(Method.method,'rCE')||strcmp(Method.method,'rCCE')
  Method.precalculateHamiltonian = false;
end

if ~isfield(Method,'getNuclearContributions')
  Method.getNuclearContributions = false;
end

% Base Units
if ~isfield(System,'joule')
  System.joule = 1; % J;
end
if ~isfield(System,'meter')
  System.meter = 1; % 1; % m.
end
if ~isfield(System,'second')
  System.second = 1; % 1; % s.
end
if ~isfield(System,'tesla')
  System.tesla = 1; % T.
end
if ~isfield(System,'kelvin')
  System.kelvin = 1; % K.
end

System.coulomb = System.joule*System.second/System.tesla/System.meter^2;
System.volt = System.joule/System.coulomb;

% Physical Constants (https://physics.nist.gov/cuu/Constants/index.html)
% constant = SI value * SI units % SI units
System.c = 299792458.0*System.meter/System.second; % m/s.
System.hbar = 1.054571800e-34*System.joule*System.second; % 1.054571800e-34; % J s.
System.h = 6.626070040e-34*System.joule*System.second; % 6.626070040e-34; % J s.
System.muN = 5.050783699e-27*System.joule/System.tesla; % J/T.
System.muB = 927.400e-26*System.joule/System.tesla; % 927.400e-26; % J/T.
System.mu0 = (4*pi*1e-7)*System.meter^3*System.tesla^2/System.joule; % 1.2566e-06; % J^-1 m^3 T^2. % 1.2566e-06
System.kB = 1.38064852e-23*System.joule/System.kelvin; % 1.38064852e-23; % J/K.
System.e = 1.6021766208e-19*System.coulomb;
System.eV = System.e*System.volt; % joule
System.barn = 1e-28*System.meter^2;
% Other Constants
System.angstrom = System.meter*1e-10; % m.
System.wavenumber = System.h*(100*System.c); % J*cm;
System.avogadro= 6.022140857e23;

% System Constants

% Set magnetic field.
if ~isfield(System,'magneticField')
  System.magneticField = 1.2*System.tesla;
end

if ~isfield(System,'inner_radius')
  System.inner_radius = 0;
end
% Set temperature.
if ~isfield(System,'temperature')
  System.temperature = 20 ; % K
end

System.kT = System.temperature*System.kB;


if ~isfield(System,'Methyl')
  System.Methyl = struct;
end
if ~isfield(System.Methyl,'include')
  System.Methyl.include = true;
end
if ~isfield(System.Methyl,'moment_of_inertia')
  System.Methyl.moment_of_inertia =  (5.3373e-47)*System.joule*System.second^2; % kg m^2.;
end
if ~isfield(System,'methyl_V3')
  System.Methyl.V3 = 86*1e-3*System.eV;
end
if ~isfield(System,'tunnel_spliting')
  System.Methyl.omega_harmonic_oscillator = sqrt(9*System.Methyl.V3/System.Methyl.moment_of_inertia);
  
  System.Methyl.instanton_action = 8*System.Methyl.moment_of_inertia*System.Methyl.omega_harmonic_oscillator/9;
  
  System.Methyl.K = 4/3*System.Methyl.omega_harmonic_oscillator^(3/2)*sqrt(System.Methyl.moment_of_inertia/pi/System.hbar);
  
  System.Methyl.tunnel_splitting = 3*System.Methyl.K*exp(-System.Methyl.instanton_action/System.hbar);
  
end
if ~isfield(System.Methyl,'temperature')
  System.Methyl.temperature = System.temperature;
end
System.Methyl.kT = System.Methyl.temperature*System.kB;

% Set electronic spin.
if ~isfield(System.Electron,'spin')
  System.Electron.spin = 1/2;
end

% Set g matrix.
if ~isfield(System.Electron,'g')
  System.Electron.g = 2.0023;
end
if ~isfield(System,'g')
  System.g = 2.0023*[1,1,1];
end
System.gMatrix_gFrame = diag(System.g);

System.omega_Larmor = System.Electron.spin*System.muB*max(max(abs(System.gMatrix_gFrame)))*System.magneticField/System.hbar;

System.Electron.partition_function = 0;
for ii = 0:(2*System.Electron.spin)
System.Electron.partition_function = System.Electron.partition_function + exp((-System.Electron.spin+ii)*System.omega_Larmor*System.hbar/System.kT);
end

System.Electron.State = zeros(1,2*System.Electron.spin+1);
for ii = 1:(2*System.Electron.spin+1)
System.Electron.State(ii) = exp((-System.Electron.spin+ii-1)*System.omega_Larmor*System.hbar/System.kT)/System.Electron.partition_function;
end
% Set pulse sequence.
if ~isfield(System,'experiment')
  System.experiment = 'Hahn';
end



% Set up time grid
if isfield(System,'timepoints') && isfield(System,'dt')
  System.Time = 0:System.dt:(System.timepoints - 1)*System.dt;
elseif isfield(System,'Time')
  System.timepoints = length(System.Time);
  System.dt = abs(System.Time(2) - System.Time(1));
  System.Time = 0:System.dt:(System.timepoints - 1)*System.dt;
  warning('System.Time is not recommended. System.timepoints and System.dt are recommeded instead.');
else
  error('System.timepoints and System.dt are required.');
end

if ~isfield(System,'t0')
  System.t0 = 0;
end
System.Time = System.Time + System.t0;

switch System.experiment
  case 'FID'
    System.dimensionality = 1;
  case 'Hahn'
    System.dimensionality = 1;
  case 'CPMG'
    System.dimensionality = 1;
    if ~isfield(System,'dt_')
      System.dt_ = System.dt;
    end
    System.Time_ = 0:System.dt_:(System.timepoints - 1)*System.dt_;
  case 'CPMG-const'
    System.dimensionality = 1;
  case 'CPMG-2D'
    System.dimensionality = 2;
  otherwise
  error('The experiment ''%s'' is not supported.',System.experiment);
end  
if ~isfield(System,'averaging')
  System.averaging = 'powder';
end

if ~isfield(System,'gridSize') || isempty(System.gridSize)
  System.gridSize = 14;
end


if ~isfield(System,'electron_Zeeman')
  System.electron_Zeeman = true;
end

if ~isfield(System,'nuclear_Zeeman')
  System.nuclear_Zeeman = true;
end

if ~isfield(System,'hyperfine')
  System.hyperfine = [true false];
end

if ~isfield(System,'nuclear_dipole')
  System.nuclear_dipole = [true true true true];
end

if ~isfield(System,'nuclear_quadrupole')
  System.nuclear_quadrupole = true;
end
if ~isfield(System,'nuclear_quadrupole_scale_e2qQh')
  System.nuclear_quadrupole_scale_e2qQh = 1;
end
if ~isfield(System,'nuclear_quadrupole_scale_eta')
  System.nuclear_quadrupole_scale_eta = 1;
end

if ~isfield(System,'theory')
  System.theory = [System.electron_Zeeman,...
    System.nuclear_Zeeman,...
    System.hyperfine(1), System.hyperfine(2), ...
    System.nuclear_dipole(1), System.nuclear_dipole(2), ...
    System.nuclear_dipole(3), System.nuclear_dipole(4), ...
    System.nuclear_quadrupole];
end

% System limiting options
if ~isfield(System,'limitToSpinHalf')
  System.limitToSpinHalf = false;
end

if ~isfield(System,'solventOnly')
  System.solventOnly = false;
end

if ~isfield(System,'D2O')
  System.D2O = false;
end

if ~isfield(System,'deuterateProtein')
  System.deuterateProtein = false;
end

if isfield(System,'deuterateAll') && islogical(System.deuterateAll) && System.deuterateAll
  System.deuterateProtein = true;
  System.D2O = true;
end


% The methods rCE and rCE do not use precomputed Hamiltonians.
if strcmp(Method.method,'rCE')||strcmp(Method.method,'rCCE')
  System.limitToSpinHalf = true;
end

if System.limitToSpinHalf
  disp('Based on the input options, the simulation will only include spin-1/2 nuclei.');
end


if isfield(Method,'mixed_eState') && Method.mixed_eState
  
  if ~isfield(System,'Detection_Operator')
    System.Detection_Operator = spinRaise(System.Electron.spin)/System.Electron.spin^2;
  end
  
  if ~isfield(System,'Pulse')
    System.Pulse = cos(pi/4)*eye(2*System.Electron.spin + 1) + 1i*sin(pi/4)*spinX(System.Electron.spin)/System.Electron.spin;
    System.Pulse(:,:,2) =cos(pi/2)*eye(2*System.Electron.spin + 1) + 1i*sin(pi/2)*spinX(System.Electron.spin)/System.Electron.spin;
  end
  
  if isfield(System,'Flip_Angles')
    ii=1;
    System.Pulse = cos(System.Flip_Angles(ii)/2)*eye(2*System.Electron.spin + 1) + 1i*sin(System.Flip_Angles(ii)/2)*spinX(System.Electron.spin)/System.Electron.spin;
    for ii = 1:length(System.Flip_Angles)
      System.Pulse(:,:,ii) = cos(System.Flip_Angles(ii)/2)*eye(2*System.Electron.spin + 1) + 1i*sin(System.Flip_Angles(ii)/2)*spinX(System.Electron.spin)/System.Electron.spin;
    end
  end
  
  if ~isfield(System,'full_Hyperfine_Tensor')
    System.full_Hyperfine_Tensor = false;
  end
  
else
  Method.mixed_eState = false;
end

if ~isfield(Method,'vectorized')
  Method.vectorized = false;
end

% save options
if ~isfield(Data,'saveLevel')
  Data.saveLevel = 0;
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
