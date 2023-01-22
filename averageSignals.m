%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [SignalMean, experiment_time, ...
  TM_powder,Order_n_SignalMean,Nuclei,statistics] ...
  = averageSignals( ...
  System,Method,Data,OutputData,...
  experiment_time, Nuclei,Clusters,...
  Alpha,Beta,gridWeight,GridInfo,...
  uncertainty,Progress)

Method.verbose = Method.verbose;
if Method.verbose
  fprintf('  Computing cluster Hamiltonians and cluster signals.\n');
end

nOrientations = numel(Alpha);
if Method.verbose
  fprintf("Calculating signal for %d orientations.\n",nOrientations);
end

nTimePoints = sum(System.nPoints);
[SignalsToCalculate,SignalMean,Order_n_SignalMean,...
  TempSignals,Temp_Order_n_Signals,graphs,...
  Signals,Order_n_Signals,Statistics,TM,AuxiliarySignal,...
  Calculate_Signal,Ori_Clusters,Progress,saveAll,numberOfSignals] ...
  = initializeSignals(System, Method, Data, OutputData,...
  nOrientations,Progress,nTimePoints);

pool = setUpPool(Method);

% Loop over orientations
if Method.parallelComputing
  
  parfor iOri = 1:numberOfSignals
    
    [TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics{iOri},...
      graphs{iOri},Ori_Clusters{iOri}] ...
      = getOrientationSignals(System,Method,Nuclei,Clusters, Alpha,Beta, ...
      iOri,Method.verbose,OutputData,Data,Data.InputData,SignalsToCalculate,...
      gridWeight,nOrientations);
    
    TempSignals{iOri} = TempSignals_;
    Temp_Order_n_Signals{iOri} = Temp_Order_n_Signals_;
    
    if saveAll || Method.getContributions
      AuxiliarySignal{iOri} = AuxiliarySignal_;
    end
  end
  
else
  
  for iOri = 1:numberOfSignals
    
    [TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,Statistics{iOri},...
      graphs{iOri},Ori_Clusters{iOri}] ...
      = getOrientationSignals(System,Method,Nuclei,Clusters, Alpha,Beta, ...
      iOri,Method.verbose,OutputData,Data,Data.InputData,SignalsToCalculate,...
      gridWeight,nOrientations);
    
    TempSignals{iOri} = TempSignals_;
    Temp_Order_n_Signals{iOri} = Temp_Order_n_Signals_;
    

      
    TempSignals{iOri} = TempSignals_;
    Temp_Order_n_Signals{iOri} = Temp_Order_n_Signals_;

    if saveAll || Method.getContributions
      AuxiliarySignal{iOri} = AuxiliarySignal_;
    end
      

    
  end
  
end


% ========================================================================
% Save reselts
% ========================================================================
[SignalMean, experiment_time, ...
  TM_powder,Order_n_SignalMean,Nuclei,statistics] ...
= saveSignalResults(Nuclei, Clusters,System, Method, Data,...
SignalMean,Order_n_SignalMean,...
  experiment_time,OutputData,uncertainty,GridInfo,...
  TempSignals,Temp_Order_n_Signals,graphs,...
  Signals,Order_n_Signals,Statistics,TM,AuxiliarySignal,...
  Calculate_Signal,Ori_Clusters,Progress,nOrientations);

if ~isempty(pool)
  delete(pool);
end

% Delete temporary files
for iOri = 1:nOrientations
  temp_file = ['temp_', OutputData, '_sig_', num2str(iOri), '.mat'] ;
  if isfile(temp_file)
    delete(temp_file);
  end
end


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [pool,pc] = setUpPool(Method)

if ~Method.parallelComputing
  pool = [];
  pc = [];
  return;
end


% remove current pool if it exists.
delete(gcp('nocreate'));

% determine number of cores available
numCores = min(feature('numcores'), Method.numberCores);

% create parallel pool
pc = parcluster('local');

if isfield(Method,'JobStorageLocation')
  pc.JobStorageLocation = Method.JobStorageLocation;
end

pool = parpool(numCores);


end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [SignalsToCalculate,SignalMean,Order_n_SignalMean,...
  TempSignals,Temp_Order_n_Signals,graphs,...
  Signals,Order_n_Signals,Statistics,TM,AuxiliarySignal,...
  Calculate_Signal,Ori_Clusters,Progress,saveAll,numberOfSignals] = ...
  initializeSignals(System, Method, Data, OutputData,...
  nOrientations,Progress,nTimePoints)

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
  SignalMean = zeros(1,nTimePoints^System.dimensionality,...
    length(System.deuteriumFraction));
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
  if Method.partialSave && isfile(temp_file) ...
      && ~(saveAll || Method.getContributions)
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
      if Method.verbose, fprintf(['Loading signal %d from ', ...
          temp_file, '.\n'],iOri); end
      
      % set values to simulation variables
      Calculate_Signal{iOri} = false;
      

      Signals{iOri} = signal;
      Order_n_Signals{iOri} = order_n;
      Statistics{iOri} = statistics;
      if Method.Ori_cutoffs
        Ori_Clusters{iOri} = iOri_Clusters;
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
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% ========================================================================
% Calculates signal for a set of orientations
% ========================================================================
function [TempSignals_, AuxiliarySignal_,Temp_Order_n_Signals_,...
  Statistics_isignal,graphs_isignal,iOri_Clusters] ...
    = getOrientationSignals(System,Method,Nuclei,inClusters, Alpha,Beta,...
    isignal,verbose,OutputData,Data,InputData,SignalsToCalculate,...
    gridWeight,iSignal_max)
     
  
  % grid point indices
  if System.newIsotopologuePerOrientation && ~Method.reparseNuclei

    if Method.useCentralSpinSystem
      error(['Error in getOrientationSignals(): ',...
        'Method.useCentralSpinSystem is not yet compatible with',...
        ' System.newIsotopologuePerOrientation.']);
    end
    Nuclei = newHydronIsotopologue(Nuclei,System);
    if verbose
      fprintf('Generated a new hydron isotopologue for orientation %d.\n',...
        isignal);
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
    Temp_Order_n_Signals_{iorder} ...
      = gridWeight(igrid)*Temp_Order_n_Signals_{iorder}/normalizing_factor;
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
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
    isotopologueStatistics.number_1H_exchangeable ...
      = Nuclei.number_1H_exchangeable;
    isotopologueStatistics.number_1H_nonExchangeable ...
      = Nuclei.number_1H_nonExchangeable;
    isotopologueStatistics.number_2H_exchangeable ...
      = Nuclei.number_2H_exchangeable;
    isotopologueStatistics.number_2H_nonExchangeable ...
      = Nuclei.number_2H_nonExchangeable;
  end
  
  StatisticsTemporary = getPairwiseStatistics(System, Method, Nuclei);
  Nuclei.Statistics = StatisticsTemporary;
  Adjacency = getAdjacencyMatrix(System, Nuclei,Method);
  Nuclei.Adjacency = Adjacency;
  
  Statistics.isotopologueStatistics = isotopologueStatistics;
  Nuclei.Statistics = Statistics;
  
  if ~Method.neighborCutoff.sizeDependent
    Ori_Clusters = findClusters_treeSearch(Nuclei,Method.order,...
      Method.extraOrder,{}, Method);
  end

  for clusterSize = 1:Method.order
    if Method.neighborCutoff.sizeDependent
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
  
  
    fprintf('Found %d orientation clusters of size %d.\n', ...
      size(Clusters{clusterSize},1),clusterSize);
  end
end

% Rotate bath spin tensors.
for inucleus = 1:Nuclei.number
  
  Nuclei.Atensor(inucleus,:) = reshape( ...
    R_pdb2lab*reshape( full(Nuclei.Atensor(inucleus,:))',3,3)*R_pdb2lab' , ...
    9,1);
  if ~System.limitToSpinHalf
    Nuclei.Qtensor(:,:,inucleus) ...
      = R_pdb2lab*Nuclei.Qtensor(:,:,inucleus)*R_pdb2lab';
    % Elementwise Qtensor manipulation used for testing.  
    % The default filter is ones(3);
    Nuclei.Qtensor(:,:,inucleus) ...
      = Nuclei.Qtensor(:,:,inucleus).*System.nuclear_quadrupole_filter;
  end
end

% Rotate the g-matrix.
if isfield(System,'gFrame')
  % Get rotation matrix.
  g2MolRotation ...
    = rotateZYZ(-System.gFrame(1),-System.gFrame(2),-System.gFrame(3));
  
  % Rotate to molecular frame.
  System.gMatrix = g2MolRotation*System.gMatrix_gFrame*g2MolRotation';
end

% Rotate the g-matrix to the lab frame.
System.gMatrix = R_pdb2lab*System.gMatrix_gFrame*R_pdb2lab';

% Get the g-value along the magnetic field direction.
System.Electron.g = System.gMatrix(3,3);

graphs = Nuclei.Adjacency; 

if Method.useMultipleBathStates
  [Nuclei.MeanFieldCoefficients, Nuclei.MeanFieldTotal] ...
    = getMeanFieldCoefficients(Nuclei,System);
else
  Nuclei.MeanFieldCoefficients = [];
  Nuclei.MeanFieldTotal = [];
end


if isfield(Method,'exportHamiltonian') && Method.exportHamiltonian
  
  % Generate Cartesian spin-spin coupling Hamiltonian.
  
  geff=System.gMatrix(3,3);
  
  [Tensors,zeroIndex] = pairwisetensors(Nuclei.Nuclear_g, Nuclei.Coordinates,...
    1:Nuclei.number,Nuclei.Atensor, System.magneticField, System.ge, ...
    geff, System.muB, System.muN, System.mu0, System.hbar,System.theory,[]);
  
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

    [Signal, AuxiliarySignal, Signals] ...
      = calculate_signal(System, Method, Nuclei,Clusters);

    Order_n_Signal = cell(1,Method.order);
    
    for ii = 1:Method.order
      Order_n_Signal{ii} = Signals(1,ii);
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
          = calculateSignal_forHDCCE(System, Method, Nuclei,Clusters);
        
        % Change all hydrons to deuteron.
        System.deuteriumFraction = 1;
        Nuclei = newHydronIsotopologue(Nuclei,System);

        % Get deuteron auxiliary signals.

        [~, Coherences_1D,Coherences_2D,SubclusterIndices_2D,~,~] ...
          = calculateSignal_forHDCCE(System, Method, Nuclei,Clusters);
        
        
        [~, ...
          ~,~] = ...
          doHydrogenIsotopologueCCE(...
          Coherences_1H,Coherences_2H,Coherences_1D,Coherences_2D, ...
          fractions, ClusterArray, ...
          SubclusterIndices_2H,SubclusterIndices_2D,...
          timepoints,dimensionality, Method.order,...
          Nuclei.numberClusters,Nuclei.Exchangeable,Nuclei.MoleculeID);
        
        System.deuteriumFraction = fractions;

      else

        [~, AuxiliarySignal_ofOrder, ~] ...
          = calculate_signal(System, Method, Nuclei,Clusters);
     
      end
      % Record the appropraite signals.
      for record_order = iorder:new_order
        
        AuxiliarySignal{record_order} = AuxiliarySignal_ofOrder{record_order};
       
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
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>





%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [SignalMean, experiment_time, ...
  TM_powder,Order_n_SignalMean,Nuclei,statistics] ...
= saveSignalResults(Nuclei, Clusters,System, Method, Data,...
SignalMean,Order_n_SignalMean,...
  experiment_time,OutputData,uncertainty,GridInfo,...
  TempSignals,Temp_Order_n_Signals,graphs,...
  Signals,Order_n_Signals,Statistics,TM,AuxiliarySignal,...
  Calculate_Signal,Ori_Clusters,Progress,nOrientations)

tempSignalIndex = 0;

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




if System.doPruneNuclei && ~Method.reparseNuclei
  statistics = Nuclei.Isotopologue;
else
  statistics.Statistics = Statistics;
  statistics.uncertainty = uncertainty;
end
Nuclei.graphs = graphs;


if System.newIsotopologuePerOrientation
  Nuclei.PowderStatistics.Isotopologue.Mean_TypeNumber = [0,0,0];
  Nuclei.PowderStatistics.Isotopologue.Mean_Instance_2H_Number = [0,0,0];
  Nuclei.PowderStatistics.Isotopologue.Mean_Instance_2H_Fraction = [0,0,0];
end

for iOri = 1:nOrientations
  if ~Method.reparseNuclei && System.newIsotopologuePerOrientation
    Nuclei.PowderStatistics.Isotopologue.Mean_TypeNumber ...
      = Nuclei.PowderStatistics.Isotopologue.Mean_TypeNumber ...
      + Statistics{iOri}.Isotopologue.TypeNumber/nOrientations;

    Nuclei.PowderStatistics.Isotopologue.Mean_Instance_2H_Number ...
      = Nuclei.PowderStatistics.Isotopologue.Mean_Instance_2H_Number ...
      + Statistics{iOri}.Isotopologue.Instance_2H_Number/nOrientations;

    Nuclei.PowderStatistics.Isotopologue.Mean_Instance_2H_Fraction ...
      = Nuclei.PowderStatistics.Isotopologue.Mean_Instance_2H_Fraction ...
      + Statistics{iOri}.Isotopologue.InstanceFraction/nOrientations;
  end
  if Method.Ori_cutoffs
    Clusters =  combineClusters(Clusters,Ori_Clusters{iOri});
  end
  SignalMean = SignalMean + Signals{iOri};
  %Signals{isignal} = abs(Signals{isignal});
  for iorder = 1:Method.order
    if min( numel(Order_n_Signals),numel(Order_n_SignalMean)) < iorder...
        || ischar(Order_n_SignalMean{iorder}) ...
        || isempty(Order_n_Signals{iorder})
      continue;
    end

    try
      Order_n_SignalMean{iorder} ...
        = Order_n_SignalMean{iorder} + Order_n_Signals{iOri}{iorder};
    catch
      fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');

      fprintf('Error in Order %d mean at signal number %d.\n',iorder,iOri);
      disp('Could not evaluate');
      disp(['Order_n_SignalMean{iorder} = ',...
        'Order_n_SignalMean{iorder} + Order_n_Signals{isignal}{iorder};']);
      fprintf(['Order_n_SignalMean{%d} = ',...
        'Order_n_SignalMean{%d} + Order_n_Signals{%d}{%d};\n'], ...
        iorder,iorder,iOri,iorder);
      if exist('Order_n_Signals','var')
        if iscell(Order_n_Signals) && length(Order_n_Signals)>=iOri

          if iscell(Order_n_Signals{iOri}) ...
              && length(Order_n_Signals{iOri}) >= iorder
            fprintf('Length of Order_n_Signals{%d}{%d} is %d.\n',...
              iOri,iorder, length(Order_n_Signals{iOri}{iorder}));

            fprintf('Length of Order_n_SignalMean{%d} is %d.\n',...
              iorder, length(Order_n_SignalMean{iorder})  );

          else
            fprintf('Order_n_Signals{%d}{%d} does not exist.\n',...
              iOri,iorder);
          end

        else
          fprintf('Order_n_Signals{%d} does not exist.\n',iOri);
        end
      else
        fprintf('Order_n_Signals does not exist.\n');
      end

      disp('Attempting to recover data...')



      if iorder == 1

        Order_n_Signals{iOri}{iorder} ...
          = gridWeight(iOri)*ones(size(experiment_time));

        disp('Recovered.');

      elseif iorder == Method.order
        Order_n_Signals{iOri}{iorder} = Signals{iOri};
        disp('Recovered.');
      else
        Progress.Order_n_Mean = false;
        disp('Failed.');
        Order_n_Signals{iOri}{iorder} = nan(size(experiment_time));
      end

      Order_n_SignalMean{iorder} ...
        = Order_n_SignalMean{iorder} + Order_n_Signals{iOri}{iorder};
      fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    end

  end
  % Signals{isignal} = Signals{isignal}/max(Signals{isignal});
end
for iSize=1:Method.order
  Nuclei.numberClusters(iSize) = size(Clusters{iSize},1);
end


% update progress
if Progress.Order_n_Mean, Progress.Order_n_Mean = true; end

if System.dimensionality==2 && min(size(SignalMean))==1
  SignalMean = reshape(SignalMean',nTimePoints,nTimePoints)';
  for iorder = 1:Method.order
    Order_n_SignalMean{iorder} ...
      = reshape(Order_n_SignalMean{iorder}',nTimePoints,nTimePoints)';
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

if Method.verbose
  fprintf('Processing\n');
end
uncertainty = [];
% decide what to save
if ~isempty(OutputData)
 
  T = array2table([experiment_time',SignalMean']);
  T.Properties.VariableNames(1:2) = {'time','signal'};
  writetable(T,[OutputData(1:end-4) , '.csv']);
  if Data.save_mat_file
    switch Data.saveLevel
      case 0
        save(OutputData,'SignalMean','Signals','TM','TM_powder',...
          'Progress','uncertainty','GridInfo','-append');

      case 1
        save(OutputData,'SignalMean','Signals','TM','TM_powder','Progress',...
          'Nuclei','Order_n_SignalMean','Order_n_Signals','uncertainty',...
          'GridInfo','AuxiliarySignal','-append');

      case 2
        save(OutputData,'-v7.3');
    end
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

if Method.verbose
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
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [Signal,AuxiliarySignal,Order_n_Signal] = ...
  doLCE(System,Method, Nuclei, Clusters, verbose)

MethylID = [];

Signal = zeros(size(System.Time));
V2 = Signal;
Order_n_Signal{1} = Signal;
tau1 = System.Time;
tau2 = System.Time';

if Method.conserveMemory
  AuxiliarySignal = 'The auxiliary signals are not saved in memory conservation mode.';
end

for isize = 2
  for icluster = 1:Nuclei.numberClusters(isize)
    
    thisCluster = Clusters{isize}(icluster,:);
    
    [tensors,zeroIndex] = pairwisetensors(Nuclei.Nuclear_g, Nuclei.Coordinates,...
      thisCluster,Nuclei.Atensor, System.magneticField, System.ge, System.gMatrix(3,3), System.muB, System.muN, System.mu0, System.hbar,System.theory,MethylID);
    
    b = -2*pi*tensors(3,3,2,3)/4; % rad/s.
    omega = 2*pi*(tensors(3,3,1,2) -  tensors(3,3,1,3))/2; % rad/s.
    
    
    switch System.experiment
      case 'Hahn'
      case 'CPMG-2D'
        V2 = - 4*(b/omega*(cos(omega.*tau1) - cos(omega.*tau2)  )).^2;
        AuxiliarySignal_  = V2;
        Signal = Signal + AuxiliarySignal_;
    end
    % initialize the nth auxiliary signal
    if ~Method.conserveMemory
      AuxiliarySignal{icluster} = AuxiliarySignal_;
    end
  end
  

  
  Signal = exp(Signal);
  if System.dimensionality ==2
    Signal = reshape(Signal.',1,[]);
  end
  Order_n_Signal{1} = Signal;
  Order_n_Signal{isize} = Signal;
  
  
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Calculate < psi | H | psi > for a given | psi >,
% where  | psi > = Kron_{n} |n>.
function [h,hTr] = getMeanFieldCoefficients(Nuclei,System)

% Only operators proportional Iz are needed.
theory = System.theory;
% useEZ       = theory(1);
useNZ       = false; %  theory(2);
useHF_SzIz  = false; %  theory(3);
% useHF_SzIxy = theory(4);
useNucA     = theory(5);
% useNucB     = theory(6);
useNucCD    = true; % theory(7);
% useNucEF    = theory(8);
useNQ       = false; %  theory(9);

nStates = max(Nuclei.nStates);
h_ = zeros(Nuclei.number,Nuclei.number,4);
h = zeros(Nuclei.number,Nuclei.number,4,nStates);
hTr = zeros(1,nStates);
% ENUM
E = 1; Z = 2; RAISE = 3; SZ = 4;


% Constants
km = System.mu0/(4*pi);
muB = System.muB;
muN = System.muN;
hbar = System.hbar;
ge = System.Electron.g;

% Loop through bath states.
for istate = 1:nStates
  
  % Loop through spin states.
  for iSpin = 1:Nuclei.number
    
    % Get the ith soin state.
    psi_i = Nuclei.ZeemanStates(istate, iSpin);
    
    % Get the ith spin.
    I = Nuclei.Spin(iSpin);
    
    % Determine the z-projection.
    MI = psi_i - I - 1;
    
    % Nuclear Zeeman
    if useNZ
      
      NuclearZeeman = MI*Nuclei.Nuclear_g(iSpin)*System.magneticField*System.muN; % J.
      NuclearZeeman = NuclearZeeman/(2*pi*hbar);
      h_(iSpin,iSpin,E) = h_(iSpin,iSpin,E) + NuclearZeeman;
      
    end
    
    % Hyperfine
    if useHF_SzIz
      
      
      if strcmp(Nuclei.Type{iSpin},'e')
        muN = -muB;
      end
      
      gni = Nuclei.Nuclear_g(iSpin);
      r = Nuclei.Coordinates(iSpin,:);
      if size(r,2)==3
        r=r';
      end
      
      n = r/norm(r);
      nnt = n*n';
      r3 = norm(r)^3;
      dd = km*ge*muB*gni*muN/r3*(eye(3)-3*nnt);
      Hyperfine = -dd/(2*pi*hbar); % Hz.
      
      h_(iSpin,iSpin,SZ) = h_(iSpin,iSpin,SZ) + Hyperfine(3,3)*MI;
      
    end
    
    % Quadrupole
    
    if useNQ && Nuclei.StateMultiplicity(iSpin) > 2
      Q_ = Nuclei.Qtensor(:,:,iSpin);
      H_nuclear_quadrupole = MI*Q_(3,3)*MI ...
        + (Q_(1,1) + Q_(2,2))*(I*(I + 1)  -1/2*MI^2);
      
      h_(iSpin,iSpin,E) = h_(iSpin,iSpin,E) ...
        + H_nuclear_quadrupole;
      
    end
    
    % Loop through other bath spins.
    for jSpin = 1:iSpin-1
      
      % Get jth spin state.
      psi_j = Nuclei.ZeemanStates(jSpin);
      
      % Get the jth spin.
      J = Nuclei.Spin(jSpin);
      
      % Determine the z-projection.
      MJ = psi_j - J -1;
      
      if useNucA || useNucCD
        
        % Get classical dipole moment.
        muNi = System.muN;
        
        % Check if the spin is an electron/
        if strcmp(Nuclei.Type{iSpin},'e')
          % Adjust to Bohr magneton.
          muNi = -muB;
        end
        
        % Get classical dipole moment.
        muNj = System.muN;
        
        % Check if the spin is an electron/
        if strcmp(Nuclei.Type{jSpin},'e')
          % Adjust to Bohr magneton.
          muNj = -muB;
        end
        
        % Get spin g-factors.
        gni = Nuclei.Nuclear_g(iSpin);
        gnj = Nuclei.Nuclear_g(jSpin);
        
        % Get inter-spin veparation vector.
        r = Nuclei.Coordinates(iSpin,:)-Nuclei.Coordinates(jSpin,:);
        
        % Set r to be a column vector.
        if size(r,2)==3
          r=r';
        end
        
        % Get direction of r.
        n = r/norm(r);
        nnt = n*n';
        r3 = norm(r)^3;
        dd = -km*gni*gnj*muNi*muNj/r3*(eye(3)-3*nnt);
        dd = dd/(2*pi*hbar); % Hz.
        
        % nucleus-nucleus secular (A term)
        if useNucA
          
          h_(iSpin,jSpin,E) = h_(iSpin,jSpin,E) + MI*dd(3,3)*MJ;
          h_(jSpin,iSpin,E) = h_(iSpin,jSpin,E); % + MJ*dd(3,3)*MI;
          
          h_(iSpin,jSpin,Z) = h_(iSpin,jSpin,Z) + MI*dd(3,3);
          h_(jSpin,iSpin,Z) = h_(iSpin,jSpin,Z);% + MJ*dd(3,3);
        end
        
        % nucleus-nucleus dipolar C and D terms
        if useNucCD
          
          cd = 1/2*(dd(1,3) - 1i*dd(2,3));
          
          h_(iSpin,jSpin,RAISE) = h_(iSpin,jSpin,RAISE) + MI*cd;
          h_(jSpin,iSpin,RAISE) = h_(jSpin,iSpin,RAISE) + MJ*cd;
          
        end
      end
    end
  end
  
  for iSpin = 1:Nuclei.number
    h_(iSpin,iSpin,Z) = sum( h_(:,iSpin,Z));
    h_(iSpin,iSpin,RAISE) = sum( h_(:,iSpin,RAISE));
  end
  hTr(istate) = trace(h_(:,:,E));
  
h(:,:,:,istate) = h_;
end
 
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [Signal,AuxiliarySignal,Order_n_Signal] = doRestrictedCCE(System, Method, Nuclei,verbose)

% initialize
Hyperfine = zeros(Nuclei.number,1);
Signal = ones(size(System.Time));
Order_n_Signal{1} = Signal;
numberNuclei=Nuclei.number;
Adjacency =  Nuclei.Adjacency;
if Method.conserveMemory
  AuxiliarySignal = 'The auxiliary signals are not saved in memory conservation mode.';
end

if verbose
  % threshold and counter for how often to display info
  verboseThreshold = 0.05*numberNuclei;
  verboseCounter = 0;
  disp('Starting restricted CCE.');
end

if ~all(Nuclei.Spin==1/2)
  warning('System has spins with I~=1/2. Restricted CCE is limited to spin-1/2 nuclei and will ignore them.');
end

% loop over all nuclei
for inucleus = 1:numberNuclei
  if (abs(Nuclei.Spin(inucleus) -0.5)>1e-6)
    continue;
  end
  
  % initialize the nth auxiliary signal 
  if ~Method.conserveMemory
    AuxiliarySignal{inucleus} = ones(size(Signal));
  end
  % Calculating hyperfine coupling
  
  gamma_e = -System.Electron.g*System.muB/System.hbar;
  Rn = norm(Nuclei.Coordinates(inucleus,:) );
  gamma_n = Nuclei.Nuclear_g(inucleus)*System.muN/System.hbar;
  cosTheta2 = cosFieldAngle([0,0,0],Nuclei.Coordinates(inucleus,:));
  cosTheta2 = cosTheta2*cosTheta2;
  
  Hyperfine(inucleus) = Nuclei.FermiContact(inucleus)-(System.mu0/4/pi)*gamma_n*gamma_e*System.hbar*(1-3*cosTheta2)*Rn^-3;
  % Calculating bath coupling
  for jnucleus = 1:inucleus-1
    
    % skip over I != 1/2
    if Nuclei.Spin(inucleus)~=1/2, continue; end
    
    if ~Adjacency(inucleus,jnucleus), continue; end
    
    
    % calculate dipolar coupling
    %{
    cosThetaSquared = (cosFieldAngle(Nuclei.Coordinates(inucleus,:),Nuclei.Coordinates(jnucleus,:)))^2;
    b = 0.25*(System.mu0/4/pi)*Nuclei.Nuclear_g(inucleus)*Nuclei.Nuclear_g(jnucleus)*System.muN^2; % J m^3.
    r = norm(Nuclei.Coordinates(inucleus,:) - Nuclei.Coordinates(jnucleus,:));
    r3 = r^3;
    b = -b*(3*cosThetaSquared - 1)/r3; % J.
    b = b/(System.hbar); % 1/s.
    
    c = ( Hyperfine(inucleus)-Hyperfine(jnucleus) )/(4*b);
    w = b*sqrt(1+c^2);
    %}
    b = Nuclei.Statistics.Nuclear_Dipole(inucleus,jnucleus)/4*2*pi;
    modDepth = Nuclei.Statistics.Modulation_Depth(inucleus,jnucleus);
    w = Nuclei.Statistics.Frequency_Pair(inucleus,jnucleus)*2*pi;
   
    AuxiliarySignal_ = 1 - modDepth * sin(w*System.Time).^4;
    Signal = Signal.*AuxiliarySignal_;
    if ~Method.conserveMemory
      AuxiliarySignal{inucleus} =AuxiliarySignal{inucleus}.*AuxiliarySignal_;
    end
  end
  if verbose
    verboseCounter = verboseCounter + 1;
    
    if verboseCounter > verboseThreshold
      verboseCounter = 0;
      fprintf('Nuclei: %d/%d, (%s).\n',inucleus,numberNuclei,datetime);
    end
    
  end
  Order_n_Signal{2} = Signal;
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function cos_theta = cosFieldAngle(ParticleCoordinates_1,ParticleCoordinates_2)
R = ParticleCoordinates_2 - ParticleCoordinates_1;
cos_theta = R*[0;0;1]/norm(R);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [Signal,AuxiliarySignal,Order_n_Signal] = doRestrictedCE(System, Method, Nuclei, r0, verbose)

nNuclei = Nuclei.number;

% initialize 
Hyperfine = zeros(nNuclei,1);
ContributionSum = zeros(size(System.Time));
Order_n_Signal{1} = ones(size(System.Time));


if Method.conserveMemory
  AuxiliarySignal = 'The auxiliary signals are not saved in memory conservation mode.';
else
  AuxiliarySignal{nNuclei} = zeros(size(System.Time));
end

if ~all(Nuclei.Spin==1/2)
  warning('System has spins with I~=1/2. Restricted CE is limited to spin-1/2 nuclei and will ignore them.');
end

fprintf('%d of %d nuclei are spin-1/2.',sum(Nuclei.Spin==1/2),nNuclei);

verboseN = 0.01*nNuclei;

for inucleus = 1:nNuclei
  
  if ~Method.conserveMemory
    ln_AuxiliarySignal_ = zeros(size(System.Time));
  end
  
  if verbose && inucleus > verboseN
    verboseN = verboseN + 0.01*nNuclei;
    fprintf('Calculating: %d/%d.\n', inucleus,nNuclei);
  end
  
  if Nuclei.Spin(inucleus)~=1/2, continue; end
  
  % Calculating hyperfine coupling
  gamma_e = -System.Electron.g*System.muB/System.hbar;
  
  Rn = norm(Nuclei.Coordinates(inucleus,:) );
  
  gamma_n = Nuclei.Nuclear_g(inucleus)*System.muN/System.hbar;
  
  cosTheta2 = cosFieldAngle([0,0,0],Nuclei.Coordinates(inucleus,:))^2;
  
  Hyperfine(inucleus) = Nuclei.Hyperfine(inucleus)- ...
    (System.mu0/4/pi)*gamma_n*gamma_e*System.hbar*(1-3*cosTheta2)*Rn^-3;
  
  % Calculating bath coupling
  for jnucleus = 1:(inucleus-1)
    if Nuclei.Spin(jnucleus)~=1/2, continue; end
    
    % Calculate inter-nuclear distance and skip if larger than threshold
    Rij = norm(Nuclei.Coordinates(inucleus,:)-Nuclei.Coordinates(jnucleus,:));
    if Rij > r0, continue; end

    
    
    
    cosThetaSquared = (cosFieldAngle(Nuclei.Coordinates(inucleus,:),Nuclei.Coordinates(jnucleus,:)))^2;
    b = 0.25*(System.mu0/4/pi)*Nuclei.Nuclear_g(inucleus)*Nuclei.Nuclear_g(jnucleus)*System.muN^2; % J m^3
    r = norm(Nuclei.Coordinates(inucleus,:) - Nuclei.Coordinates(jnucleus,:));
    b = -b*(3*cosThetaSquared - 1)/r^3; % J
    b = b/System.hbar; % J -> rad/s
    
    c = ( Hyperfine(inucleus)-Hyperfine(jnucleus) )/(4*b); % unitless
    
    omega = 2*b*sqrt(1+c^2);
    
    if abs(c)>1e9
      error('c is larger than 1e9!');
    end
    
    % find the contribution from iNuc and jNuc to to ln( Signal )
    ln_AuxiliarySignal = -Nuclei.Abundance(inucleus)*Nuclei.Abundance(jnucleus)*( (c/(1+c^2)) * (cos(omega*System.Time) -1) ).^2;
    
    % add to running sum of ln( Signal ).
    ContributionSum = ContributionSum +ln_AuxiliarySignal;

    if ~Method.conserveMemory
      % add to running sum of ln( AuxiliarySignal )
      ln_AuxiliarySignal_ = ln_AuxiliarySignal_ + ln_AuxiliarySignal;
    end
    
  end
  
  if ~Method.conserveMemory
    % exponentiate to turn AuxiliarySignal into a true auxiliary signal
    AuxiliarySignal{inucleus} = exp(ln_AuxiliarySignal_);
  end
  
end

Signal = exp(ContributionSum);
Order_n_Signal{2} = Signal;

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>