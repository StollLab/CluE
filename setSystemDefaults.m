
% ========================================================================
% New Function
% ========================================================================

function [System,Method, Data, statistics] = setSystemDefaults(System,Method, Data)

if ~isfield(Method,'errorTolerance')
  Method.errorTolerance = 1e-9;
end
if ~isfield(Method,'useCentralSpinSystem')
  Method.useCentralSpinSystem = false;
end

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
  Method.getClusterContributions = false;
  Method.getNuclearStatistics = false;
  % The following behaviors are are used regardless of settings,
  % so the settings are modified to reflect what is done.
  
  Method.precalculateHamiltonian = false;
  
end
if ~isfield(Method,'getNuclearStatistics')
  Method.getNuclearStatistics = false;
end
if ~isfield(Method,'gpu')
  Method.gpu = false;
end

if ~isfield(Method,'exportHamiltonian')
  Method.exportHamiltonian = false;
end
if ~isfield(Method,'exportClusters')
  Method.exportClusters = false;
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
if ~isfield(Method,'slurm')
  Method.slurm = false;
end
if ~isfield(Method,'numberCores')
  Method.numberCores = inf;
end
% Set verbosity.
if ~isfield(Method,'verbose')
  Method.verbose = false;
end

% cutoff criteria
if ~isfield(Method,'order')
  Method.order = 2;
end
if ~isfield(Method,'allowHDcoupling')
  Method.allowHDcoupling = false;
end
if ~isfield(Method,'Criteria') || isempty(Method.Criteria)
  Method.Criteria = {'dipole'};
end

num_criteria = numel(Method.Criteria);
for ii = 1:num_criteria
  switch Method.Criteria{ii}
    case 'methyl only'
      if ~isfield(Method.cutoff,'methylOnly')
        Method.cutoff.methylOnly = true(1,Method.order);
      end
      break;
      
    case 'methyl coupled only'
      if ~isfield(Method.cutoff,'methylCoupledOnly')
        Method.cutoff.methylCoupledOnly = true(1,Method.order);
      end
      break;
      
    case 'distance'
      if ~isfield(Method.cutoff,'rMax')
        Method.cutoff.rMax = inf;
      end
      if ~isfield(Method.cutoff,'rMin')
        Method.cutoff.rMin = 0;
      end
      break;
  end
end



zer = zeros(1,Method.order);
if ~isfield(Method,'cutoff') 
  Method.cutoff.modulation = zer;
  Method.cutoff.dipole = zer;
  Method.cutoff.dipoleHalf = zer;
  Method.cutoff.dipoleOne = zer;
  Method.cutoff.bAmax = zer;
  Method.cutoff.maxAmax = inf;
  Method.cutoff.minAmax = zer;
  Method.cutoff.max_distance = inf + zer;
  Method.cutoff.min_distance = zer;
  Method.cutoff.hyperfine_sup = inf + zer;
  Method.cutoff.hyperfine_inf = zer;
  Method.cutoff.methylOnly = false(1,Method.order);
  Method.cutoff.methylCoupledOnly = false(1,Method.order);
end
if ~isfield(Method.cutoff,'bAmax') 
  Method.cutoff.bAmax = zer;
end
if ~isfield(Method.cutoff,'dipole') 
  Method.cutoff.dipole = zer;
end
if ~isfield(Method.cutoff,'dipoleHalf') 
  Method.cutoff.dipoleHalf = zer;
end
if ~isfield(Method.cutoff,'dipoleOne') 
  Method.cutoff.dipoleOne = zer;
end
if ~isfield(System,'radius_nonSpinHalf')
  System.radius_nonSpinHalf = System.radius;
end
if ~isfield(Method.cutoff,'radius_nonSpinHalf') 
  Method.cutoff.radius_nonSpinHalf = System.radius_nonSpinHalf;
else
  System.radius_nonSpinHalf = Method.cutoff.radius_nonSpinHalf;
end

if ~isfield(Method,'lock_bAmax')
  Method.lock_bAmax = false;
end
if ~isfield(Method,'Ori_cutoffs')
  Method.Ori_cutoffs = false;
end
if ~isfield(Method,'reparseNuclei')
  Method.reparseNuclei = false;
end
if ~isfield(Method.cutoff,'methylOnly')
  Method.cutoff.methylOnly = false(1,Method.order);
end
if numel(Method.cutoff.methylOnly) < Method.order
  n_ = numel(Method.cutoff.methylOnly);
  Method.cutoff.methylOnly(n_:Method.order) = Method.cutoff.methylOnly(n_);
end

if ~isfield(Method.cutoff,'methylCoupledOnly')
  Method.cutoff.methylCoupledOnly = false(1,Method.order);
end
if numel(Method.cutoff.methylCoupledOnly) < Method.order
  n_ = numel(Method.cutoff.methylCoupledOnly);
  Method.cutoff.methylCoupledOnly(n_:Method.order) = Method.cutoff.methylCoupledOnly(n_);
end
if numel(Method.cutoff.dipole) < Method.order
  n_ = numel(Method.cutoff.dipole);
  Method.cutoff.dipole(n_:Method.order) = Method.cutoff.dipole(n_);
end
if numel(Method.cutoff.dipoleHalf) < Method.order
  n_ = numel(Method.cutoff.dipoleHalf);
  Method.cutoff.dipoleHalf(n_:Method.order) = Method.cutoff.dipoleHalf(n_);
end
if numel(Method.cutoff.dipoleOne) < Method.order
  n_ = numel(Method.cutoff.dipoleOne);
  Method.cutoff.dipoleOne(n_:Method.order) = Method.cutoff.dipoleOne(n_);
end
if isfield(Method.cutoff, 'maxAmax') && numel(Method.cutoff.maxAmax) < Method.order
  n_ = numel(Method.cutoff.maxAmax);
  Method.cutoff.maxAmax(n_:Method.order) = Method.cutoff.maxAmax(n_);
end
if isfield(Method.cutoff, 'minAmax') && numel(Method.cutoff.minAmax) < Method.order
  n_ = numel(Method.cutoff.minAmax);
  Method.cutoff.minAmax(n_:Method.order) = Method.cutoff.minAmax(n_);
end
if isfield(Method.cutoff, 'minimum_frequency') ...
 && numel(Method.cutoff.minimum_frequency) < Method.order
  n_ = numel(Method.cutoff.minimum_frequency);
  Method.cutoff.minimum_frequency(n_:Method.order) ...
  = Method.cutoff.minimum_frequency(n_);
end
if isfield(Method.cutoff, 'modulation') ...
 && numel(Method.cutoff.modulation) < Method.order
  n_ = numel(Method.cutoff.modulation);
  Method.cutoff.modulation(n_:Method.order) = Method.cutoff.modulation(n_);
end
if isfield(Method.cutoff, 'DeltaHyperfine') ...
 && numel(Method.cutoff.DeltaHyperfine) < Method.order
  n_ = numel(Method.cutoff.DeltaHyperfine);
  Method.cutoff.DeltaHyperfine(n_:Method.order) ...
   = Method.cutoff.DeltaHyperfine(n_);
end

if numel(Method.cutoff.bAmax) < Method.order
  n_ = numel(Method.cutoff.bAmax);
  Method.cutoff.bAmax(n_:Method.order) = Method.cutoff.bAmax(n_);
end

if norm(Method.cutoff.dipole-Method.cutoff.dipole(1)) > 0
  Method.cutoff.sizeDependent = true;
else
  Method.cutoff.sizeDependent = false;
end

% Radius to load nuclei in to.
if ~isfield(System,'load_radius')
  System.load_radius = System.radius;
end
System.load_radius = max(System.load_radius, System.radius);

% cluster order max
if ~isfield(Method,'order')
  Method.order = 2;
end
if ~isfield(Method,'clusterization')
  Method.clusterization = 'tree-search';
end
if ~isfield(Method,'combineClusters')
  Method.combineClusters = false;
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
    Fraction = ones(1,Method.order);
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
if ~isfield(Method,'useMultipleBathStates')
  Method.useMultipleBathStates = false;
end
if ~isfield(Method,'useInterlacedClusters')
  Method.useInterlacedClusters = false;
end

if ~isfield(Method,'emptyClusterSetsOkay')
  Method.emptyClusterSetsOkay = false;
end
% if any(Method.useInterlacedClusters(:)) && Method.useMultipleBathStates
%   error('Method.useInterlacedClusters and Method.useMultipleBathStates are not compatiple.')
% end
if (size(Method.useInterlacedClusters,1) ~= size(Method.useInterlacedClusters,2)) || (size(Method.useInterlacedClusters,1) < Method.order)

  if any(Method.useInterlacedClusters(:))
  M_ = eye(Method.order)>0;
  nrow_ = min(Method.order, size(Method.useInterlacedClusters,1));
  ncol_ = min(Method.order, size(Method.useInterlacedClusters,1));
  M_(1:nrow_,1:ncol_) = Method.useInterlacedClusters;
  Method.useInterlacedClusters = M_;
  else
    Method.useInterlacedClusters = false(Method.order);
  end
end
if ~isfield(Method,'extraOrder')
  Method.extraOrder = Method.order;
end
% number of product states to average over
if ~isfield(System,'nStates')
  System.nStates = ones(1,Method.order);
end
if numel(System.nStates) < Method.order
  System.nStates(end+1:Method.order) = 1;
end

if~isfield(Method,'includeAllSubclusters')
  Method.includeAllSubclusters = false;
end


if ~isfield(Method,'precalculateHamiltonian')
  Method.precalculateHamiltonian = false;
  % If there is not enough memory to save the entire Hamiltonian at once
  % setting precalculateHamiltonian to false may allow the calculation to
  % proceed.
end

% allowing for alternate inputs
if ~isfield(Method,'method')
  Method.method = 'CCE';
end
if strcmp(Method.method,'restrictedCE'),  Method.method = 'rCE';  end
if strcmp(Method.method,'restrictedCCE'),  Method.method = 'rCCE';  end

% The methods rCE and rCE do not use precomputed Hamiltonians.
if strcmp(Method.method,'rCE')||strcmp(Method.method,'rCCE')
  Method.precalculateHamiltonian = false;
end

if ~isfield(Method,'getNuclearContributions')
  Method.getNuclearContributions = false;
end
if ~isfield(Method,'getNuclearSpinContributions')
  Method.getNuclearSpinContributions = false;
end
if ~isfield(Method,'getClusterContributions')
  Method.getClusterContributions = false;
end
if ~isfield(Method,'getUncertainty')
  Method.getUncertainty = false;
end
Method.getContributions = Method.getNuclearContributions || ...
  Method.getNuclearSpinContributions || Method.getClusterContributions ...
  || Method.getUncertainty;
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
System.kg = System.joule*(System.second /System.meter)^2; 

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
System.ge = 2.00231930436256;% https://physics.nist.gov/cgi-bin/cuu/Value?gem 2020-02-08
System.m1H = 1.007825/System.avogadro*System.kg/1000; % J. Emsley, The Elements, Oxford Chemistry Guides (Oxford Univ. Press, New York, NY, 1995).
System.epsilon0 = (8.8541878128e-12)*System.coulomb^2/System.volt/System.meter;
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
  System.Methyl.include = false;
end
if ~isfield(System.Methyl,'method')
  System.Methyl.method = 2;
end
if ~isfield(System.Methyl,'methylMethylCoupling')
  System.Methyl.methylMethylCoupling = false;
end
if ~isfield(System.Methyl,'numberExtraProtons')
  System.Methyl.numberExtraProtons = 0;
end
if ~isfield(Method.cutoff,'methylCoupledOnlyNumber')
  if System.Methyl.method==1
    Method.cutoff.methylCoupledOnlyNumber = 1;
  else
    Method.cutoff.methylCoupledOnlyNumber = 3;
  end
end
if ~isfield(System.Methyl,'moment_of_inertia')
  System.Methyl.moment_of_inertia =  (5.3373e-47)*System.joule*System.second^2; % kg m^2.;
end
if ~isfield(System.Methyl,'V3') && ~isfield(System.Methyl,'tunnel_splitting')
  System.Methyl.V3 = 86*1e-3*System.eV;
end
if ~isfield(System.Methyl,'lockRotors')
  System.Methyl.lockRotors = false;
end
if ~isfield(System.Methyl,'tunnel_splitting')
  System.Methyl.omega_harmonic_oscillator = sqrt(9*System.Methyl.V3/System.Methyl.moment_of_inertia); % rad/s
  
%   System.Methyl.instanton_action = 8*System.Methyl.moment_of_inertia*System.Methyl.omega_harmonic_oscillator/9; % J*s
  System.Methyl.instanton_action = 8*sqrt(System.Methyl.moment_of_inertia*System.Methyl.V3)/3; % J*s
  
%   System.Methyl.K = 4/3*System.Methyl.omega_harmonic_oscillator^(3/2)...
%     *sqrt(System.Methyl.moment_of_inertia/pi/System.hbar); % rad/s
  System.Methyl.K = sqrt(48/pi/System.hbar*sqrt(System.Methyl.V3^3/System.Methyl.moment_of_inertia));
  
  System.Methyl.K = System.Methyl.K/2/pi; % rad/s -> Hz.
  System.Methyl.tunnel_splitting = 3*System.Methyl.K...
    *exp(-System.Methyl.instanton_action/System.hbar); % Hz
  
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

if ~isfield(System,'spinCenter')
  System.spinCenter = 'unknown';
end
% Set pulse sequence.
if ~isfield(System,'experiment')
  System.experiment = 'Hahn';
end
if ~isfield(System,'nPulses') 
  
  if strcmp(System.experiment,'CP_N') || strcmp(System.experiment,'Uhrig_N')
    error(['Please specify the number of pi pulses via "System.nPulses".']) ;
  end
  System.nPulses = "parameter not needed";
 
end


if ~isfield(System,'t0')
  System.t0 = 0;
end
if ~isfield(System,'dt2')
  System.dt2 = 0;
end
if ~isfield(System,'Ndt') || System.Ndt > System.timepoints
  System.Ndt = System.timepoints;
end
% Set up time grid
if System.t0  ~=0
  error('Use of t0 is not recomended.');
end
if isfield(System,'timepoints') && isfield(System,'dt')
  if System.t0 > 0
    System.Time = zeros(1,System.timepoints);
    System.Time(2:end) = ...
      System.t0 + (0:System.dt:(System.timepoints - 2)*System.dt);
  else
    System.Time(1:System.Ndt) = 0:System.dt:(System.Ndt - 1)*System.dt;
    System.Time(System.Ndt+1:System.timepoints) = System.Time(System.Ndt) + (System.dt2:System.dt2:(System.timepoints - System.Ndt)*System.dt2);
  end

elseif isfield(System,'Time')
  System.timepoints = length(System.Time);
  System.dt = abs(System.Time(2) - System.Time(1));
  System.Time = 0:System.dt:(System.timepoints - 1)*System.dt;
  warning('System.Time is not recommended. System.timepoints and System.dt are recommeded instead.');
else
  error('System.timepoints and System.dt are required.');
end

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
  case 'CP_N'
    System.dimensionality = 1;
  case 'Uhrig_N'
    System.dimensionality = 1;
  case 'Hahn-TR' 
    System.dimensionality = 1;
  otherwise
  error('The experiment ''%s'' is not supported.',System.experiment);
end  
if ~isfield(System,'averaging')
  System.averaging = 'powder';
end

if ~isfield(System,'gridSize') || isempty(System.gridSize)
  System.gridSize = 1;
end
if ~isfield(System,'randomOrientation')
  System.randomOrientation = false;
end
if ~isfield(Method,'r_min')
  Method.r_min = 0.1*System.meter*1e-10; % m.
end


% Define theory.
if isfield(System,'Theory')
  System.theory = any(System.Theory);
end

if ~isfield(System,'theory')
  
  if ~isfield(System,'electron_Zeeman')
    System.electron_Zeeman = true;
  end
  
  if ~isfield(System,'nuclear_Zeeman')
    System.nuclear_Zeeman = true;
  end
  
  if ~isfield(System,'hyperfine')
    System.hyperfine = [true true];
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
  if ~isfield(System,'nuclear_quadrupole_filter')
    System.nuclear_quadrupole_filter = ones(3);
  end
  if ~isfield(System,'useMeanField')
    System.useMeanField = false;
  end

  System.theory = [System.electron_Zeeman,...
    System.nuclear_Zeeman,...
    System.hyperfine(1), System.hyperfine(2), ...
    System.nuclear_dipole(1), System.nuclear_dipole(2), ...
    System.nuclear_dipole(3), System.nuclear_dipole(4), ...
    System.nuclear_quadrupole, ...
    System.useMeanField];
else
  
  System.electron_Zeeman    = System.theory(1);
  System.nuclear_Zeeman     = System.theory(2);
  System.hyperfine          = System.theory(3:4);
  System.nuclear_dipole     = System.theory(5:8);
  System.nuclear_quadrupole = System.theory(9);
  System.useMeanField       = System.theory(10);
  
  if ~isfield(System,'nuclear_quadrupole_scale_e2qQh')
    System.nuclear_quadrupole_scale_e2qQh = 1;
  end
  if ~isfield(System,'nuclear_quadrupole_scale_eta')
    System.nuclear_quadrupole_scale_eta = 1;
  end
  if ~isfield(System,'nuclear_quadrupole_filter')
    System.nuclear_quadrupole_filter = ones(3);
  end
end


% Define cluster size specific theories.
if ~isfield(System,'Theory')
  System.Theory = ones(Method.order,length(System.theory)).*System.theory;  
elseif size(System.Theory ,1) < Method.order
  defined_order_ = size(System.Theory ,1);
  System.Theory(defined_order_+1:Method.order,:) = ...
    ones(Method.order-defined_order_,length(System.theory))...
    .*System.Theory(defined_order_,:);
end

System.nStates(~System.Theory(1:length(System.nStates),10)') = 1;

if ~isfield(System,'useThermalEnsemble')
  System.useThermalEnsemble = ~Method.useMultipleBathStates;
end
if Method.useMultipleBathStates && System.useThermalEnsemble
  error('The setting of Method.useMultipleBathStates= useThermalEnsemble = true is not supported.');
end


if ~isfield(System,'solventOnly')
  System.solventOnly = false;
end

if ~isfield(System,'D2O')
  System.D2O = false;
end

if ~isfield(System,'spinHalfOnly')
  System.spinHalfOnly = false;
end
% System limiting options
if ~isfield(System,'limitToSpinHalf')
  System.limitToSpinHalf = Method.reparseNuclei && System.spinHalfOnly;
end
if ~isfield(System,'deuterateProtein')
  System.deuterateProtein = false;
end

if isfield(System,'deuterateAll') && islogical(System.deuterateAll) && System.deuterateAll
  System.deuterateProtein = true;
  System.D2O = true;
end

if ~isfield(System,'defaultExchangability')
  System.defaultExchangability = true;
end
if ~isfield(System,'TMguess')
  % very rough
  System.TMguess = (5+45*System.D2O)*1e-6;
end

if ~isfield(System,'RandomEnsemble') || ~isfield(System.RandomEnsemble,'include')
System.RandomEnsemble.include = false;
end 
if System.RandomEnsemble.include && ~isfield(System.RandomEnsemble,'concentration')
  error('Please specify System.RandomEnsemble.concentration [=] mol/L.');
end
if ~isfield(System.RandomEnsemble,'radius')
  System.RandomEnsemble.radius = System.radius;
end
if ~isfield(System.RandomEnsemble,'innerRadius')
  System.RandomEnsemble.innerRadius = getVanDerWaalsRadius('O');
end
if ~isfield(System.RandomEnsemble,'Type')
  System.RandomEnsemble.Type = 'H';
end
if ~isfield(System.RandomEnsemble,'sphereRadius')
  System.RandomEnsemble.sphereRadius = getVanDerWaalsRadius(System.RandomEnsemble.Type)*System.meter;
end
if ~isfield(System.RandomEnsemble,'Exchangeable')
  System.RandomEnsemble.Exchangeable = true;
end
if ~isfield(System.RandomEnsemble,'isSolvent')
  System.RandomEnsemble.isSolvent = true;
end
if ~isfield(System.RandomEnsemble,'isWater')
  System.RandomEnsemble.isWater = false;
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


if ~isfield(System,'RF')
  System.RF.B1x = 0;
  System.RF.B1y = 0;
  System.RF.nuRF = 0;
  System.RF.use = false;
end
if ~isfield(System.RF, 'B1x')
  System.RF.B1x = 0;
  System.RF.use = true;
end
if ~isfield(System.RF, 'B1y')
  System.RF.B1y = 0;
  System.RF.use = true;
end
if ~isfield(System.RF, 'nuRF')
  System.RF.nuRF = 0;
  System.RF.use = true;
end
if System.RF.use && any(any( System.Theory(:,7:8)))
  error('Only coherence order 0 dipole-dipole coupling is allowed with RF.');
end

if ~isfield(System,'pdbTranslation')
    System.pdbTranslation = [];
end
if ~isfield(System,'pdbRotate')
    System.pdbRotate = false;
end

if System.pdbRotate && ~isfield(System,'pdbOrigin')
    System.pdbOrigin = [0,0,0];
end
if System.pdbRotate && ~isfield(System,'pdbAlpha')
    System.pdbAlpha = 0;
end
if System.pdbRotate && ~isfield(System,'pdbBeta')
    System.pdbBeta = 0;
end
if System.pdbRotate && ~isfield(System,'pdbGamma')
    System.pdbGamma = 0;
end

% save options
Data.path2CluE  = mfilename('fullpath');
Data.path2CluE = Data.path2CluE(1:end-17);
if ~strcmp(Data.path2CluE(end-5:end),[filesep 'CluE' filesep])
  error('Could not establish path to CluE.')
end  
if ~isfield(Data,'OutputData')
  Data.OutputData = '';
end
if ~isfield(Data,'writeSpinPDB')
  Data.writeSpinPDB = true;
end
if ~isfield(Data,'saveLevel')
  Data.saveLevel = 0;
end 
if ~isfield(Data,'ClusterData')
  Data.ClusterData = '';
  Data.exitOnFailedLoad = false;
else
  Data.exitOnFailedLoad = true;
end
if ~isfield(Data,'outPDBoptions')
  Data.outPDBoptions.Honly = false;
end

if strcmp(Method.method,'HD-CCE')
  if isempty(Data.ClusterData)
    error('Please supply a set of clusters.')
  end
  Data.exitOnFailedLoad = true;
  
  if Method.order > 2
    error('HD-CCE is only implemented to 2-clusters.')
  end
  
end

if ~isfield(System,'deuteriumFraction')
  System.deuteriumFraction = 1;
end
if ~isfield(System,'deuteriumFraction_nonExchangeable')
  System.deuteriumFraction_nonExchangeable = System.deuteriumFraction;
end
if ~isfield(System,'newIsotopologuePerOrientation')
  if isempty(Data.ClusterData) && (System.deuteriumFraction >= 1 || System.deuteriumFraction <= 0)
    System.newIsotopologuePerOrientation = false;
  else
    System.newIsotopologuePerOrientation = true;
  end
end
if ~isfield(System,'doPruneNuclei')
  System.doPruneNuclei = false;
end

if ~isfield(System,'HydrogenExchange')
  System.HydrogenExchange = 'OH';
end  

if ~isfield(System,'isUnitCell')
  System.isUnitCell = true;
end

statistics.parameters.radius = System.radius;
statistics.parameters.neighborCutoffCriteria = Method.Criteria;
statistics.parameters.neighborCutoff = Method.cutoff;
statistics.parameters.gridSize = System.gridSize;
end
