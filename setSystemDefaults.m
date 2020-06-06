
% ========================================================================
% New Function
% ========================================================================

function [System,Method, Data] = setSystemDefaults(System,Method, Data)

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

% Set verbosity.
if ~isfield(Method,'verbose')
  Method.verbose = false;
end

% cutoff criteria
if ~isfield(Method,'order')
  Method.order = 2;
end
if ~isfield(Method,'Criteria') || isempty(Method.Criteria)
  Method.Criteria = {'dipole'};
end
if ~isfield(Method,'cutoff') 
  zer = zeros(1,Method.order);
  Method.cutoff.modulation = zer;
  Method.cutoff.dipole = zer;
  Method.cutoff.max_distance = inf + zer;
  Method.cutoff.min_distance = zer;
  Method.cutoff.hyperfine_sup = inf + zer;
  Method.cutoff.hyperfine_inf = zer;
end

if numel(Method.cutoff.dipole) < Method.order
  n_ = numel(Method.cutoff.dipole);
  Method.cutoff.dipole(n_:Method.order) = Method.cutoff.dipole(n_);
end

% cluster order max
if ~isfield(Method,'order')
  Method.order = 2;
end
if ~isfield(Method,'clusterization')
  Method.clusterization = 'tree-search';
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

% number of product states to average over
if ~isfield(System,'nStates')
  System.nStates = ones(1,Method.order);
end
if numel(System.nStates) < Method.order
  System.nStates(end+1:Method.order) = 1;
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
System.ge = 2.00231930436256;% https://physics.nist.gov/cgi-bin/cuu/Value?gem 2020-02-08
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

if ~isfield(System,'spinCenter')
  System.spinCenter = 'unknown';
end
% Set pulse sequence.
if ~isfield(System,'experiment')
  System.experiment = 'Hahn';
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
  otherwise
  error('The experiment ''%s'' is not supported.',System.experiment);
end  
if ~isfield(System,'averaging')
  System.averaging = 'powder';
end

if ~isfield(System,'gridSize') || isempty(System.gridSize)
  System.gridSize = 1;
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


if ~isfield(System,'useThermalEnsemble')
  System.useThermalEnsemble = ~System.useMeanField;
end
if System.useMeanField && System.useThermalEnsemble
  error('The setting of useMeanField = useThermalEnsemble = true is not supported.');
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
Data.path2CluE  = mfilename('fullpath');
Data.path2CluE = Data.path2CluE(1:end-17);
if ~strcmp(Data.path2CluE(end-5:end),[filesep 'CluE' filesep])
  error('Could not establish path to CluE.')
end  
if ~isfield(Data,'OutputData')
  Data.OutputData = '';
end
if ~isfield(Data,'saveLevel')
  Data.saveLevel = 0;
end 
if ~isfield(Data,'overwriteLevel')
  Data.overwriteLevel = 0;
end
if ~isfield(Data,'ClusterData')
  Data.ClusterData = '';
end

if ~isfield(System,'deuteriumFraction')
  System.deuteriumFraction = 1;
end

if ~isfield(System,'newIsotopologuePerOrientation')
  if isempty(Data.ClusterData) && (System.deuteriumFraction >= 1 || System.deuteriumFraction <= 0)
    System.newIsotopologuePerOrientation = false;
  else
    System.newIsotopologuePerOrientation = true;
  end
end

end
