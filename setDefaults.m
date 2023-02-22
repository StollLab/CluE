% Supplements the user-provided System, Method and Data structures with default
% values of fields that the user did not provide explicitly.


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [System,Method,Data,statistics] = setDefaults(System,Method,Data)

% Method: define defaults for all fields
defaultMethod.method = 'CCE';
defaultMethod.order = 2;
defaultMethod.errorTolerance = 1e-9;
defaultMethod.useCentralSpinSystem = false;
defaultMethod.conserveMemory = false;
defaultMethod.getNuclearStatistics = false;
defaultMethod.exportHamiltonian = false;
defaultMethod.exportClusters = false;
defaultMethod.propagationDomain='time-domain';
defaultMethod.partialSave = true;
defaultMethod.clear_partialSave = true;
defaultMethod.parallelComputing = false;
defaultMethod.parfor_over_clusters = false;
defaultMethod.numberCores = inf;
defaultMethod.verbose = false;
defaultMethod.allowHDcoupling = false;
defaultMethod.combineClusters = false;
defaultMethod.order_lower_bound = 1;
defaultMethod.shuffle = true;
defaultMethod.graphCriterion = 'connected';
defaultMethod.seed = 42;
defaultMethod.startSpin = 0;
defaultMethod.endSpin = inf;
defaultMethod.record_clusters = false;
defaultMethod.useMultipleBathStates = false;
defaultMethod.useInterlacedClusters = false;
defaultMethod.emptyClusterSetsOkay = false;
defaultMethod.getNuclearContributions = false;
defaultMethod.getNuclearSpinContributions = false;
defaultMethod.getClusterContributions = false;
defaultMethod.getUncertainty = false;
defaultMethod.extraOrder = Method.order;
defaultMethod.includeAllSubclusters = false;
defaultMethod.r_min = 0.1e-10; % m
defaultMethod.mixed_eState = false;
defaultMethod.useMethylPseudoParticles = false;
defaultMethod.fullyConnectedMethyls = false;
defaultMethod.writeClusterStatistics = true;


% Method: Add defaults for fields missing in user-provided structure
Method = supplementdefaults(Method,defaultMethod);

if Method.conserveMemory
  % Auxiliary signals are not saved so getNuclearContributions will fail.
  Method.getNuclearContributions = false;
  Method.getClusterContributions = false;
  Method.getNuclearStatistics = false;
  % The following behaviors are are used regardless of settings,
  % so the settings are modified to reflect what is done.  
  Method.precalculateHamiltonian = false;
  
end

Method = setNeighborCutoffs(Method);

if ~isfield(Method,'vertexCutoff')
  Method.vertexCutoff = struct();
end

if ~isfield(Method.vertexCutoff,'radius_nonSpinHalf') 
  Method.vertexCutoff.radius_nonSpinHalf = System.radius;
end

% Radius to load nuclei in to.
if ~isfield(System,'load_radius')
  System.load_radius = System.radius;
end
System.load_radius = max(System.load_radius, System.radius);

% cluster order max
Method.order_lower_bound = max(1,floor(Method.order_lower_bound));
Method.order_lower_bound= min(Method.order,Method.order_lower_bound);

% Monte Carlo option
if ~isfield(Method,'MonteCarlo')
    Method.MonteCarlo.use = false;
end

if ~isfield(Method.MonteCarlo,'Cluster_Limit')
  Method.MonteCarlo.Cluster_Limit = inf*(1:Method.order); 
elseif length(Method.MonteCarlo.Cluster_Limit) < Method.order
  Cluster_Limit = inf*(1:Method.order);
  Cluster_Limit(1:length(Method.MonteCarlo.Cluster_Limit)) ...
   = Method.MonteCarlo.Cluster_Limit;

  for ii = length(Method.MonteCarlo.Cluster_Limit):Method.order
    Cluster_Limit(ii) ...
     = Method.MonteCarlo.Cluster_Limit(length(Method.MonteCarlo.Cluster_Limit));
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
  Increment(1:length(Method.MonteCarlo.Increment)) ...
   = length(Method.MonteCarlo.Increment);
  for ii = length(Method.MonteCarlo.Increment):Method.order
    Increment(ii) = ...
     Method.MonteCarlo.Increment(length(Method.MonteCarlo.Increment));
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

if (size(Method.useInterlacedClusters,1) ~=...
    size(Method.useInterlacedClusters,2)) ...
    || (size(Method.useInterlacedClusters,1) < Method.order)

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
% number of product states to average over
if ~isfield(System,'nStates')
  System.nStates = ones(1,Method.order);
end
if numel(System.nStates) < Method.order
  System.nStates(end+1:Method.order) = 1;
end

if ~isfield(Method,'precalculateHamiltonian')
  Method.precalculateHamiltonian = false;
  % If there is not enough memory to save the entire Hamiltonian at once
  % setting precalculateHamiltonian to false may allow the calculation to
  % proceed.
end

% The methods rCE and rCE do not use precomputed Hamiltonians.
if strcmp(Method.method,'rCE')||strcmp(Method.method,'rCCE')
  Method.precalculateHamiltonian = false;
end

Method.getContributions = Method.getNuclearContributions || ...
  Method.getNuclearSpinContributions || Method.getClusterContributions ...
  || Method.getUncertainty;

% Fundamental constants
% (see https://physics.nist.gov/cuu/Constants/index.html)
% (all in SI units)
defaultSystem.c = 299792458;  % speed of light in vacuum, m/s
defaultSystem.h = 6.626070040e-34;  % Planck constant, J s
defaultSystem.hbar = 1.054571800e-34;  % reduced Planck constant, J s
defaultSystem.muN = 5.050783699e-27;  % nuclear magneton, J/T
defaultSystem.muB = 927.400e-26;  % Bohr magneton, J/T
defaultSystem.mu0 = 4*pi*1e-7;  % magnetic constant, J^-1 m^3 T^2
defaultSystem.epsilon0 = 8.8541878128e-12;  % electric constant
defaultSystem.kB = 1.38064852e-23; % Boltzmann constant, J/K
defaultSystem.e = 1.6021766208e-19; % elementary charge, C
defaultSystem.eV = 1.6021766208e-19; % electron-volt, J
defaultSystem.barn = 1e-28;  % atomic unit of quadrupole moment, m^2
defaultSystem.angstrom = 1e-10; % m
defaultSystem.avogadro = 6.022140857e23;  % Avogadro number
defaultSystem.ge = 2.00231930436256;  % g value of free electron
defaultSystem.u = 1.66053906660e-27;  % unified atomic mass unit, kg
defaultSystem.m1H = 1.007825*defaultSystem.u; % mass of hydrogen atom, kg

defaultSystem.spinCenter = 'unknown';
defaultSystem.magneticField = 1.2;  % T
defaultSystem.inner_radius = 0;  % m
defaultSystem.temperature = 20 ;  % K
defaultSystem.dt = [];
defaultSystem.nPoints = [];
defaultSystem.solventOnly = false;
defaultSystem.D2O = false;
defaultSystem.spinHalfOnly = false;
defaultSystem.limitToSpinHalf = Method.reparseNuclei && System.spinHalfOnly;
defaultSystem.deuterateProtein = false;
defaultSystem.defaultExchangability = true;
defaultSystem.pdbTranslation = [];
defaultSystem.pdbRotate = false;
defaultSystem.randomOrientation = false;
defaultSystem.TMguess = (5+45)*1e-6;  % s, very rough
defaultSystem.doPruneNuclei = false;
defaultSystem.HydrogenExchange = 'OH';
defaultSystem.isUnitCell = true;
defaultSystem.g = 2.0023*[1,1,1];
defaultSystem.deuteriumFraction = 1;
defaultSystem.CPMG_const_time = 0;

System = supplementdefaults(System,defaultSystem);

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
if ~isfield(Method.neighborCutoff,'methylCoupledOnlyNumber')
  if System.Methyl.method==1
    Method.neighborCutoff.methylCoupledOnlyNumber = 1;
  else
    Method.neighborCutoff.methylCoupledOnlyNumber = 3;
  end
end
if ~isfield(System.Methyl,'moment_of_inertia')
  System.Methyl.moment_of_inertia = 5.3373e-47; % kg m^2
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

% Set electronic spin
if ~isfield(System.Electron,'spin')
  System.Electron.spin = 1/2;
end

% Set g matrix
if ~isfield(System.Electron,'g')
  System.Electron.g = 2.0023;
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

if isempty(System.dt)
  error('Provide a time step in System.dt.');
end
if numel(System.dt)<2
  System.dt(2) = 0;
end
if isempty(System.nPoints)
  error('Provide the number of time points in System.nPoints.');
end
if numel(System.nPoints)<2
  System.nPoints(2) = 0;
end

% Set up time grid
if isfield(System,'nPoints') && isfield(System,'dt')
  N = sum(System.nPoints);
  System.Time(1:System.nPoints(1)) = (0:System.nPoints(1)-1)*System.dt(1);
  if numel(System.nPoints)>1
    System.Time(System.nPoints(1)+1:N) = System.Time(System.nPoints(1)) + (1:N - System.nPoints(1))*System.dt(2);
  end
elseif isfield(System,'Time')
  System.nPoints = length(System.Time);
  System.dt = abs(System.Time(2) - System.Time(1));
  System.Time = (0:System.nPoints-1)*System.dt;
  warning('System.Time is not recommended. System.nPoints and System.dt are recommeded instead.');
end

% Set pulse sequence.
if ~isfield(System,'experiment')
  System.experiment = 'Hahn';
end
if strcmp(System.experiment,'CP_N') || strcmp(System.experiment,'Uhrig_N')
  if ~isfield(System,'nPulses') 
    error('Please specify the number of pi pulses via "System.nPulses".');
  end
else
  System.nPulses = [];
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
    System.Time_ = 0:System.dt_:(System.nPoints - 1)*System.dt_;
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

% Define theory
if isfield(System,'Theory')
  System.theory = any(System.Theory,1);
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

% Define cluster size specific theories
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

if isfield(System,'deuterateAll') && islogical(System.deuterateAll) && System.deuterateAll
  System.deuterateProtein = true;
  System.D2O = true;
end

% Random ensemble
if ~isfield(System,'RandomEnsemble')
  System.RandomEnsemble = struct;
end
if ~isfield(System.RandomEnsemble,'include')
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
  System.RandomEnsemble.sphereRadius = getVanDerWaalsRadius(System.RandomEnsemble.Type);
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

% The methods rCE and rCE do not use precomputed Hamiltonians
if strcmp(Method.method,'rCE') || strcmp(Method.method,'rCCE')
  System.limitToSpinHalf = true;
end

if System.limitToSpinHalf
  disp('Based on the input options, the simulation will only include spin-1/2 nuclei.');
end

if Method.mixed_eState
  
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

% Save options
%-------------------------------------------------------------------------------
Data.path2CluE = fileparts(which(mfilename));
defaultData.OutputData = '';
defaultData.writeSpinPDB = false;
defaultData.saveLevel = 0;
defaultData.outPDBoptions.Honly = false;
defaultData.ClusterData = '';
defaultData.save_mat_file = false;
Data = supplementdefaults(Data,defaultData);

Data.exitOnFailedLoad = ~isempty(Data.ClusterData);

if strcmp(Method.method,'HD-CCE')
  if isempty(Data.ClusterData)
    error('Please supply a set of clusters.')
  end
  Data.exitOnFailedLoad = true;
  if Method.order > 2
    error('HD-CCE is only implemented to 2-clusters.')
  end
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

% Set statistics structure
statistics.parameters.radius = System.radius;
statistics.parameters.neighborCutoffCriteria = Method.Criteria;
statistics.parameters.neighborCutoff = Method.neighborCutoff;
statistics.parameters.gridSize = System.gridSize;

end  % end of main function
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function out = supplementdefaults(user,defaults)
% Supplements the provided user structure with fields and values from
% the structure defaults and returns the result in out.
out = defaults;
fields = fieldnames(user);
for f = 1:numel(fields)
  Name = fields{f};
  out.(Name) = user.(Name);
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function Method = setNeighborCutoffs(Method)

if isfield(Method,'Criteria')
  error(['Error in setMethodCutoffs(): ',...
      'Method.Criteria is deprecated as an external variable.']);
end
Method.Criteria = cell(0);

% Supplement missing cutoff criteria
if ~isfield(Method,'neighborCutoff') 
  Method.neighborCutoff = struct;
end

zer = zeros(1,Method.order);

if ~isfield(Method.neighborCutoff,'bAmax') 
  Method.neighborCutoff.bAmax = zer;
else
   Method.Criteria{end+1} = 'bAmax';
end

if ~isfield(Method.neighborCutoff,'DeltaHyperfine') 
  Method.neighborCutoff.DeltaHyperfine = zer;
else
   Method.Criteria{end+1} = 'delta hyperfine';
end

if ~isfield(Method.neighborCutoff,'dipole') 
  Method.neighborCutoff.dipole = zer;
else
   Method.Criteria{end+1} = 'dipole';
end

if ~isfield(Method.neighborCutoff,'dipoleHalf') 
  Method.neighborCutoff.dipoleHalf = zer;
else
   Method.Criteria{end+1} = 'dipoleHalf';
end

if ~isfield(Method.neighborCutoff,'dipoleOne') 
  Method.neighborCutoff.dipoleOne = zer;
else
   Method.Criteria{end+1} = 'dipoleOne';
end

if ~isfield(Method.neighborCutoff,'modulation') 
  Method.neighborCutoff.modulation = zer;
else
   Method.Criteria{end+1} = 'modulation';
end

if ~isfield(Method.neighborCutoff,'modDepthFreq4') 
  Method.neighborCutoff.modDepthFreq4 = zer;
else
   Method.Criteria{end+1} = 'Hahn: k*omega^4';
end

if ~isfield(Method.neighborCutoff,'Hahn_1us') 
  Method.neighborCutoff.Hahn_1us = zer;
else
   Method.Criteria{end+1} = 'Hahn: k*sin^4(omega*1us)';
end

if ~isfield(Method.neighborCutoff,'rMax') 
  Method.neighborCutoff.rMax = inf*ones(1,Method.order);
else
   Method.Criteria{end+1} = 'rMax';
end

if ~isfield(Method.neighborCutoff,'rMin') 
  Method.neighborCutoff.rMin = zer;
else
   Method.Criteria{end+1} = 'rMin';
end

if ~isfield(Method.neighborCutoff,'methylCoupledOnly')
  Method.neighborCutoff.methylCoupledOnly = false(1,Method.order);
else
   Method.Criteria{end+1} = 'methyl coupled only';
end

if ~isfield(Method.neighborCutoff,'methylOnly')
  Method.neighborCutoff.methylOnly = false(1,Method.order);
else
   Method.Criteria{end+1} = 'methyl only';
end

if ~isfield(Method,'Ori_cutoffs')
  Method.Ori_cutoffs = false;
end
if ~isfield(Method,'reparseNuclei')
  Method.reparseNuclei = false;
end

if numel(Method.neighborCutoff.methylOnly) < Method.order
  n_ = numel(Method.neighborCutoff.methylOnly);
  Method.neighborCutoff.methylOnly(n_:Method.order) = Method.neighborCutoff.methylOnly(n_);
end

if numel(Method.neighborCutoff.methylCoupledOnly) < Method.order
  n_ = numel(Method.neighborCutoff.methylCoupledOnly);
  Method.neighborCutoff.methylCoupledOnly(n_:Method.order) = Method.neighborCutoff.methylCoupledOnly(n_);
end
if numel(Method.neighborCutoff.dipole) < Method.order
  n_ = numel(Method.neighborCutoff.dipole);
  Method.neighborCutoff.dipole(n_:Method.order) = Method.neighborCutoff.dipole(n_);
end
if numel(Method.neighborCutoff.dipoleHalf) < Method.order
  n_ = numel(Method.neighborCutoff.dipoleHalf);
  Method.neighborCutoff.dipoleHalf(n_:Method.order) = Method.neighborCutoff.dipoleHalf(n_);
end
if numel(Method.neighborCutoff.dipoleOne) < Method.order
  n_ = numel(Method.neighborCutoff.dipoleOne);
  Method.neighborCutoff.dipoleOne(n_:Method.order) = Method.neighborCutoff.dipoleOne(n_);
end
if isfield(Method.neighborCutoff, 'maxAmax') && numel(Method.neighborCutoff.maxAmax) < Method.order
  n_ = numel(Method.neighborCutoff.maxAmax);
  Method.neighborCutoff.maxAmax(n_:Method.order) = Method.neighborCutoff.maxAmax(n_);
end
if isfield(Method.neighborCutoff, 'minAmax') && numel(Method.neighborCutoff.minAmax) < Method.order
  n_ = numel(Method.neighborCutoff.minAmax);
  Method.neighborCutoff.minAmax(n_:Method.order) = Method.neighborCutoff.minAmax(n_);
end
if isfield(Method.neighborCutoff, 'minimum_frequency') ...
 && numel(Method.neighborCutoff.minimum_frequency) < Method.order
  n_ = numel(Method.neighborCutoff.minimum_frequency);
  Method.neighborCutoff.minimum_frequency(n_:Method.order) ...
  = Method.neighborCutoff.minimum_frequency(n_);
end
if isfield(Method.neighborCutoff, 'modulation') ...
 && numel(Method.neighborCutoff.modulation) < Method.order
  n_ = numel(Method.neighborCutoff.modulation);
  Method.neighborCutoff.modulation(n_:Method.order) = Method.neighborCutoff.modulation(n_);
end
if isfield(Method.neighborCutoff, 'DeltaHyperfine') ...
 && numel(Method.neighborCutoff.DeltaHyperfine) < Method.order
  n_ = numel(Method.neighborCutoff.DeltaHyperfine);
  Method.neighborCutoff.DeltaHyperfine(n_:Method.order) ...
   = Method.neighborCutoff.DeltaHyperfine(n_);
end

if isfield(Method.neighborCutoff, 'modDepthFreq4') ...
 && numel(Method.neighborCutoff.modDepthFreq4) < Method.order
  n_ = numel(Method.neighborCutoff.modDepthFreq4);
  Method.neighborCutoff.modDepthFreq4(n_:Method.order) ...
   = Method.neighborCutoff.modDepthFreq4(n_);
end
if isfield(Method.neighborCutoff, 'Hahn_1us') ...
 && numel(Method.neighborCutoff.Hahn_1us) < Method.order
  n_ = numel(Method.neighborCutoff.Hahn_1us);
  Method.neighborCutoff.Hahn_1us(n_:Method.order) ...
   = Method.neighborCutoff.Hahn_1us(n_);
end

if numel(Method.neighborCutoff.bAmax) < Method.order
  n_ = numel(Method.neighborCutoff.bAmax);
  Method.neighborCutoff.bAmax(n_:Method.order) = Method.neighborCutoff.bAmax(n_);
end

if norm(Method.neighborCutoff.dipole-Method.neighborCutoff.dipole(1)) > 0
  Method.neighborCutoff.sizeDependent = true;
else
  Method.neighborCutoff.sizeDependent = false;
end
end % of function
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

