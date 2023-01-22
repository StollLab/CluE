%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function tests = test_calculate_signal()
tests = functiontests(localfunctions);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function test_2cce(testCase)
Method.Ori_cutoffs = true;

System.experiment = 'Hahn';

% The gridSize must be in {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230,
% 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,
% 3074, 3470, 3890, 4334, 4802, 5294, 5810}.
System.gridSize =  1;

pow = 8;
t_us = 10;

System.nPoints = 100;
System.dt = [t_us]./2./System.nPoints.*1e-6; % s


%==========================================================================
% General Setting
%==========================================================================

OutputData = 'SIM_TEMPO_3CCE';


disp(OutputData)
Data.OutputData = OutputData;

Data.InputData = 'assets/TEMPO_3gly_1npr_center_1001_box50A.pdb';


% save options :
Data.saveLevel = 1;
%==========================================================================
% System Settings
%==========================================================================

% experiment choics: FID, Hahn, CP
% The restricted cluster exansions only work for Hahn echo simulations.

System.spinCenter = 'TEMPO';

System.averaging = 'powder';

System.Electron.Coordinates = {28, 29};
System.X = {28, 29};
System.Y = {1,19};

System.temperature = 20; % K
System.Methyl.include    = true;
System.Methyl.method = 2;
System.magneticField  = 1.2; % T.
System.temperature = 20; % K
kB = 1.38064852e-23;
eV = 1.6021766208e-19;
System.Methyl.include    = true;
%                  eZ    nZ    HF1   HF2    ddA   ddB  ddCD  ddEF  NQI meanField
System.Theory = [ true, true, true, true, true, true, true, true, true, false];

System.g = [2.0097, 2.0064,2.0025];
%==========================================================================
% Method Settings
%==========================================================================
Method.useCentralSpinSystem = true;

% verbosity option: true, false
Method.verbose = true;

% parallel computing
Method.parallelComputing = false;
Method.parfor_over_clusters = false;


%==========================================================================
% Run Cluster Expansion
%==========================================================================
Data.overwriteLevel = 2;
Method.partialSave = false;

Method.order = 3;

System.Methyl.include = true;
Method.useMethylPseudoParticles = true;
  
System.radius = 6e-10; %14e-10;
Method.neighborCutoff.dipole =         10^3.2;
Method.neighborCutoff.DeltaHyperfine = 10^4.4;


nuT = 0*1e3; % Hz

System.Methyl.tunnel_splitting = nuT;
System.particleOptions = {...
  'methyl','TEM', 'tunnelSplitting', 80e3, ...
  'methyl','!TEM', 'tunnelSplitting', System.Methyl.tunnel_splitting, ...
  'nitrogen','all', 'active', true...
  };


Data.OutputData = ['SIM_TEMPO_80_kHz_3gly_1nPrOH_dipole_all_',...
  num2str(Method.order),'CCE_r14A_b_1.6_kHz_DeltaA_25_kHz_ori',...
  num2str(System.gridSize),'_nuT_',...
  num2str(nuT*1e-3), 'kHz'];


Method.use_new_calculate_signals = false;
signal_ref = CluE(System,Method,Data);

Method.use_new_calculate_signals = true;
[signal,t] = CluE(System,Method,Data);

abs_delta_signal = abs(signal - signal_ref);
max_err = max(abs_delta_signal); 
assert(max_err<1e-12);


end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>