% Cluster expansion test script

%==========================================================================
% General Setting
%==========================================================================
clear

Data.InputData = 'TEMPO_100K_dt100ps_01.pdb';
Data.OutputData = '';
Data.saveLevel = 0;

%==========================================================================
% System Settings
%==========================================================================
System.experiment = 'Hahn';
% averaging choices: none, powder, xy
System.averaging = 'powder';
System.gridSize = 1;

% radius from the electron spin to the edge of the system, [m]
System.radius = 7e-10; % m; % converges at 1.7 nm, but 0.7 nm shows a reasonable decay curve, but with high TM.
System.inner_radius = 0e-10; % m.

% time points per delay period
System.timepoints = 2^13;%11; %1e3 + 1;
System.nitrogen = true;
%time step size [s]
% System.dt = 5.0e-9; % s.
total_time = 150e-6; % s.
System.dt = total_time/System.timepoints/2; % s.
%electron coordinate choices
% [ n ] coordinates of the nth atom from the pdb file
% [ m, n ] mean coordinates of the mth and nth atoms from the pdb file
% [ x, y, z ] spatial 3-vector
System.Electron.Coordinates = {28, 29};

System.magneticField  = 1.2; % T.

% deuterium options
%System.deuterateAll = true;
System.deuterateProtein = false;
System.D2O = true;

System.electron_Zeeman = true;
System.nuclear_Zeeman = true;
System.nuclear_dipole = [true true false false]; % [A, B, CD, EF]
System.hyperfine = [true true]; % [zz, zx+zy]
System.nuclear_quadrupole = true;

System.Methyl.include = false;

System.g = [2.0097, 2.0064,2.0025];

System.Flip_Angles = pi/180*[90,180];
System.Detection_Operator = [0,1;0,0];

%==========================================================================
% Method Settings
%==========================================================================

% cluster mehod choices CE, CCE, restrictedCE, restrictedCCE
Method.method = 'CCE';

% maximum cluster size
Method.order = 2;
Method.order_lower_bound = 1;
Method.divisions = 'numSpins';
% maximum nucleus-nucleus coupling distance
% Method.Criteria = {'neighbor','modulation','dipole','minimum-frequency'};
Method.Criteria = {'dipole'};
Method.cutoff.distance = 3e-10;
Method.cutoff.modulation = 0*1e-2;
Method.cutoff.dipole = 1e3;
% Method.cutoff.minimum_frequency = 1/total_time;

% Method.propagationDomain = 'frequency-domain';
Method.propagationDomain = 'time-domain';
% Method.propagationDomain = 'frequency-domain --matlab';

% verbosity option: true, false
Method.verbose = true;

% estimate how much each nucleus contribute to decoherence
% will fail if the coherence never reaches 1/e
% requires CCE or restrictedCCE
Method.getNuclearContributions = false;

% save Hamiltonian as a cell array of 2x2 coupling matrices (bool)
Method.exportHamiltonian = false;

% parallel computing
Method.parallelComputing = false;

% Hamiltonian Type
%Method.HamiltonianType = 'spinTensor';
Method.partialSave = false;
%precalculate Hamiltonian or calculate Hamiltonians clusterwise
Method.precalculateHamiltonian = false;

Method.shuffle = true;

Method.conserveMemory = false;

Method.mixed_eState = false;

Method.vectorized = true;

Method.MonteCarlo.use =false;
Method.MonteCarlo.Cluster_Limit = [inf, 1000, 2000, 2000];
Method.MonteCarlo.Increment = [inf,500,500,500];
Method.MonteCarlo.Threshold =1e-2;
%==========================================================================
% Convergence Options
%==========================================================================
options.converge.radius = true;
options.threshold.radius = 1e-0;
options.delta.radius = 1e-10; % m
options.limit.radius = 12e-10; %m

options.converge.r0 = false;
options.threshold.r0 = 1e-3;
options.delta.r0 = 0.5e-10; % m
options.limit.r0 = 9e-10; %m

options.converge.modulation = false;
options.threshold.modulation = 1e-2;
options.delta.modulation = -0.5;
options.limit.modulation = -9;

options.converge.dipole = true;
options.threshold.dipole = 1e-0;
options.delta.dipole = -0.2;
options.limit.dipole = -9;

% options.converge.minimum_frequency = false;
% options.threshold.minimum_frequency = 1e-3;
% options.delta.minimum_frequency = -0.2;
% options.limit.minimum_frequency = -9;

options.converge.powder = true;
options.threshold.powder = 1e-1;
options.limit.grid_points = 6;

options.verbose = true;

[SignalMean, twotau, TM_powder] = nuclear_spin_diffusion(System,Method,Data);

plot(twotau*1e6,abs(SignalMean/SignalMean(1)),'-','color',[0,1,1]*0.6,'linewidth',1.5);
% xlim([0,10])
xlabel('2\tau (\mus)');
ylabel('v/v_{0}');
set(gca,'fontsize',12);
grid on;  zoom on; %axis tight;
