% Cluster expansion test script

%==========================================================================
% General Setting
%==========================================================================
clear;
clc;

oldpath = path;
path('../',oldpath);
 
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
System.radius = 10e-10; % m; % converges at 1.7 nm, but 0.7 nm shows a reasonable decay curve, but with high TM.
System.inner_radius = 0e-10; % m.

% time points per delay period
System.timepoints = 2^7;%11; %1e3 + 1;
System.nitrogen = true;
%time step size [s]
% System.dt = 5.0e-9; % s.
total_time = 20e-6; % s.
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
System.D2O = false;

System.electron_Zeeman = true;
System.nuclear_Zeeman = true;
System.nuclear_dipole = [true true false false]; % [A, B, CD, EF]
System.hyperfine = [true false]; % [zz, zx+zy]
System.nuclear_quadrupole = true;
%System.nuclear_quadrupole_filter =[0,0,0;0,0,0;0,0,1]; 
System.useMeanField = false;
System.Methyl.include = false;

System.g = [2.0097, 2.0064,2.0025];

System.nStates = [1,1]; 
%==========================================================================
% Method Settings
%==========================================================================

% cluster mehod choices CE, CCE, restrictedCE, restrictedCCE
Method.method = 'CCE';

% maximum cluster size
Method.order = 2;
Method.order_lower_bound = 1;

% maximum nucleus-nucleus coupling distance
% Method.Criteria = {'neighbor','modulation','dipole','minimum-frequency'};
Method.Criteria = {'dipole'};
Method.cutoff.dipole = [0, 0, 0, 0, 0];

Method.propagationDomain = 'time-domain';

% verbosity option: true, false
Method.verbose = false;

% parallel computing
Method.parallelComputing = false;

Method.partialSave = false;


[System, Tensors, Nuclei,Clusters] = setUpSystem(System,Data);

%==========================================================================
% Convergence Options
%==========================================================================
options.converge.radius = true;
options.threshold.radius = 1e-2;
options.delta.radius = 1e-10; % m
options.limit.radius = 20e-10; %m

options.converge.dipole = true;
options.threshold.dipole = 1e-2;
options.delta.dipole = -0.2; 
options.limit.dipole = -9; 
 
options.converge.powder = true;
options.threshold.powder = 1e-2;
options.limit.grid_points = 5810;

options.verbose = true;
options.doPlot = true;
%==========================================================================
% Run Cluster Expansion
%==========================================================================
clf;
hold on;

disp('Calculating convergence.');
[parameters,System,Method,Data] = converge_parameters(System,Method,Data, options);
disp(parameters);
save('parameter.mat','parameters');

Data.OutputData = ['SIM_dh8TEMPO_D2O_converged'];
System.radius = parameters.radius;
Method.cutoff.dipole = parameters.dipole;
System.gridSize = parameters.gridSize;


%==========================================================================
% Run simulation
%==========================================================================
[SignalMean, twotau, TM_powder,order_b_signals,Nuclei] = CluE(System,Method,Data);


%--------------------------------------------------------------------------
%% Plot.
%--------------------------------------------------------------------------
plot(twotau*1e6,abs(SignalMean/SignalMean(1)),'-','color','black','linewidth',2);

plot(twotau*1e6,real(SignalMean),'-','linewidth',1.5);
xlabel('2\tau (\mus)');
ylabel('coherence');
set(gca,'fontsize',12);
grid on;  zoom on; 
fontsize = 24;
set(gca,'fontsize',fontsize);
hold on;
box on