% Cluster expansion test script

%==========================================================================
% General Setting
%==========================================================================
% clear

oldpath = path;
path('../',oldpath);

Data.InputData = './TEMPO_100K_dt100ps_01.pdb';
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
System.radius = 8e-10; % m; 
System.inner_radius = 0e-10; % m.

% time points per delay period
System.timepoints = 2^7;%11; %1e3 + 1;
System.nitrogen = true;
%time step size [s]
% System.dt = 5.0e-9; % s.
total_time = 30e-6; % s.
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
Method.verbose = true;

% parallel computing
Method.parallelComputing = false;

Method.partialSave = false;


[System, Tensors, Nuclei,Clusters] = setUpSystem(System,Data);

%==========================================================================
%% Run simulations
%==========================================================================
tic
[SignalMean, twotau, TM_powder,order_b_signals,Nuclei] = nuclear_spin_diffusion(System,Method,Data);
t_CCE = toc


fh = fopen('path_to_spinach.txt');
allLines = textscan(fh,'%s','whitespace','','delimiter','\n');
fclose(fh);
path2spinach = allLines{1}{1};

oldpath = path;

tic
fid = spinach_hahn_echo(System,Data,path2spinach);
t_spinach = toc;

path(oldpath);

fprintf('\nTime for CCE = %d s.\n',t_CCE);
fprintf('Time for spinach = %d s.\n',t_spinach);

%--------------------------------------------------------------------------
%% Plot.
%--------------------------------------------------------------------------
clf;
v_spinach = fid./fid(1);


subplot(2,1,1)
plot(twotau*1e6,abs(SignalMean),'-','linewidth',3, 'color', discrete_color_map(5));
hold on;
plot(twotau*1e6,abs(v_spinach),'--','linewidth',3,'color', discrete_color_map(2));
xlabel('2\tau (\mus)');
ylabel('coherence');
set(gca,'fontsize',12);
grid on;  zoom on; 
fontsize = 24;
set(gca,'fontsize',fontsize);


subplot(2,1,2)
dt = twotau(2)-twotau(1);
nt = size(twotau,2);
nu  = linspace(0,1/dt,nt);
F = fft(SignalMean);
FID = fft(v_spinach);
plot(nu*1e-6,abs(F),'-','linewidth',3, 'color', discrete_color_map(5));
hold on;
plot(nu*1e-6,abs(FID),'--','linewidth',3, 'color', discrete_color_map(2));
% set(gca,'xscale','log');
xlabel('\nu (MHz)');
grid on;  zoom on; 
set(gca,'fontsize',fontsize);
hold on;



