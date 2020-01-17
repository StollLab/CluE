% Cluster expansion test script

%==========================================================================
% General Setting
%==========================================================================
clear

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
System.radius = 7e-10; % m; % converges at 1.7 nm, but 0.7 nm shows a reasonable decay curve, but with high TM.
System.inner_radius = 0e-10; % m.

% time points per delay period
System.timepoints = 2^8;%11; %1e3 + 1;
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


System.hyperfine = [true true]; % [zz, zx+zy]
System.nuclear_dipole = [true true true true]; % [A, B, CD, EF]
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
Method.cutoff.dipole = 0;

Method.propagationDomain = 'time-domain';

% verbosity option: true, false
Method.verbose = true;

% parallel computing
Method.parallelComputing = false;

Method.partialSave = false;



%==========================================================================
% Run simulation
%==========================================================================

[SignalMean, twotau, TM_powder] = nuclear_spin_diffusion(System,Method,Data);



%--------------------------------------------------------------------------
%% Plot.
%--------------------------------------------------------------------------
% plot(twotau*1e6,abs(SignalMean/SignalMean(1)),'-','color',[0,1,1]*0.6,'linewidth',1.5);

subplot(2,1,1)
plot(twotau*1e6,real(SignalMean),'-','linewidth',1.5);
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
plot(nu*1e-6,abs(F),'-o','linewidth',1.5);
xlabel('\nu (MHz)');
grid on;  zoom on; 
set(gca,'fontsize',fontsize);