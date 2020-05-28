% Cluster expansion test script

%==========================================================================
% General Setting
%==========================================================================
% clear;
clc;

oldpath = path;
path('../',oldpath);
 
Data.InputData = 'TEMPO_100K_dt100ps_01.pdb';
Data.OutputData = '';
Data.saveLevel = 0;

%==========================================================================
% System Settings
%==========================================================================
System.experiment = 'CPMG';
% averaging choices: none, powder, xy
% System.averaging = 'powder';
System.averaging = 'Nitroxide_Wband_Weights';
System.gridSize = 14;

% radius from the electron spin to the edge of the system, [m]
System.radius = 10e-10; % m; % converges at 1.7 nm, but 0.7 nm shows a reasonable decay curve, but with high TM.
System.inner_radius = 0;%e-10; % m.

% time points per delay period
System.timepoints = 2^7;%11; %1e3 + 1;
System.nitrogen = false;
%time step size [s]
% System.dt = 5.0e-9; % s.
total_time = 25e-6; % s.
System.dt = total_time/System.timepoints/2; % s.
%electron coordinate choices
% [ n ] coordinates of the nth atom from the pdb file
% [ m, n ] mean coordinates of the mth and nth atoms from the pdb file
% [ x, y, z ] spatial 3-vector
System.Electron.Coordinates = {28, 29};
System.X = {28, 29};
System.Y = {1,19};
System.magneticField  = 3.37; % T.


% deuterium options
%System.deuterateAll = true;
System.deuterateProtein = false;
System.D2O = false;
% System.deuteriumFraction = 0.5;

%{
System.electron_Zeeman = true;
System.nuclear_Zeeman = true;
System.nuclear_dipole = [true true true true]; % [A, B, CD, EF]
System.hyperfine = [true true]; % [zz, zx+zy]
System.nuclear_quadrupole = true;
System.useMeanField = false;
%}
System.Methyl.include = false;
%                  eZ    nZ    HF1   HF2    ddA   ddB  ddCD  ddEF  NQI meanField
System.Theory = [ true, true, true, true, true, true, true, true, true, false; ... % 1-clusters
                  true, true, true, true, true, true, true, true, true, false];    % 2-clusters
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
Method.cutoff.dipole = 10^4; % Hz

Method.propagationDomain = 'time-domain';

% verbosity option: true, false
Method.verbose = true;

% parallel computing
Method.parallelComputing = false;

Method.partialSave = false;


% [System, Tensors, Nuclei,Clusters] = setUpSystem(System,Data);

%==========================================================================
%% Run simulation
%==========================================================================
Method.exportClusters = false;
% Data.OutputData = 'tempfile';
% Data.ClusterData = 'Clusters.mat';

% System.nuclear_quadrupole_filter = diag([1,1,1]);
[SignalMean, twotau, TM_powder,order_n_signals,Nuclei] = CluE(System,Method,Data);



%--------------------------------------------------------------------------
%% Plot.
%--------------------------------------------------------------------------
% plot(twotau*1e6,abs(SignalMean/SignalMean(1)),'-','color',[0,1,1]*0.6,'linewidth',1.5);
clf
subplot(2,1,1)

hold on
% plot(twotau*1e6,real(SignalMean),'-','linewidth',1.5,'color','blue');
if Method.order>2
  plot(twotau*1e6,real(order_n_signals{2}),'-','linewidth',1.5);
end
if Method.order>3
  plot(twotau*1e6,real(order_n_signals{3}),'-','linewidth',1.5);
end
if Method.order>4
  plot(twotau*1e6,real(order_n_signals{4}),'-','linewidth',1.5);
end
if Method.order>5
   plot(twotau*1e6,real(order_n_signals{5}),'-','linewidth',1.5);
end
plot(twotau*1e6,imag(SignalMean),'-','linewidth',1.5,'color','red');
plot(twotau*1e6,abs(SignalMean),'--','linewidth',3,'color','blue');
xlabel('2\tau (\mus)');
ylabel('coherence');
set(gca,'fontsize',12);
grid on;  zoom on; 
fontsize = 24;
set(gca,'fontsize',fontsize);
hold on;
subplot(2,1,2)
dt = twotau(2)-twotau(1);
nt = size(twotau,2);
nu  = linspace(0,1/dt,nt);
F = fft(SignalMean);
semilogx(nu*1e-6,abs(F),'-o','linewidth',1.5);
xlabel('\nu (MHz)');
grid on;  zoom on; 
set(gca,'fontsize',fontsize);
hold on;