% Cluster expansion test script

%==========================================================================
% General Setting
%==========================================================================
clc
clear
oldpath = path;
path('../',oldpath);

Data.InputData = 'TEMPO_GLY_100K_1001.pdb';
Data.OutputData = '';
Data.overwriteLevel = 2;

%==========================================================================
% System Settings
%==========================================================================
System.experiment = 'Hahn';
System.spinCenter = 'TEMPO';
System.gridSize = 1;

% radius from the electron spin to the edge of the system
System.radius = 12e-10; % m

% Whether to include nitrogen or carbon in the Hamiltonian.
System.nitrogen = true;
System.carbon = false;

% Time steps and number of points
System.dt = [0.1905/2, 2.25]*1e-6; % time steps, s
System.nPoints = [64 0]; % number of points

System.Electron.Coordinates = {28, 29};
System.X = {28, 29};
System.Y = {1,19};

% electron g value
System.g = [2.0097,2.0064,2.0025];

% applied magnetic field
System.magneticField  =1.2; % T.


%==========================================================================
% Method Settings
%==========================================================================

% cluster method choices CE, CCE, rCE, rCCE
Method.method = 'CCE';

% maximum cluster size
Method.order = 2;
Method.neighborCutoff.dipole = [10^3,10^3,10^3]; % Hz

% verbosity option: true, false
Method.verbose = true;


%==========================================================================
%% Run simulation
%==========================================================================
Options = struct('System',System,'Method',Method,'Data',Data);
[SignalMean, twotau] = CluE(System,Method,Data);
Options.Data.InputData = 'TEMPO_TIP4P_connect.pdb';
SignalMean0 = CluE(Options);

%--------------------------------------------------------------------------
%% Plot
%--------------------------------------------------------------------------
clf
fontsize = 24; 
% subplot(2,1,1)

hold on

plot(twotau*1e6,abs(SignalMean),'-','linewidth',3,'color','red');
plot(twotau*1e6,abs(SignalMean0),'-','linewidth',3,'color','blue');


xlabel('2\tau (\mus)');
ylabel('coherence');
set(gca,'fontsize',12);
grid on;  zoom on; 

set(gca,'fontsize',fontsize);

% subplot(2,1,2)
% dt = twotau(2)-twotau(1);
% nt = size(twotau,2);
% nu  = linspace(0,1/dt,nt);
% F = fft(SignalMean);
% 
% semilogx(nu*1e-6,abs(F),'-','linewidth',3,'color',[0,0,0]);
% hold on;
% semilogx(nu*1e-6,real(F),'-','linewidth',1.5,'color',[0,0,1]);
% semilogx(nu*1e-6,imag(F),'-','linewidth',1.5,'color',[1,0,0]);
% legend('Abs','Re','Im')
% xlabel('\nu (MHz)');
% grid on;  zoom on; 
% set(gca,'fontsize',fontsize);
% hold on;
