% Cluster expansion test script

%==========================================================================
% General Setting
%==========================================================================
clc
clear
% Data.ClusterData = "Clusters.mat";
oldpath = path;
path('../',oldpath);

Data.InputData = 'TEMPO_GLY_100K_1001.pdb';
% Data.InputData = 'TEMPO_100K_dt100ps_01.pdb'
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

% electron coordinate choices
% { n } coordinates of the nth atom from the pdb file
% { m, n } mean coordinates of the mth and nth atoms from the pdb file
% [ x, y, z ] spatial 3-vector
System.Electron.Coordinates = {28, 29};
System.X = {28, 29};
System.Y = {1,19};

% electron g value
System.g = [2.0097,2.0064,2.0025];

% applied magnetic field
System.magneticField  =1.2; % T.
% System.temperature = 20; % K

% deuterate non-solvents?
System.deuterateProtein = false;

% deuterate solvents?
System.D2O = false;
System.deuteriumFraction = 1;

% Method.useInterlacedClusters = true;
% Method.useMultipleBathStates = true;
% System.nStates = [1,1];
% theory selection 
%                  eZ    nZ    HF1   HF2    ddA   ddB  ddCD  ddEF  NQI meanField
System.Theory = [ true, true, true, true, true, true, true, true, true, false;  ... % 1-clusters
                  true, true, true, true, true, true, true, true, true, false; ... % 2-clusters
                  true, true, true, false, true, true, true, true, true, false; ... % 3-clusters
                  true, true, true, false, true, true, true, true, true, false; ... % 4-clusters
                  true, true, true, false, true, true, true, true, true, false; ... % 5-clusters
                  true, true, true, false, true, true, true, true, true, false];    % 6-clusters

%==========================================================================
% Method Settings
%==========================================================================

% cluster method choices CE, CCE, rCE, rCCE
Method.method = 'CCE';

% maximum cluster size
Method.order = 2;
% Method.extraOrder = 3;
% Method.useInterlacedClusters = [0,0;0,1]>0;
% System.nStates = ones(1,Method.order)*2^10;
% System.nStates(1) = 1;
% System.useThermalEnsemble = true;
% maximum nucleus-nucleus coupling distance;
Method.Criteria = {'dipole'};
Method.neighborCutoff.dipole = [10^3,10^3,10^3]; % Hz

% verbosity option: true, false
Method.verbose = true;

% parallel computing?
Method.parallelComputing = false;

% save each orientation in temporary files until the sim completes
Method.partialSave = true;

%==========================================================================
%% Run simulation
%==========================================================================

[SignalMean, twotau,~,Order_n_SignalMean] = CluE(System,Method,Data);
Data.InputData = 'TEMPO_TIP4P_connect.pdb';
SignalMean0 = CluE(System,Method,Data);

%--------------------------------------------------------------------------
%% Plot
%--------------------------------------------------------------------------
clf
fontsize = 24; 
subplot(2,1,1)

hold on
TM = getTM(twotau*1e6,abs(SignalMean))
TM0  = getTM(twotau*1e6,abs(SignalMean0))
2*(TM-TM0)/(TM+TM0)
plot(twotau*1e6,abs(SignalMean),'-','linewidth',3,'color','red');
plot(twotau*1e6,abs(SignalMean0),'-','linewidth',3,'color','blue');
% plot(twotau*1e6,abs(Order_n_SignalMean{Method.order-1}),'-','linewidth',3,'color','red');
% if max(abs(SignalMean(:) ) )>1
%   set(gca,'yscale','log');
% end
% plot(twotau*1e6,imag(SignalMean),'-','linewidth',1.5,'color','red');
% plot(twotau*1e6,real(SignalMean),'-','linewidth',1.5,'color','blue');
% if ~System.D2O, xlim([0,10]); end

xlabel('2\tau (\mus)');
ylabel('coherence');
set(gca,'fontsize',12);
grid on;  zoom on; 

set(gca,'fontsize',fontsize);

subplot(2,1,2)
dt = twotau(2)-twotau(1);
nt = size(twotau,2);
nu  = linspace(0,1/dt,nt);
F = fft(SignalMean);

semilogx(nu*1e-6,abs(F),'-','linewidth',3,'color',[0,0,0]);
hold on;
semilogx(nu*1e-6,real(F),'-','linewidth',1.5,'color',[0,0,1]);
semilogx(nu*1e-6,imag(F),'-','linewidth',1.5,'color',[1,0,0]);
legend('Abs','Re','Im')
xlabel('\nu (MHz)');
grid on;  zoom on; 
set(gca,'fontsize',fontsize);
hold on;
