% Cluster expansion test script

%==========================================================================
% General Setting
%==========================================================================
% clear;
clc;
clear;
% Data.ClusterData = "Clusters.mat";
oldpath = path;
path('../',oldpath);
 
Data.InputData = 'TEMPO_GLY_100K_1001.pdb';
% Data.InputData = 'TEMPO_100K_dt100ps_01.pdb'
Data.overwriteLevel = 2;

%==========================================================================
% System Settings
%==========================================================================
System.experiment = 'CPMG-2D';
System.spinCenter = 'TEMPO';
System.gridSize = 1;

% radius from the electron spin to the edge of the system, [m]
System.radius = 14e-10; % m;

% Do not include nitrogen in the Hamiltonian.
System.nitrogen = false;
System.carbon = false;

% Number of time points 
System.nPoints = [2^7 0];

% Time step sizes
System.dt = [5/sum(System.nPoints)*1e-6 2.25e-6]; % s

% electron coordinate choices
% { n } coordinates of the nth atom from the pdb file
% { m, n } mean coordinates of the mth and nth atoms from the pdb file
% [ x, y, z ] spatial 3-vector
System.Electron.Coordinates = {28, 29};
System.X = {28, 29};
System.Y = {1,19};

% electron g value
System.g = [2.0097, 2.0064,2.0025];

% applied magnetic field
System.magneticField  =1.2; % T.
% System.temperature = 20; % K



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

% cluster mehod choices CE, CCE, restrictedCE, restrictedCCE
Method.method = 'LCE';

% maximum cluster size
Method.order = 2;
% Method.extraOrder = 3;
% Method.useInterlacedClusters = [0,0;0,1]>0;
% System.nStates = ones(1,Method.order)*2^10;
% System.nStates(1) = 1;
% System.useThermalEnsemble = true;
% maximum nucleus-nucleus coupling distance;
Method.Criteria = {'dipole'};
Method.cutoff.dipole = [10^3.4,10^3.4,10^3.4]; % Hz

% verbosity option: true, false
Method.verbose = true;

% parallel computing?
Method.parallelComputing = false;

% save each orientation in temporary files until the simm completes?
Method.partialSave = false;

%==========================================================================
%% Run simulation
%==========================================================================

% [SignalMean, twotau, TM_powder,order_n_signals,Nuclei, uncertainty] = CluE(System,Method,Data);

[v, fourtau,~,Order_n_SignalMean] = CluE(System,Method,Data);


%--------------------------------------------------------------------------
%% Plot.
%--------------------------------------------------------------------------
clf
fontsize = 24; 
lw = 1;
tau = fourtau/4*1e6;
fig=pcolor(tau,fliplr(tau),real(flipud(v)));
set(fig, 'EdgeColor', 'none')

cb = colorbar;
cb.Limits=[0 1];

hold on
contour(tau,fliplr(tau),real(flipud(v)), 'color', 'black','linewidth',lw)

plot(tau,tau,'--','color','red','linewidth',lw)

try
  colormap(color_map('no-flip','viridis'))
end

xlabel('\tau_{1} (μs)');
ylabel('\tau_{2} (μs)');
ylabel(cb,'coherence');
pbaspect([1 1 1])


zoom on; 

set(gca,'fontsize',fontsize);

