% Cluster expansion test script
%==========================================================================
% General Setting
%==========================================================================
clear;
clc;

oldpath = path;
path('../',oldpath);
 
Data.InputData = 'TEMPO_100K_dt100ps_01.pdb';

%==========================================================================
% System Settings
%==========================================================================
System.experiment = 'Hahn';
System.averaging = 'powder';
System.gridSize = 6;

% radius from the electron spin to the edge of the system, [m]
System.radius = 12e-10; % m; 


System.timepoints = 2^7;
System.dt = 0.1905/2*1e-6; % s


% electron coordinate choices
% { n } coordinates of the nth atom from the pdb file
% { m, n } mean coordinates of the mth and nth atoms from the pdb file
% [ x, y, z ] spatial 3-vector
System.Electron.Coordinates = {28, 29};
System.X = {28, 29};
System.Y = {1,19};
System.magneticField  =1.2; % T.


% deuterium options
%ystem.deuterateProtein = false;
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
System.Theory = [ true, true, true, true, true, true, true, true, true, false;  ... % 1-clusters
                  true, true, true, false, true, true, true, true, true, false; ... % 2-clusters
                  true, true, true, false, true, true, true, true, true, false; ... % 3-clusters
                  true, true, true, false, true, true, true, true, true, false; ... % 4-clusters
                  true, true, true, false, true, true, true, true, true, false; ... % 5-clusters
                  true, true, true, false, true, true, true, true, true, false];    % 6-clusters
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
Method.Ori_cutoffs = false;
Method.cutoff.dipole = 10^3; % Hz

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
clf
fontsize = 24;
if strcmp(Method.method,'count clusters')

plot(abs(SignalMean),'o--','linewidth',3,'color','blue');
xlabel('count');
ylabel('cluster size');
set(gca,'fontsize',12);
grid on;  zoom on; 
fontsize = 24;
set(gca,'fontsize',fontsize);
elseif strcmp(System.experiment,'CPMG-2D')
 time = twotau;
 V = SignalMean;
contour(time/4*1e6,fliplr(time/4)*1e6,real(flipud(V)))
colormap(color_map('-flip','black-body'));
% colormap(color_map('-no-flip','smooth-cool-warm'));
% colormap(color_map('no-flip','smooth-cool-warm'));
cmax = mma(V);
caxis([-cmax cmax])
% imagesc(time,time,V)

lw_2d = 1;
hold on
cb=colorbar;
imagesc(time/4*1e6,fliplr(time/4)*1e6,real(flipud(V)));
contour(time/4*1e6,fliplr(time/4)*1e6,real(flipud(V)),'color',discrete_color_map(8),'linewidth',lw_2d);
plot(time/4*1e6,time/4*1e6,'--','color',discrete_color_map(5),'linewidth',lw_2d);
% plot(twotau1*1e6,twotau2*1e6,'--','color','black','linewidth',lw_2d);

xlabel('\tau_{1} (\mus)');
ylabel('\tau_{2} (\mus)');

%     ylabel(cb,'V(\tau_{1},\tau_{2})/sup_{\tau_{1}}(V(\tau_{1}|\tau_{2}))');
ylabel(cb,'coherence');

set(gca,'fontsize',fontsize);   
else  
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
plot(twotau*1e6,abs(SignalMean),'-','linewidth',3,'color','black');
plot(twotau*1e6,imag(SignalMean),'-','linewidth',1.5,'color','red');
plot(twotau*1e6,real(SignalMean),'-','linewidth',1.5,'color','blue');


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
xlabel('\nu (MHz)');
grid on;  zoom on; 
set(gca,'fontsize',fontsize);
hold on;
end
