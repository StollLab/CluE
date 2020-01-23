% Cluster expansion test script
%function  [SignalMean, twotau, Signals,TM, TM_powder, Progress] = sim_TEMPO()

%==========================================================================
% General Setting
%==========================================================================
%clear;
% add path to nuclear_spin_diffusion
oldpath = path;
path('../',oldpath);

% path to input pdb
%    Data.InputData = './TEMPO_xyz.pdb';
%  Data.InputData = './TEMPO_TIP4p_10A_conect.pdb';
%Data.InputData = './TEMPO_TIP4P.pdb';
% Data.InputData = 'TEMPO_TIP4P_Rsys_10A_pymol.pdb';
Data.InputData = 'TEMPO_100K_dt100ps_01.pdb';

% Data.InputData = '/home/kudarizaka/uw/stoll/spectral_diffusion/gromacs/TEMPO/rep_01b/TEMPO_TIP4P.pdb';

% name of output file
Data.OutputData = ['SIM_TEMPO'];

% save options :
Data.saveLevel = 1;

%==========================================================================
% System Settings
%==========================================================================

% experiment choics: FID, Hahn, CP
% The restricted cluster exansions only work for Hahn echo simulations.
System.experiment = 'Hahn';

% averaging choices: none, powder, xy
System.averaging = 'powder';
% The gridSize must be in {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230,
% 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,
% 3074, 3470, 3890, 4334, 4802, 5294, 5810}.
System.gridSize = 1 ;

% radius from the electron spin to the edge of the system, [m]
System.radius = 10e-10; % m; % converges at 1.7 nm, but 0.7 nm shows a reasonable decay curve, but with high TM.
System.inner_radius = 0e-10; % m.
% time points per delay period
System.timepoints = 2^10;%11; %1e3 + 1;
System.nitrogen = true;
%time step size [s]
% System.dt = 5.0e-9; % s.
total_time = 15e-6; % s.
System.dt = total_time/System.timepoints/2; % s.
%electron coordinate choices
% [ n ] coordinates of the nth atom from the pdb file
% [ m, n ] mean coordinates of the mth and nth atoms from the pdb file
% [ x, y, z ] spatial 3-vector
System.Electron.Coordinates = {28, 29};
%r_N = [27.360  27.600  28.360]*1e-10; % N
% r_O = [27.140  26.510  29.200]*1e-10; % O
% System.Electron.Coordinates = (r_N + r_O)/2;
% System.Electron.Coordinates = [3.8200   10.5940    7.9120]*1e-10; % m.
% electron g tensor
System.g = 2.0023*[1,1,1];

% magnetic field [T]
System.magneticField  = 1.2; %1.2; % T.

% deuterium options
%System.deuterateAll = true;

System.Methyl.include = false;

System.hyperfine = [true false]; % [zz, zx+zy]
System.nuclear_dipole = [true true false false]; % [A, B, CD, EF]

% System.g = 2.0023*[1,1,1];

%  System.g = [2.007,2.005,2.0078];
% System.g = [2.007,2.0078,2.005];
% System.g = [2.005,2.0078,2.007];
System.g = [2.0097, 2.0064,2.0025];

System.Flip_Angles = pi/180*[90,180];
System.Detection_Operator = [0,1;0,0];

%==========================================================================
% Method Settings
%==========================================================================

% cluster mehod choices CE, CCE, restrictedCE, restrictedCCE
Method.method = 'CCE';
% Method.method = 'count clusters';
% maximum cluster size
% CE <= 4  % CCE any  % restrictedCE = 2 % restrictedCCE = 2
Method.order = 2;
Method.order_lower_bound = 1;
Method.divisions = 'numSpins';
% maximum nucleus-nucleus coupling distance
% Method.Criteria = {'neighbor','modulation','dipole','minimum-frequency'};
Method.Criteria = {'dipole'};
% Method.r0 = 3e-10; % m. converges at 0.9 nm
Method.cutoff.distance = Method.r0;
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
Method.exportHamiltonian = true;

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
%==========================================================================
% Run Cluster Expansion
%==========================================================================
 
% disp('Calculating convergence.');
% [parameters,System,Method,Data] = converge_parameters(System,Method,Data, options);
%  disp(parameters);

disp('Calculating coherence.');
System.deuterateProtein = false;
% tic
[SignalMean, twotau, TM_powder,Order_n_SignalMean,Nuclei]
[SignalMean, twotau, TM_powder] = nuclear_spin_diffusion(System,Method,Data);
out = load(Data.OutputData);
% toc
% [SignalMean_h18, twotau_h18, TM_powder_h18] = nuclear_spin_diffusion(System,Method,Data);
% time_h18 = toc;
% System.deuterateProtein = true;
% [SignalMean_d18, twotau_d18, TM_powder_d18] = nuclear_spin_diffusion(System,Method,Data);
% 
% if false
%   v = cell(6,1);
%   V = zeros(1,System.timepoints);
%   for ii=0:5
%     ang = 180+ii;
%     System.Flip_Angles = pi/180*[90,ang];
%     [SignalMean, twotau, TM_powder] = nuclear_spin_diffusion(System,Method,Data);
%     v{ii+1}=SignalMean;
%     V = V + v{ii+1}*sin(ang*pi/180);
%   end
%   V =V/max(abs(V));
%   SignalMean = V;
% end
%==========================================================================
%% Plot Data
%==========================================================================

doplot = true;
if doplot
  if Method.getNuclearContributions && ~isnan(TM_powder)
    figure();
    %load('SIM_TEMPO_NuclearSpinContribution.mat')
    nSpins =length(NuclearContribution.spin_distance);
    dr = 5.0e-11;
    r = 0:dr:System.radius;
    weights = zeros(size(r));
    for ii = 1:nSpins
      nr = floor(NuclearContribution.spin_distance(ii)/dr);
      weights(nr+1) = weights(nr+1) + NuclearContribution.powder_contribution(ii);
    end
    
    histogram('BinEdges',[r,r(end)+dr]*1e10,'BinCounts',abs(weights),'FaceColor',discrete_color_map(5),'FaceAlpha',1)
    
    title('Spin Distance vs. Importance');
    xlabel('r (angstroms)'); ylabel('blame factor');
    set(gca,'fontsize',12); axis tight; grid on;  zoom on;
    saveas(gcf,'sim_TEMPO_blameFactor_vs_r.png');
  end
  
  if false %~System.Methyl.include %strcmp(Method.propagationDomain,'time-domain')
    figure();
    plot(twotau*1e6,SignalMean,'-','color','black','linewidth',1.5);
    % hold on;
  else
    CloneFig(1,2); hold on;
    % xlim([0,100]);
    % plot(twotau*1e6,SignalMean,'color','black','linewidth',1.5);
%     plot(twotau*1e6,SignalMean*2728/5185,'--','color','red','linewidth',1.5);

%     plot(twotau*1e6,SignalMean,'-','color',discrete_color_map(5),'linewidth',1.5);

%     plot(twotau*1e6,real(SignalMean),'-','color','blue','linewidth',1.5);
%     plot(twotau*1e6,imag(SignalMean),'-','color','red','linewidth',1.5);
%      plot(twotau*1e6,abs(SignalMean),'-','color',[1,0,1]*0.9,'linewidth',1.5);
     plot(twotau*1e6,abs(SignalMean/SignalMean(1)),'-','color',[0,1,1]*0.6,'linewidth',1.5);
%      plot(-twotau*1e6,abs(SignalMean/SignalMean(1)) + 0.00*Method.order,'-','color',[1,0,1]*0.9,'linewidth',1.5);
     
  end
  % xlim([0,10])
  title('TEMPO Coherence Decay');
  xlabel('2\tau (\mus)'); ylabel('v/v_{0}');
  set(gca,'fontsize',12);  
  grid on;  zoom on; %axis tight;
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
  print(['sim_TEMPO_decay'],'-dpng','-r0');
  set(gcf,'WindowStyle','Docked');
end
disp(TM_powder);
%end