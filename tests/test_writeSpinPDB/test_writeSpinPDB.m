% function test_writeSpinPDB
% Sam Jahn 2021-04-01
% This script runs two titration series.  One, from a fully deuterated matrix
% to one only fully deuterated on OH bonded hydron,
% while the other replaces CH bonded hydrons with protons.
clear
% System.limitToSpinHalf = true% DO NOT USE
% System.spinHalfOnly = true;
% Cluster expansion test script
% function  sim_TEMPO_IMC()
%%
HW1 = [26.476   3.014   1.911];
HW2 = [25.846   3.554   0.651];
deltaH = HW2-HW1;
rH = norm(deltaH)

% -----------------------------------------------------------------
% 100% Protiation
% powder = 14 .
% r = 1.100000e-09 A. b = 3.981072e+03 Hz. bAmax = -Inf. nOri = 14. 
% TM2 = 5.016431e-06 s.
% eta = 1.917554e-02.
% -----------------------------------------------------------------
rHalf = 16e-10; % m
bHalf = 10^3.6; % Hz.
% -----------------------------------------------------------------
% 100% Deuteration
% powder = 86 .
% r = 1.600000e-09 A. b = 6.309573e+01 Hz. bAmax = -Inf. nOri = 86. 
% TM2 = 1.947080e-04 s.
% eta = 1.369883e-02.
% -----------------------------------------------------------------
rOne = 16e-10; % m
bOne = 10^1.8; % Hz.
%%
%==========================================================================
% General Setting
%==========================================================================
%  clear;
% add path to nuclear_spin_diffusion
oldpath = path;
numCores = feature('numcores');
if numCores > 10
  isHyak = true;
  disp('Running on mox.')
  path('~/app/CluE',oldpath);
  System.gridSize =  86;
  % parallel computing
  useParallelComputing = true;
  Method.parallelComputing = useParallelComputing;
  prename = 'SIMMC_';
  threshold = 5e-3;
  Data.overwriteLevel = 0;
  N = 8;
else
  isHyak = false;
  disp('Running on laptop.')
  path('../../../../CluE',oldpath);
  options.doPlot = true;
  System.gridSize =  1;
  % parallel computing
  useParallelComputing = false;
  Method.parallelComputing = useParallelComputing;
  prename = 'TEST_';
  threshold = 1e-1;
  Data.overwriteLevel = 2;
  N = 4;
end
% path('~/apps/CluE',oldpath);
Data.InputData = '../TEMPO_Gly_70A.pdb';
% Data.InputData = './TEMPO.pdb';
% save options :
Data.saveMore = true;
Data.saveLevel = 1;
Method.getNuclearSpinContributions = false;
%==========================================================================
% System Settings
%==========================================================================

% experiment choics: FID, Hahn, CP
% The restricted cluster exansions only work for Hahn echo simulations.
System.experiment = 'Hahn';
System.spinCenter = 'TEMPO';
System.D2O = true;
% System.deuteriumFraction = 0.01; %1-0.09341/100;
System.deuterateProtein = true;

% averaging choices: none, powder, xy
System.averaging = 'powder';
% The gridSize must be in {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230,
% 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,
% 3074, 3470, 3890, 4334, 4802, 5294, 5810}.


System.randomOrientation = false;
% radius from the electron spin to the edge of the system, [m]
rmin = 16e-10; % m.
System.radius = rmin;
% time points per delay period
System.carbon = false;
System.nitrogen = false;
%time step size [s]
% System.dt = 5.0e-9; % s.
% total_time = 250e-6; % s.
% System.dt = total_time/System.timepoints/2; % s.

pow = 8;
t_us = 10;
T_us = 400;

System.timepoints = 2^pow; %1e3 + 1;
System.Ndt = 2^(pow-1);
System.dt = t_us/2/System.Ndt*1e-6; % s
System.dt2 = (T_us - t_us)/2/(System.timepoints-System.Ndt)*1e-6; % s



%electron coordinate choices
% [ n ] coordinates of the nth atom from the pdb file
% [ m, n ] mean coordinates of the mth and nth atoms from the pdb file
% [ x, y, z ] spatial 3-vector

System.Electron.Coordinates = {28, 29};
System.X = {28, 29};
System.Y = {1,19}; 

%r_N = [27.360  27.600  28.360]*1e-10; % N
% r_O = [27.140  26.510  29.200]*1e-10; % O
% System.Electron.Coordinates = (r_N + r_O)/2;
% System.Electron.Coordinates = [3.8200   10.5940    7.9120]*1e-10; % m.
% electron g tensor

% magnetic field [T]
% System.magneticField  = 0.345; % T.
% System.magneticField  = 3.3; % T.
System.magneticField  = 1.2; % T.
System.temperature = 20; % K

System.Methyl.include    = false;


%                  eZ    nZ    HF1   HF2    ddA   ddB  ddCD  ddEF  NQI meanField
System.Theory = [ true, true, true, true, true, true, true, true, true, false; ... % 1-clusters
                  true, true, true, false, true, true, true, true, true, false];    % 2-clusters

System.g = [2.0097, 2.0064,2.0025];
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

Method.partialSave = true;
%precalculate Hamiltonian or calculate Hamiltonians clusterwise
Method.precalculateHamiltonian = false;

Method.conserveMemory = false;

Method.gpu = false;


Method.Ori_cutoffs = true;
Method.Criteria = {'same g','dipoleHalf','dipoleOne','radius_nonSpinHalf'};%,'spin 1/2 only'};
Method.cutoff.dipoleHalf = bHalf;
Method.cutoff.dipoleOne = bOne;

Method.cutoff.radius_nonSpinHalf = rOne; % m.

%% mixed isotopic composition ================================================

System.randomOrientation = true;
System.newIsotopologuePerOrientation = true;
System.doPruneNuclei = true;
dN = 2;

options.saveEveryN = dN;
Data.OutputData = '';
options.maxN = inf;
options.newProgress = [];

if System.doPruneNuclei
  System.gridSize = 1;
  Method.partialSave = false;
  options.parallelComputing = useParallelComputing; %Method.parallelComputing;
  Method.parallelComputing = false;
  N = 2*numCores;
  dN = 2*numCores;
end

if isHyak
  % allH_frac0 =[0.17357,0.51058,0.7533,0.83085,0.90699,0.93008,0.94211,...
  %   0.95014,0.95428,0.965,0.97339,0.98195,0.99];

   EEfrac = [0.0115    0.5106    0.7533    0.8308    0.9070    0.9301 ...
     0.9421    0.9501    0.9543    0.9650    0.9734    0.9820    0.9900];

  allH_frac = [EEfrac,1-10.^(-2.25:-0.25:-5)];
  %allH_frac = [0.0115    0.5106    0.8308    0.9070];

  OH_frac  = [0.92204 ,0.93163 , 0.95582 , 0.96485 , 0.98055];
  %OH_frac  = [];

  CH_frac =  [0.95272 ,0.96166 , 0.97571 , 0.97878 , 0.98595];
  %CH_frac =  [0.95272 ,0.96166 , 0.97571 , 0.97878];
else
  allH_frac = [0.99]; % [0.0115];
  CH_frac =  [0.95272];
  OH_frac  = [];
end

frac_ = 0.99;

% pdf info
N_water = 7469;
N_glycerol = 1500;

N_OH = 2*N_water + 3*N_glycerol;
N_CH = 5*N_glycerol;
N_H = N_OH + N_CH;

xOH = N_OH/N_H;
xCH = N_CH/N_H;

% ENUM
CH = 1;
OH = 2;
ALL = 3;


Method.useCentralSpinSystem = true;
Method.reparseNuclei = Method.useCentralSpinSystem;
Method.conserveMemory = true;

Data.outPDBoptions.Honly = true;


%{
for ii = 1:numel(CH_frac)
  
  %Method.cutoff.dipoleOne = bOne*frac_; 
  %Method.cutoff.dipoleHalf = bHalf*Hfrac_;
  %System.radius = max(rmin, rHalf*Hfrac_^(-1/3));
  
  % frac_tot*N_H = frac_OH*N_OH + frac_CH*N_CH
  % frac_OH = (frac_tot*N_H -frac_CH*N_CH) /N_OH.   
  System.deuteriumFraction = frac_;
  System.deuteriumFraction_nonExchangeable = ...
  (CH_frac(ii)*N_H -frac_*N_OH)/N_CH;
  
  
  fracOH = 1 - System.deuteriumFraction;
  fracCH = 1 - System.deuteriumFraction_nonExchangeable;
    
  System.particleOptions = {...
      'nitrogen','TEM', 'active', false, ...
      'hydrogen','TEM', 'abundance', 0, ...
      'hydrogen','TEM', 'active', false, ...
      'hydrogen','all','switchParticle', 'void',...
      'hydrogen','all','extraCellSwitchParticle', 'void',...
      '1H_exchangeable','all', 'abundance', fracOH, ... % all exchangeable hydrons
      '1H_nonExchangeable','!TEM', 'abundance', fracCH ... % all nonExchangeable hydron except those on except TEMPO
      };

  N_ = fracOH*N_OH + fracCH*N_CH;
  r_ = rHalf*(N_H/N_)^(1/3);

  System.radius = max(rmin, r_);
  Method.cutoff.dipoleOne = bOne*(...
    System.deuteriumFraction*N_OH ...
    + System.deuteriumFraction_nonExchangeable*N_CH...
    )/N_H;
  Method.cutoff.dipoleHalf = bHalf*(N_/N_H);

  
  savefile = [prename, 'Honly_d18TEMPO_', ...
    'all', num2str(1000-N_/N_H*1000),'pm_', ...
    'OH_', num2str(System.deuteriumFraction*1000),'pm_', ...
    'CH_', num2str(System.deuteriumFraction_nonExchangeable*1000) , 'm.mat'];

  disp(savefile);

  isotopeMonteCarlo(System,Method,Data, savefile,N,dN, threshold,options);
  
end

for ii = 1:numel(OH_frac)
  System.deuteriumFraction = ...
  (OH_frac(ii)*N_H -frac_*N_CH)/N_OH;

  System.deuteriumFraction_nonExchangeable = frac_;
  
  fracOH = 1 - System.deuteriumFraction;
  fracCH = 1 - System.deuteriumFraction_nonExchangeable;

  System.particleOptions = {...
    'nitrogen','TEM', 'active', false, ...
    'hydrogen','TEM', 'abundance', 0, ...
    'hydrogen','TEM', 'active', false, ...
    'hydrogen','all','switchParticle', 'void',...
    'hydrogen','all','extraCellSwitchParticle', 'void',...
    '1H_exchangeable','all', 'abundance', fracOH, ... % all exchangeable hydrons
    '1H_nonExchangeable','!TEM', 'abundance', fracCH ... % all nonExchangeable hydron except those on except TEMPO
    };
  N_ = fracOH*N_OH + fracCH*N_CH;
  r_ = rHalf*(N_H/N_)^(1/3);

  System.radius = max(rmin, r_);
  Method.cutoff.dipoleOne = bOne*(...
    System.deuteriumFraction*N_OH ...
    + System.deuteriumFraction_nonExchangeable*N_CH...
    )/N_H;
  Method.cutoff.dipoleHalf = bHalf*(N_/N_H);

  savefile = [prename, 'Honly_d18TEMPO_', ...
    'all', num2str(1000 - N_/N_H*1000),'pm_', ...
    'OH_', num2str(System.deuteriumFraction*1000),'pm_', ...
    'CH_', num2str(System.deuteriumFraction_nonExchangeable*1000) , 'm.mat'];

    disp(savefile);

    isotopeMonteCarlo(System,Method,Data, savefile,N,dN, threshold,options);
  
end
%}
for ii = 1:numel(allH_frac)
  System.deuteriumFraction = allH_frac(ii);

  System.deuteriumFraction_nonExchangeable = allH_frac(ii);

  fracOH = 1 - System.deuteriumFraction;
  fracCH = 1 - System.deuteriumFraction_nonExchangeable;

  System.particleOptions = {...
    'nitrogen','TEM', 'active', false, ...
    'hydrogen','TEM', 'abundance', 0, ...
    'hydrogen','TEM', 'active', false, ...
    'hydrogen','all','switchParticle', 'void',...
    'hydrogen','all','extraCellSwitchParticle', 'void',...
    '1H_exchangeable','all', 'abundance', fracOH, ... % all exchangeable hydrons
    '1H_nonExchangeable','!TEM', 'abundance', fracCH ... % all nonExchangeable hydron except those on except TEMPO
    };
  N_ = fracOH*N_OH + fracCH*N_CH;
  r_ = rHalf*(N_H/N_)^(1/3);

  System.radius = max(rmin, r_);
  Method.cutoff.dipoleOne = bOne*(...
    System.deuteriumFraction*N_OH ...
    + System.deuteriumFraction_nonExchangeable*N_CH...
    )/N_H;
  Method.cutoff.dipoleHalf = bHalf*(N_/N_H);

  
  Data.OutputData = [prename, 'Honly_d18TEMPO_', ...
    'all', num2str(1000 - N_/N_H*1000),'pm_', ...
    'OH_', num2str(System.deuteriumFraction*1000),'pm_', ...
    'CH_', num2str(System.deuteriumFraction_nonExchangeable*1000) , 'm.mat'];

    disp(Data.OutputData);
    tic
    [System, Method, Data,statistics] = setSystemDefaults(System,Method,Data);
    pdb = parsePDBfile(Data.InputData,System.angstrom);
    [Nuclei, System]= centralSpinSystem(System,Method,Data,pdb);
    toc
end
% exit();







% end
