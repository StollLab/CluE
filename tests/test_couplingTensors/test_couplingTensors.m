% Cluster expansion test script

%==========================================================================
% General Setting
%==========================================================================
clear;
clc;
clf;

oldpath = path;
path('../../',oldpath);
 
Data.InputData = './TEMPO.pdb';

Data.overwriteLevel = 2;

%==========================================================================
% System Settings
%==========================================================================
System.experiment = 'Hahn';
System.spinCenter = 'TEMPO';
System.gridSize = 1;

% radius from the electron spin to the edge of the system, [m]
System.radius = 12e-10; % m;

% Do not include nitrogen in the Hamiltonian.
System.nitrogen = true;
System.carbon = false;

% Number of time points
System.nPoints = [2^6 0];

% time step sizes
System.dt = [0.1905/2*1e-6 2.25e-6]; % s

% electron coordinate choices
% { n } coordinates of the nth atom from the pdb file
% { m, n } mean coordinates of the mth and nth atoms from the pdb file
% [ x, y, z ] spatial 3-vector
System.Electron.Coordinates = {28, 29};


% electron g value
System.g = 2.00231930436256;

% applied magnetic field
System.magneticField  = 1.2; % T.
% System.temperature = 20; % K


Method.useCentralSpinSystem = true;
System.particleOptions = {...
  '1H','TEM', 'abundance', 1 ...
  };

% theory selection 
%                  eZ    nZ    HF1   HF2    ddA   ddB  ddCD  ddEF  NQI meanField
System.Theory = [ true, true, true, true, true, true, true, true, true, false;  ... % 1-clusters
                  true, true, true, false, true, true, true, true, true, false; ... % 2-clusters
                  ];

%==========================================================================
% Method Settings
%==========================================================================

% cluster mehod choices CE, CCE, restrictedCE, restrictedCCE
Method.method = 'CCE';

% maximum cluster size
Method.order = 2;
Method.Criteria = {'dipole'};
Method.cutoff.dipole = 0; % Hz

% verbosity option: true, false
Method.verbose = true;

% parallel computing?
Method.parallelComputing = false;

% save each orientation in temporary files until the simm completes?
Method.partialSave = true;

%==========================================================================
%% Make System
%==========================================================================

[System, Method, Data,statistics] = setDefaults(System,Method,Data);
System.pdb = parsePDBfile(Data.InputData, System.angstrom);
[Nuclei, System] = centralSpinSystem(System,Method,Data,System.pdb);


thisCluster = 1:Nuclei.number;

[tensors,zeroIndex] = pairwisetensors_gpu(Nuclei.Nuclear_g, ...
  Nuclei.Coordinates,thisCluster,Nuclei.Atensor,System.magneticField, ...
  System.ge,System.g,...
  System.muB, System.muN, System.mu0, System.hbar,System.theory,...
  0,0,0, ...
  [], []);

%%
N = Nuclei.number+1 + 4;
Tvertex = zeros(N,1);
NN = N*N/2 - N/2;
Tneighbor = zeros(NN,1);
counter = 0;
for iSpin = 0:(Nuclei.number)
  Tvertex(iSpin+1) = tensors(3,3,iSpin+1,iSpin+1);
  for jSpin = (iSpin+1):Nuclei.number
    counter = counter + 1;
    Tneighbor(counter) = tensors(3,3,iSpin+1,jSpin+1);
  end
end
Tvertex = sort(Tvertex); 
T = array2table(Tvertex);
T.Properties.VariableNames(1) = {'vertex'};
writetable(T,'Tvertex.csv');

Tneighbor = sort(Tneighbor);
T = array2table(Tneighbor);
T.Properties.VariableNames(1) = {'Tneighbor'};
writetable(T,'Tneighbor.csv');
%==========================================================================
%% Run simulation
%==========================================================================

% [SignalMean, twotau,~,Order_n_SignalMean] = CluE(System,Method,Data);
% plot(twotau*1e6,real(SignalMean))
% 
