% Cluster expansion test script
function testCouplingTensors()

inputList = {'./methyl.pdb','./TEMPO.pdb','./TEMPO.pdb'};
outputList = {'methyl','h18TEMPO','d18TEMPO'};
abundance = [1,1,0];
ePos = { [5,0,0]*1e-10, {28,29}, {28,29}  };
nZeros = [1,4,4];

for ii=2
  generateCouplingCSV(inputList{ii},outputList{ii},...
    abundance(ii),nZeros(ii),ePos{ii} );
end

end

function generateCouplingCSV(inPDB,outStr, Habundance,nZeros,ePos )
oldpath = path;
path('../../',oldpath);
 
Data.InputData = inPDB;

Data.overwriteLevel = 2;

%==========================================================================
% System Settings
%==========================================================================
System.experiment = 'Hahn';
System.gridSize = 1;

% radius from the electron spin to the edge of the system, [m]
System.radius = 1; % m;

% Do not include nitrogen in the Hamiltonian.
System.nitrogen = true;
System.carbon = false;
System.isUnitCell = false;

% number of timepoints
System.nPoints = [2^6 0];

% Time steps sizes
System.dt = [0.1905/2*1e-6 2.25e-6]; % s

% electron coordinate choices
% { n } coordinates of the nth atom from the pdb file
% { m, n } mean coordinates of the mth and nth atoms from the pdb file
% [ x, y, z ] spatial 3-vector
System.Electron.Coordinates = ePos;


% electron g value
System.g = 2.00231930436256;

% applied magnetic field
System.magneticField  = 1.2; % T.
% System.temperature = 20; % K


Method.useCentralSpinSystem = true;
System.particleOptions = {...
  '1H','all', 'abundance', Habundance ...
  %'N','all', 'active', false...
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

%
N = Nuclei.number+1 + nZeros;
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
writetable(T,['Tvertex_',outStr,'.csv']);

Tneighbor = sort(Tneighbor);
T = array2table(Tneighbor);
T.Properties.VariableNames(1) = {'Tneighbor'};
writetable(T,['Tneighbor_',outStr,'.csv']);
end
