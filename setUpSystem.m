
function [System, Tensors, Nuclei,Clusters] = setUpSystem(System,Data)

tic
Method = [];
% set defaults base on specified parameters and for unspecified parameters
[System, Method, Data] = setSystemDefaults(System,Method,Data);

if isfield(System,'Alpha')
  Alpha = System.Alpha;
else
  Alpha = 0;
end
if isfield(System,'Beta')
  Beta = System.Gamma;
else
  Beta = 0;
end  
if isfield(System,'Gamma')
  Gamma = System.Gamma;
else
  Gamma = 0;
end

% Get input data.
if ~isfield(Data,'InputData') || ~ischar(Data.InputData)
  error(['Data.InputData is not a valid file identifier.']);
end
InputData = Data.InputData;
if ~isfile(InputData)
  error(['Could not find the specified input file ', InputData, '.']);
end



% Set verbosity.
verbose = Method.verbose;

if verbose, fprintf('Setting up...\n');  end

if min( (InputData(end-3:end)) == '.pdb')
  Nuclei = parseNuclei(System, Method, InputData);
  System.Electron.Coordinates = [0,0,0];
else
  error('Input data not recognized')
end
Method.order = min(Method.order,Nuclei.number);

if verbose, fprintf('Setup initialized %i nuclei.\n', Nuclei.number); end

% Get rotation matrix from PDB frame to lab frame, via Euler angles
R_pdb2lab = rotateZYZ(Alpha,Beta,Gamma);

% Rotate nuclear coordinates.
Nuclei.Coordinates = Nuclei.Coordinates*R_pdb2lab';

% Rotate nuclear quadrupole tensors.
for inucleus = 1:Nuclei.number  
  Nuclei.Qtensor(:,:,inucleus) = R_pdb2lab*Nuclei.Qtensor(:,:,inucleus)*R_pdb2lab';
  % Elementwise Qtensor manipulation used for testing.  The default filer is ones(3); 
  Nuclei.Qtensor(:,:,inucleus) = Nuclei.Qtensor(:,:,inucleus).*System.nuclear_quadrupole_filter;
end

% Rotate the g-matrix.
if isfield(System,'gFrame')
  % Get rotation matrix.
  g2MolRotation = rotateZYZ(-System.gFrame(1),-System.gFrame(2),-System.gFrame(3));
  
  % Rotate to molecular frame.
  System.gMatrix = g2MolRotation*System.gMatrix_gFrame*g2MolRotation';
end

% Rotate the g-matrix to the lab frame.
System.gMatrix = R_pdb2lab*System.gMatrix_gFrame*R_pdb2lab';


Tensors = pairwisetensors(System,Nuclei,[1:Nuclei.number]);
 
% ========================================================================
% Compile list of connected clusters
% ========================================================================
Clusters = [];
% Loop over cluster sizes, start at the largest (most time consuming) size
for clusterSize = Method.order:-1:1
  
  
  if verbose, fprintf('Finding clusters of size %d.\n', clusterSize); end
  
  Clusters{clusterSize} = findClusters(Nuclei, clusterSize);
  
  % save the number of clusters of each size for export to the user
  Nuclei.numberClusters(clusterSize) = size(Clusters{clusterSize},1);
  
  if verbose
    fprintf('  Found %d clusters of size %d.\n', size(Clusters{clusterSize},1),clusterSize);
  end
  
end


% Determine the output file name.
if isfield(Data,'OutputData') && ~isempty(Data.OutputData)
  OutputData = Data.OutputData;
  if ~strcmp(OutputData(end-3:end),'.mat')
    OutputData = [OutputData, '.mat'];
  end
else
  OutputData = '';
end

% Save if filename is provided.
if ~isempty(OutputData)
  save(OutputData,'System','Tensors','Nuclei','Clusters','InputData');
end
toc
end