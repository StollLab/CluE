% Read PDB file
function pdbData = parsePDBfile(filename, angstrom)

if nargin < 2
  angstrom = 1e-10; % meters.
end

pdbData = struct( ...
  ... % Crystal Information
  'a', -1, ...
  'b', -1, ...
  'c', -1, ...
  'alpha', 90, ...
  'beta', 90, ...
  'gamma', 90, ...
  'sGroup', '', ...
  'Z', -1, ...
  ... %s
  ... % Atom Information
  'number', 0, ...                                          
  'serial', []);

  pdbData.name = {};
  pdbData.altLoc =  {};                                      
  pdbData.resName =  {};                               
  pdbData.chainID =  {};
  pdbData.resSeq =  [];
  pdbData.iCode =  {};
  pdbData.x =  [];
  pdbData.y =  [];
  pdbData.z =  [];
  pdbData.occupancy = [];
  pdbData.tempFactor = [];
  pdbData.element =  {};
  pdbData.charge =  {};
  %
  % Conectivity Information
  pdbData.connections =  sparse([]);
  pdbData.maxConnections =  0;


% Open pdb file.
fh = fopen(filename);

% Load in content of pdb file.
allLines = textscan(fh,'%s','whitespace','','delimiter','\n');

% Close pdb file.
fclose(fh);

% Remove extra lat=yer of cell arrays.
allLines = allLines{1};

% Count number of lines in pdb.
nLines = numel(allLines);


% Loop over lines of the pdb to count the number of atoms.
for iline = 1:nLines
  
  % Get the ith line.
  line_ = allLines{iline};
  
  if strncmp(line_,'ATOM',4) || strncmp(line_,'HETATM',6)
    % Update count.
    pdbData.number = pdbData.number + 1;
  end
  
end

% Initialize variables.
pdbData.serial      = zeros(pdbData.number,1);
pdbData.name        = cell(pdbData.number,1);
pdbData.altLoc      = cell(pdbData.number,1);
pdbData.resName     = cell(pdbData.number,1);
pdbData.chainID     = cell(pdbData.number,1);
pdbData.resSeq      = zeros(pdbData.number,1);
pdbData.iCode       = cell(pdbData.number,1);
pdbData.x           = zeros(pdbData.number,1);
pdbData.y           = zeros(pdbData.number,1);
pdbData.z           = zeros(pdbData.number,1);
pdbData.occupancy   = zeros(pdbData.number,1);
pdbData.tempFactor  = zeros(pdbData.number,1);
pdbData.element     = cell(pdbData.number,1);
pdbData.charge      = cell(pdbData.number,1);

pdbData.connections = sparse(pdbData.number,pdbData.number);

n = 0;

% Loop over lines of the pdb.
for iline = 1:nLines
  
  % Get the ith line.
  line_ = allLines{iline};
  
  % Check if line is a atom specifier.
  if strncmp(line_,'CRYST1',6)
    % Parse information about unit cell

    pdbData.a = str2double(line_(7:15))*angstrom; % angstrom -> m.
    pdbData.b = str2double(line_(16:24))*angstrom; % angstrom -> m.
    pdbData.c = str2double(line_(25:33))*angstrom; % angstrom -> m.
    
    pdbData.alpha = str2double(line_(34:40))*pi/180; % degree -> rad.
    pdbData.beta = str2double(line_(41:47))*pi/180; % degree -> rad.
    pdbData.gamma = str2double(line_(48:54))*pi/180; % degree -> rad.
    try
    pdbData.sGroup = strtrim( line_(56:66) );
    catch
      pdbData.sGroup = 1;
      disp('Could not find space group info in pdb file, setting sGroup to 1.');
    end
    try
      pdbData.Z = str2double( line_(67:70) );
    catch
      disp('Could not find the pdb number polymeric chains in a unit cell.');
      pdbData.Z = -1;
    end
    
  elseif strncmp(line_,'ATOM',4) || strncmp(line_,'HETATM',6)
    
    % Update count.
    n = n + 1;    
    
    % Parse information about atom.   
    pdbData.serial(n)     = str2double(line_(7:11));
    pdbData.name{n}       = strtrim( line_(13:16) );
    pdbData.altLoc{n}     = strtrim( line_(17) );
    pdbData.resName{n}    = strtrim( line_(18:20) );
    if isempty(pdbData.resName{n})
      pdbData.resName{n} = 'NoResName';
    end
    pdbData.chainID{n}    = strtrim( line_(22) );
    pdbData.resSeq(n)     = str2double(line_(23:26));
    pdbData.iCode{n}      = strtrim( line_(27) );
    pdbData.x(n)          = str2double(line_(31:38))*angstrom;
    pdbData.y(n)          = str2double(line_(39:46))*angstrom;
    pdbData.z(n)          = str2double(line_(47:54))*angstrom;
    pdbData.occupancy(n)  = str2double(line_(55:60));
    pdbData.tempFactor(n) = str2double(line_(61:66));
    pdbData.element{n}    = strtrim( line_(77:78) );
    pdbData.charge{n}     = strtrim( line_(79:80) );
    
  % Check if line contains connection data.  
  elseif strncmp(line_,'CONECT',6)
    % Initialize connetion array.
    if isempty(pdbData.connections)
      pdbData.connections = cell(1,iNucleus);
    end
    if any(line_(7:end)=='*')
      continue
    end
    
    currentLine = strtrim(line_);
    currentLine = reshape(currentLine(7:end),5,[]).';
    
    numConnect  = size(currentLine,1);
    if numConnect < 2
      continue;
    end
    referenceNucleus = str2double(currentLine(1,:));
    

    for ii = 2:numConnect
      index = str2double(currentLine(ii,:));
      pdbData.connections(index,referenceNucleus) = 1;
      pdbData.connections(referenceNucleus,index) = 1;
    end
     
  end
  
end


if strcmp(pdbData.sGroup,'P -1')
  addMirrorPDB();
end

return;

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function addMirrorPDB()

  pdbData.serial      = [pdbData.serial; pdbData.serial + pdbData.serial(end)];
  pdbData.name        = [pdbData.name(:);pdbData.name(:)];
  pdbData.altLoc      = [pdbData.altLoc(:);pdbData.altLoc(:)];
  pdbData.resName     = [pdbData.resName(:);pdbData.resName(:)];
  pdbData.chainID     = [pdbData.chainID(:); pdbData.chainID(:)];
  pdbData.resSeq      = [pdbData.resSeq; pdbData.resSeq + pdbData.resSeq(end)];
  pdbData.iCode       = [pdbData.iCode(:);pdbData.iCode(:)];
  pdbData.x           = [pdbData.x; -pdbData.x];
  pdbData.y           = [pdbData.y; -pdbData.y];
  pdbData.z           = [pdbData.z; -pdbData.z];
  pdbData.occupancy   = [pdbData.occupancy; pdbData.occupancy];
  pdbData.tempFactor  = [pdbData.occupancy; pdbData.tempFactor];
  pdbData.element     = [pdbData.element(:); pdbData.element(:)];
  pdbData.charge      = [pdbData.charge(:); pdbData.charge(:)];
  
  N = pdbData.number;
  pdbData.number = 2*N;
  pdbData.connections(N+1:2*N,N+1:2*N) = pdbData.connections;

  
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end






