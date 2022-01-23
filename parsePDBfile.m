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
    pdbData.serial(n)     = sscanf(line_(7:11),'%u'); 
    pdbData.name{n}       = strtrim( line_(13:16) );
    pdbData.altLoc{n}     = strtrim( line_(17) );
    pdbData.resName{n}    = strtrim( line_(18:20) );
    if isempty(pdbData.resName{n})
      pdbData.resName{n} = 'NoResName';
    end
    pdbData.chainID{n}    = strtrim( line_(22) );
    pdbData.resSeq(n)     = str2double(line_(23:26));
    pdbData.iCode{n}      = strtrim( line_(27) );
    pdbData.x(n)          = sscanf(line_(31:38),'%f')*angstrom;
    pdbData.y(n)          = sscanf(line_(39:46),'%f')*angstrom;
    pdbData.z(n)          = sscanf(line_(47:54),'%f')*angstrom;
    pdbData.occupancy(n)  = sscanf(line_(55:60),'%f');
    pdbData.tempFactor(n) = sscanf(line_(61:66),'%f');
    pdbData.element{n}    = strtrim( line_(77:78) );
    pdbData.charge{n}     = strtrim( line_(79:80) );
    
  % Check if line contains connection data.  
  elseif strncmp(line_,'CONECT',6)
    % Initialize connection array.
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

    referenceNucleus =  sscanf(line_(7:11),'%u');

    for ii = 2:numConnect
      index  =  sscanf(line_(2+5*ii:6+5*ii),'%u'); 
      pdbData.connections(index,referenceNucleus) = 1;
    end

  end
  
end

testConnectios = pdbData.connections - pdbData.connections';
if max(abs(testConnectios(:) ))>0
  error(['Error in parsePDBfile(): ', ...
    'connections matrix is not symmetric.']);
end

if strcmp(pdbData.sGroup,'P -1')
  addMirrorPDB();
elseif strcmp(pdbData.sGroup,'R3c') || strcmp(pdbData.sGroup,'H 3 c')
  addR3c()
end

pdbData.number_resSeq = numel(unique(pdbData.resSeq));

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
  pdbData.tempFactor  = [pdbData.tempFactor; pdbData.tempFactor];
  pdbData.element     = [pdbData.element(:); pdbData.element(:)];
  pdbData.charge      = [pdbData.charge(:); pdbData.charge(:)];
  
  N = pdbData.number;
  pdbData.number = 2*N;
  pdbData.connections(N+1:2*N,N+1:2*N) = pdbData.connections;

  
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function addR3c() 
  Ncopy = 6;
  NNcopy = [Ncopy,1];
  N = pdbData.number;
  pdbData.number = Ncopy*N;
  connections = pdbData.connections;
  pdbData.connections = sparse(pdbData.number,pdbData.number);
  for icon = 1:Ncopy
    pdbData.connections((icon-1)*N+1:icon*N,(icon-1)*N+1:icon*N) = connections;
  end

  serialN = pdbData.serial(end);
  resN = pdbData.resSeq(end);
  pdbData.serial      = repmat(pdbData.serial,NNcopy);
  pdbData.serial      = pdbData.serial ...
    +serialN.*floor([0:numel(pdbData.serial)-1]'./serialN);
  pdbData.name        = repmat(pdbData.name(:),NNcopy);
  pdbData.altLoc      = repmat(pdbData.altLoc(:), NNcopy);
  pdbData.resName     = repmat(pdbData.resName(:), NNcopy);
  pdbData.chainID     = repmat(pdbData.chainID(:), NNcopy);
  pdbData.resSeq      = repmat(pdbData.resSeq,NNcopy);
  pdbData.resSeq      = pdbData.resSeq ...
    + resN.*floor([0:numel(pdbData.resSeq)-1]'./resN);
  pdbData.iCode       = repmat(pdbData.iCode(:), NNcopy);
  pdbData.occupancy   = repmat(pdbData.occupancy, NNcopy);
  pdbData.tempFactor  = repmat(pdbData.tempFactor, NNcopy);
  pdbData.element     = repmat(pdbData.element(:), NNcopy);
  pdbData.charge      = repmat(pdbData.charge(:), NNcopy);
  
%  1>   x           , y           , z
%  2>   1/3 + x     , 2/3 + y     , 2/3 + z
%  3>   2/3 + x     , 1/3 + y     , 1/3 + z
%  4>   -y          , x - y       , z
%  5>   -x + y      , -x          , z
%  6>   1/3 - y     , 2/3 + x - y , 2/3 + z
%  7>   1/3 - x + y , 2/3 - x     , 2/3 + z
%  8>   2/3 - y     , 1/3 + x - y , 1/3 + z
%  9>   2/3 - x + y , 1/3 - x     , 1/3 + z
% 10>   -y          , -x          , 1/2 + z
% 11>   1/3 - y     , 2/3 - x     , 1/6 + z
% 12>   2/3 - y     , 1/3 - x     , 5/6 + z
% 13>   -x + y      , y           , 1/2 + z
% 14>   x           , x - y       , 1/2 + z
% 15>   1/3 - x + y , 2/3 + y     , 1/6 + z
% 16>   1/3 + x     , 2/3 + x - y , 1/6 + z
% 17>   2/3 - x+y   , 1/3 + y     , 5/6 + z
% 18>   2/3 + x     , 1/3 + x - y , 5/6 + z

% TO DO: DOUBLE CHECK x0, y0, AND z0 ASSIGNMENTS,
x0 = pdbData.a;
x = pdbData.x/x0;

y0 = pdbData.b;
y = pdbData.y/y0;
z0 = pdbData.c;
z = pdbData.z/z0;

% pdbData.x = [       x; 1/3+x;   2/3+x;    -y;  -x+y; 1/3-y; ...
%               1/3-x+y; 2/3-y; 2/3-x+y;    -y; 1/3-y; 2/3-y; ...    
%                  -x+y    ; x; 1/3-x+y; 1/3+x; 2/3-x+y; 2/3+x].*x0;  

% pdbData.x = [       x; 1/3+x;   2/3+x;    -y;    -x+y; 1/3-y;...
%               1/3-x+y; 2/3-y; 2/3-x+y;    -y;   1/3-y; 2/3-y;...
%                  -x+y;     x; 1/3-x+y; 1/3+x; 2/3-x+y; 2/3+x].*x0;
% 
% pdbData.y = [      y;  2/3+y; 1/3+y;     x-y;    -x; 2/3+x-y;...
%               2/3-x; 1/3+x-y; 1/3-x;      -x; 2/3-x; 1/3-x;...
%                   y;     x-y; 2/3+y; 2/3+x-y; 1/3+y; 1/3+x-y].*y0;
% 
% pdbData.y = [     y;   2/3+y; 1/3+y;     x-y;    -x; 2/3+x-y;...
%               2/3-x; 1/3+x-y; 1/3-x;      -x; 2/3-x;   1/3-x;...
%                   y;     x-y; 2/3+y; 2/3+x-y; 1/3+y; 1/3+x-y;].*y0;

% pdbData.z = [     z; 2/3+z; 1/3+z;     z;     z; 2/3+z;...
%               2/3+z; 1/3+z; 1/3+z; 1/2+z; 1/6+z; 5/6+z;...
%               1/2+z; 1/2+z; 1/6+z; 1/6+z; 5/6+z; 5/6+z].*z0;

% pdbData.z = [     z; 2/3+z; 1/3+z;     z;     z; 2/3+z;...
%               2/3+z; 1/3+z; 1/3+z; 1/2+z; 1/6+z; 5/6+z;...
%               1/2+z; 1/2+z; 1/6+z; 1/6+z; 5/6+z; 5/6+z].*z0;


pdbData.x = [ x; z; y; 1/2+x; 1/2+z; 1/2+y].*x0;
pdbData.y = [ y; x; z; 1/2+z; 1/2+y; 1/2+x].*y0;
pdbData.z = [ z; y; x; 1/2+y; 1/2+x; 1/2+z].*z0;
% 
% pdbData.x = [pdbData.x;-pdbData.x];
% pdbData.y = [pdbData.y;-pdbData.y];
% pdbData.z = [pdbData.z;-pdbData.z];

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end






