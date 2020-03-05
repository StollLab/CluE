% Read PDB file
%
% Input:
%   filename   PDB file name
%   System     structure with field
% Output:
%   Coordinates Nx3 array of coordinates (meters)
%   Type        cell array of character string with isotope information
%   UnitCell
%   Connected   cell array of vectors containing connectivity information
%   Indices_nonWater  list of nuclei that are not part of a water molecule
%   pbdID       pdb ID
%   numberH     [nProtons nDeuterons nHydrogensTotal]

function [Coordinates,Type,UnitCell,Connected, Indices_nonWater, pdbID,numberH, isSolvent] = parsePDB(filename,System)

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

% Initialize variables.
Connected = {};

isSolvent = false(1,nLines);
Coordinates = zeros(nLines,3);
pdbID = zeros(1,nLines);
Type = cell(1,nLines);

UnitCell.isUnitCell = false;
numberH = [0,0,0];
iNucleus = 0;

% Loop over lines of the pdb.
for iline = 1:nLines
  
  % Get the ith line.
  line_ = allLines{iline};
  
  % Check if line is a atom specifier.
  if strncmp(line_,'ATOM',4) || strncmp(line_,'HETATM',6)
    
    % Parse information about atom.
    
    % Update count.
    iNucleus = iNucleus + 1;
    
    % Get pdb ID #.
    pdbID(iNucleus) = sscanf(line_(7:12),'%f');
    
    % Get coordinates in meters.
    Coordinates(iNucleus,:) = sscanf(line_(31:54),'%f %f %f')*System.angstrom; % angstrom -> m
    
    % Get name.
    ResidueName_ = line_(18:20);
    
    % Get element.
    Element_ = strtrim(line_(77:78));
    
    % Determine if atom is not part of the solvent.
    isSolvent(iNucleus) = true;
    if ~strcmp(ResidueName_,'WAT')  && ~strcmp(ResidueName_,'SOL')
      isSolvent(iNucleus) = false;
    end
    
    % Determine if the atom is hydrogen.
    if strcmp(Element_,'H')
      
      % Replace H with D depending on options.
      if strcmp(ResidueName_,'WAT') || strcmp(ResidueName_,'SOL')
        if System.D2O, Element_ = 'D'; end
      else
        if System.deuterateProtein, Element_ = 'D'; end
        if System.solventOnly, Element_ = 'null'; end
      end
    end
    
    % Count hydrons.
    if strcmp(Element_,'H')
      numberH(1) = numberH(1) + 1;
      numberH(3) = numberH(3) + 1;
    elseif strcmp(Element_,'D')
      numberH(2) = numberH(2) + 1;
      numberH(3) = numberH(3) + 1;
    end
    
    % Record elemental type.
    Type{iNucleus} = Element_;
    
  % Check if line contains connection data.  
  elseif strncmp(line_,'CONECT',6)
    % Initialize connetion array.
    if isempty(Connected)
      Connected = cell(1,iNucleus);
    end
    if any(line_(7:end)=='*')
      continue
    end
    l = strtrim(line_);
    l = reshape(l(7:end),5,[]).';
    numbers = str2num(l).';
    if isempty(numbers), continue; end
    referenceNucleus = numbers(1);
    ConnectedNuclei = numbers(2:end);
    
    Connected{referenceNucleus} = unique([Connected{referenceNucleus},ConnectedNuclei]);
    
    for index_ = Connected{referenceNucleus}
      if index_ == referenceNucleus
        continue
      end
      if index_>nLines
        disp('CONECT data is not readable: check that whitspaces separate each number.');
        continue
      end
      if ~isempty(Connected{index_})
        Connected{index_} = [Connected{index_},referenceNucleus];
      else
        Connected{index_} = referenceNucleus;
      end
      Connected{index_} = unique(Connected{index_});
    end
    
  elseif strncmp(line_,'CRYST1',6)
    % Parse information about unit cell
    
    UnitCell.isUnitCell = true;
    values_ = sscanf(line_(7:54),'%f %f %f %f %f %f');
    UnitCell.ABC = values_(1:3)*System.angstrom; % angstrom -> m
    UnitCell.Angles = values_(4:6)*pi/180; % degree -> rad
    
  end
  
end

% Removed unused array slots.
isSolvent(iNucleus+1:end) = [];
Coordinates(iNucleus+1:end,:) = [];
pdbID(iNucleus+1:end) = [];
Type(iNucleus+1:end) = [];
Indices_nonWater = find(~isSolvent);
end
