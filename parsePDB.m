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

function [Coordinates,Type,UnitCell,Connected, Indices_nonWater, pdbID,numberH] = parsePDB(filename,System)

fh = fopen(filename);
allLines = textscan(fh,'%s','whitespace','','delimiter','\n');
fclose(fh);
allLines = allLines{1};

nLines = numel(allLines);

Connected = {};
Indices_nonWater = [];
UnitCell.isUnitCell = false;
numberH = [0,0,0];
iNucleus = 0;
for iline = 1:nLines
  
  line_ = allLines{iline};
  
  if strncmp(line_,'ATOM',4) || strncmp(line_,'HETATM',6)
    % Parse information about atom
    iNucleus = iNucleus + 1;
    pdbID(iNucleus) = sscanf(line_(7:12),'%f');
    Connected{iNucleus} = [];
    Coordinates(iNucleus,:) = sscanf(line_(31:54),'%f %f %f')*System.angstrom; % angstrom -> m
    ResidueName_ = line_(18:20);
    Element_ = strtrim(line_(77:78));
    if ~strcmp(ResidueName_,'WAT')  && ~strcmp(ResidueName_,'SOL')
      if isempty(Indices_nonWater)
        Indices_nonWater = iNucleus;
      else
        Indices_nonWater = [Indices_nonWater,iNucleus];
      end
    end
    if strcmp(Element_,'H')
      if strcmp(ResidueName_,'WAT') || strcmp(ResidueName_,'SOL')
        if System.D2O, Element_ = 'D'; end
      else
        if System.deuterateProtein, Element_ = 'D'; end
        if System.solventOnly, Element_ = 'null'; end
      end
    end
    if strcmp(Element_,'H')
      numberH(1) = numberH(1) + 1;
      numberH(3) = numberH(3) + 1;
    elseif strcmp(Element_,'D')
      numberH(2) = numberH(2) + 1;
      numberH(3) = numberH(3) + 1;
    end
    Type{iNucleus} = Element_;
    
  elseif strncmp(line_,'CONECT',6)
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

end
