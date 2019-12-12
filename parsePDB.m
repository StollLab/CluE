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
      continue;
    end
    
    % ConnectedNuclei = sscanf(line_(7:end),'%f %f %f %f %f');
    
    if length(line_) > 27
      ConnectedNuclei(5) = sscanf(line_(27:31),'%f');
    end
    if length(line_) > 22
      ConnectedNuclei(4) = sscanf(line_(22:26),'%f');
    end
    if length(line_) > 17
      ConnectedNuclei(3) = sscanf(line_(17:21),'%f');
    end
    if length(line_) > 12
      ConnectedNuclei(2) = sscanf(line_(12:16),'%f');
    end
    if length(line_) > 7
      ConnectedNuclei(1) = sscanf(line_(7:11),'%f');
    end
    
    if isempty(ConnectedNuclei)
      continue;
    end
    referenceNuclei = ConnectedNuclei(1);
    if referenceNuclei > iline % ignore connections where connections have no whitespaces.
      continue;
    end
    %     try
    if ~isempty(ConnectedNuclei)
      if ~isempty(Connected{referenceNuclei})
        Connected{referenceNuclei} = [Connected{referenceNuclei},ConnectedNuclei];
      else
        Connected{referenceNuclei} = ConnectedNuclei;
      end
      Connected{referenceNuclei} = unique(Connected{referenceNuclei});
      
      for index_ = Connected{referenceNuclei}
        if index_ == referenceNuclei
          continue;
        end
        if index_>nLines
          disp('CONECT data is not readable: check that whitspaces separate each number.');
          continue
        end
        if ~isempty(Connected{index_})
          Connected{index_} = [Connected{index_},referenceNuclei];
        else
          Connected{index_} = referenceNuclei;
        end
        Connected{index_} = unique(Connected{index_});
      end
      %     catch
      %       Connected{referenceNuclei} = ConnectedNuclei';
    end
  elseif strncmp(line_,'CRYST1',6)
    % Parse information about unit cell
    
    UnitCell.isUnitCell = true;
    UnitCell.ABC = sscanf(line_(7:33),'%f %f %f')*System.angstrom; % angstrom -> m
    UnitCell.Angles = sscanf(line_(34:54),'%f %f %f')*pi/180;
    
  end
  
end

end
