% Load pdb file into structure
%
%  Nuclei = parseNuclei(System,Method,datafile)
%
% Input:
%   System       structure with fields for the spin system
%   Method       structure with fields for the method
%   pdbFileName  name of PDB file

function [Nuclei, System]= parseNuclei(System,Method,Data,pdbFileName)

% set values to unspecified fields
System = setIsotopeDefaults(System,Method);
spinCenter = System.spinCenter;

Nuclei.nStates = System.nStates; 
Nuclei.number = 0;

nuT = System.Methyl.tunnel_splitting;

% Define spin operators.
%Nuclei = 
defineSpinOperators(); %Nuclei,System, Method);

% Copy graphCriterion to Nuclei.
Nuclei.graphCriterion = Method.graphCriterion;

% scale volume
scaleFactor = System.scale;

Nuclei.dataSource = pdbFileName;

% open data file
if ~strcmp(pdbFileName,'System.RandomEnsemble.include')
[pdbCoordinates,Type,UnitCell,Connected,Indices_nonSolvent,pdbID,MoleculeID,...
  numberH,isSolvent,isWater,Exchangeable,VanDerWaalsRadii] ...
  = parsePDB(pdbFileName,System);
else
  pdbCoordinates = [];
  isSolvent = [];
  Type = {};
  UnitCell.isUnitCell = false;
  numberH = [0,0,0];
end
Nuclei.isSolvent = isSolvent;

% Set electron coordinates in the pdb frame.
Nuclei.Electron_pdbCoordinates = getElectronCoordinates(System,pdbCoordinates,pdbID);

% number of PDB entries used
Npdb = length(Type);

if System.Methyl.include
  [Methyl_Data,pdbCoordinates,Type,UnitCell,Connected,Indices_nonSolvent] = ...
    findMethyls(System, pdbCoordinates,Type,UnitCell,Connected,...
    Indices_nonSolvent);
  Nuclei.Methyl_Data = Methyl_Data;
else
  Nuclei.Methyl_Data = [];
end

% Get nuclear quadrupole parameters.
if System.nuclear_quadrupole
  Nuclei.warnings.setQuadrupoleTensor = false;  
  % Water Quadrupole Values
  % Edmonds, D. T.; Mackay, A. L. 
  % The Pure Quadrupole Resonance of the Deuteron in Ice. 
  % Journal of Magnetic Resonance (1969) 1975, 20 (3), 515–519. 
  % https://doi.org/10.1016/0022-2364(75)90008-6.
  %   eta = 0.112*ones(1,Npdb*numberUnitCells);
  %   e2qQh = 213.4e3*ones(1,Npdb*numberUnitCells); % Hz
  
end


% Decide whether or not to add randomly distributed hard spheres.
if System.RandomEnsemble.include
  
  [pdbCoordinates,Type,Connected,Indices_nonSolvent,pdbID,MoleculeID,...
    numberH,isSolvent,isWater,Exchangeable,VanDerWaalsRadii] = ...
    addRandomSpins(System, Nuclei.Electron_pdbCoordinates, ...
    pdbCoordinates,Type,Connected,Indices_nonSolvent,pdbID,...
    MoleculeID,numberH,isSolvent,isWater,Exchangeable,VanDerWaalsRadii);
  
  % number of PDB entries used
  Npdb = length(Type);
  
  if length(pdbID) ~= npdb || length(Indices_nonSolvent) ~= npdb ...
      || length(VanDerWaalsRadii) ~= npdb
    error('Inconsistent output from addRandomSpins().')
  end
end



System.UnitCell = UnitCell;
cellshifts = getCellShifts(pdbCoordinates, Nuclei.Electron_pdbCoordinates,...
  System);
numberUnitCells = size(cellshifts,1);

if ~System.limitToSpinHalf
  Nuclei.quadrupole2lab = zeros(3,3,numberH(2)*numberUnitCells);
  Nuclei.Qtensor = zeros(3,3,Npdb*numberUnitCells);
  Nuclei.quadrupoleXaxis = zeros(numberH(2)*numberUnitCells,3);
  Nuclei.quadrupoleYaxis = zeros(numberH(2)*numberUnitCells,3);
  Nuclei.quadrupoleZaxis = zeros(numberH(2)*numberUnitCells,3);
end

num_ = Npdb*numberUnitCells;

% Nuclei.hyperfine2lab = sparse(3*3,num_);
Nuclei.Atensor = sparse(num_,9);
Nuclei.FermiContact = sparse(num_,1);
Nuclei.Azz = sparse(num_,1);

Nuclei.number_1H_exchangeable = 0;
Nuclei.number_1H_nonExchangeable = 0;
Nuclei.number_2H_exchangeable = 0;
Nuclei.number_2H_nonExchangeable = 0;

% loop over x unit cell spacings
iNuc = uint32(0);
for uc = 1:numberUnitCells
  
  % 3-vector offset to put each nucleus in the correct unit cell
  Delta_R = cellshifts(uc,:);
  
  ElectronCenteredCoordinates = ...
    scaleFactor*(pdbCoordinates + Delta_R - Nuclei.Electron_pdbCoordinates);
  
  % loop over all nuclei
  for inucleus = 1:size(Type,2)
    type = Type{inucleus};
    
    % get nuclear connection data
    try
      Conect = Connected{inucleus};
    catch
      Conect = {};
    end
    
    % set nuclear coordinates relative to the electron
    NuclearCoordinates = ElectronCenteredCoordinates(inucleus,:);
    
    % skip if the electron-nuclear separation is over the set cutoff
    if norm(NuclearCoordinates)>System.load_radius*scaleFactor
      continue
    end
    if norm(NuclearCoordinates)<System.inner_radius*scaleFactor
      continue
    end
    
    %[Nuclei,iNuc] = 
    setNucleus(); %Nuclei,type, inucleus,iNuc, Delta_R,...
     % System,Method,isSolvent,Exchangeable,Type,ElectronCenteredCoordinates, ...
     % NuclearCoordinates,pdbCoordinates,Conect,MoleculeID,isWater,spinCenter);
    
  end
  
end


% translate origin to electron
System.Electron.Coordinates = [0,0,0];

% get number of nuclei
try
  Nuclei.number = uint32(size(Nuclei.Index,2));
catch
  warning('No nuclear spins found.')
  return
end 

% Define methyl tunneling for APAYDIN CLOUGH methyl coupling.
% J . PHYS. c (PROC. PHYS. SOC.), 1968, SER. 2, VOL. 1. PRINTED IN GREAT BRITAIN
% Nuclear magnetic resonance line shapes of methyl groups
% undergoing tunnelling rotation

Nuclei.methylTunnelingSplitting = sparse( ...
  double(Nuclei.number),double(Nuclei.number));
for iNuc = 1:Nuclei.number
  if strcmp(Nuclei.Type{iNuc}, 'CH3')
    Nuclei.methylTunnelingSplitting(iNuc+1,iNuc+2) = nuT;
    Nuclei.methylTunnelingSplitting(iNuc+1,iNuc+3) = nuT;
    Nuclei.methylTunnelingSplitting(iNuc+2,iNuc+3) = nuT;
  end
end


% Rotate coordinates if requested by user via System.X/Y/Z
% [Nuclei,System] = 
setOrientation(); %Nuclei,System, pdbCoordinates);


% Clear excess entries.
% Nuclei = 
cleanUpNuclei(); %Nuclei,System,Method,Npdb);

% Get coupling statistics.
% Nuclei = 
computeNuclearInteractions(); %Nuclei,System, Method,scaleFactor);

if Data.writeSpinPDB
  try
    writeSpinPDB(Nuclei,ones(1,Nuclei.number),...
      [Data.OutputData, '_spinSystem.pdb']);
  end
end

% Clean
Nuclei.State = [];

if ~Method.getNuclearStatistics
  Nuclei.Statistics = [];
  Nuclei.DistanceMatrix = [];
end

if ~Method.getNuclearContributions
  Nuclei.PDBCoordinates = [];
  Nuclei.Element = [];
end

if System.newIsotopologuePerOrientation
  Nuclei.MoleculeIDunique = unique(Nuclei.MoleculeID);
else
  Nuclei.MoleculeID = [];
  Nuclei.Connected = [];
  Nuclei.isWater = [];
end    

checkNuclei(Nuclei);
% end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function ... % Nuclei = 
  cleanUpNuclei()%Nuclei,System,Method,Npdb)

if ~System.limitToSpinHalf
  Nuclei.quadrupole2lab(:,:,Nuclei.number+1:end) = [];
  Nuclei.Qtensor(:,:,Nuclei.number+1:end) = [];
  Nuclei.quadrupoleXaxis(Nuclei.number+1:end,:) = [];
  Nuclei.quadrupoleYaxis(Nuclei.number+1:end,:) = [];
  Nuclei.quadrupoleZaxis(Nuclei.number+1:end,:) = [];
end

Nuclei.isSolvent(Nuclei.number+1:end) = [];
% Nuclei.hyperfine2lab(:,Nuclei.number+1:end) = [];
Nuclei.Atensor(Nuclei.number+1:end,:) = [];
Nuclei.FermiContact(Nuclei.number+1:end) = [];
Nuclei.Azz(Nuclei.number+1:end) = [];



if System.doPruneNuclei
  if System.newIsotopologuePerOrientation && ~Method.reparseNuclei
    Nuclei = newHydronIsotopologue(Nuclei,System);
    System.newIsotopologuePerOrientation = false;
  end
  keep = Nuclei.Spin == 1/2 | ...
    vecnorm(Nuclei.Coordinates') <= Method.cutoff.radius_nonSpinHalf(1);
  
  
  
  oldIndex = Nuclei.Index;
  newIndex = oldIndex;
  cumsum_keep = cumsum(keep);
  newIndex(keep)  = cumsum_keep(keep);
  newIndex(~keep) = 0;
  newIndex(end:Npdb)=0;
  Nuclei.Index = 1:sum(keep);
  Nuclei.Type = Nuclei.Type(keep);
  Nuclei.Element = Nuclei.Element(keep);
  for iNuc = 1:Nuclei.number
    Nuclei.Connected{iNuc} = newIndex(Nuclei.Connected{iNuc});
  end
  Nuclei.Connected = Nuclei.Connected(keep);
  
  Nuclei.Spin = Nuclei.Spin(keep); % hbar
  Nuclei.StateMultiplicity = Nuclei.StateMultiplicity(keep);
  Nuclei.Nuclear_g = Nuclei.Nuclear_g(keep);
  Nuclei.Coordinates = Nuclei.Coordinates(keep,:);
  Nuclei.PDBCoordinates = Nuclei.PDBCoordinates(keep,:);
  Nuclei.MoleculeID = Nuclei.MoleculeID(keep);
  Nuclei.Exchangeable = Nuclei.Exchangeable(keep);
  Nuclei.NumberStates = Nuclei.NumberStates(keep);
  Nuclei.valid = Nuclei.valid(keep);
  Nuclei.isWater = Nuclei.isWater(keep);
  Nuclei.Abundance = Nuclei.Abundance(keep);
  
  Nuclei.isSolvent = Nuclei.isSolvent(keep);
  
  if ~System.limitToSpinHalf
    Nuclei.quadrupole2lab = Nuclei.quadrupole2lab(:,:,keep);
    Nuclei.Qtensor = Nuclei.Qtensor(:,:,keep);
    Nuclei.quadrupoleXaxis = Nuclei.quadrupoleXaxis(keep,:);
    Nuclei.quadrupoleYaxis = Nuclei.quadrupoleYaxis(keep,:);
    Nuclei.quadrupoleZaxis = Nuclei.quadrupoleZaxis(keep,:);
  end
%     Nuclei.hyperfine2lab = Nuclei.hyperfine2lab(:,keep);
    Nuclei.Atensor = Nuclei.Atensor(keep,:);
    Nuclei.FermiContact = Nuclei.FermiContact(keep);
    Nuclei.Azz = Nuclei.Azz(keep);

  % get number of nuclei
  try
    Nuclei.number = uint32(size(Nuclei.Index,2));
  catch
    warning('No nuclear spins remaining after pruning.')
    return
  end
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Rotate coordinates if requested by user via System.X/Y/Z
%-------------------------------------------------------------------------------
function ... %[Nuclei,System] = 
    setOrientation()%Nuclei,System, pdbCoordinates)
  
containsXYZ = [isfield(System,'X') isfield(System,'Y') isfield(System,'Z')];
if sum(containsXYZ)>=2
  
  % if nucleus index is given in X/Y/Z, calculate vector from electron to that
  % nucleus
  if containsXYZ(3) 
    if iscell(System.Z) && numel(System.Z)==2
      z1_ = System.Z{1};
      z2_ = System.Z{2};
      System.Z = pdbCoordinates(z2_,:) - pdbCoordinates(z1_,:);
    elseif length(System.Z)==1
      System.Z = pdbCoordinates(System.Z,:) - Electron_pdbCoordinates;
    end
  end
  if containsXYZ(2)
    if iscell(System.Y) && numel(System.Y)==2
      y1_ = System.Y{1};
      y2_ = System.Y{2};
      System.Y = pdbCoordinates(y2_,:) - pdbCoordinates(y1_,:);
    elseif length(System.Y)==1
      System.Y = pdbCoordinates(System.Y,:) - Electron_pdbCoordinates;
    end
  end
  if containsXYZ(1)
    if iscell(System.X) && numel(System.X)==2
      x1_ = System.X{1};
      x2_ = System.X{2};
      System.X = pdbCoordinates(x2_,:) - pdbCoordinates(x1_,:);
    elseif length(System.X)==1
      System.X = pdbCoordinates(System.X,:) - Electron_pdbCoordinates;
    end
  end
  
  % find the x and z unit vectors if one was not specified
  if ~containsXYZ(1)
    System.X = cross(System.Y,System.Z);
  elseif ~containsXYZ(3)
    System.Z = cross(System.X,System.Y);
  end
  
  % rotate system
  Rotation = alignCoordinates(System.X,System.Z);
  Nuclei.Coordinates = Nuclei.Coordinates*Rotation';
  

end

if System.randomOrientation
  
  % Generate random Euler angles.
  alpha_ = rand()*2*pi;
  beta_ = acos( 2*rand(1) - 1);
  gamma_ = rand()*2*pi;
  
  % Get rotation matrix.
  Rotation = rotateZYZ(alpha_,beta_,gamma_);
  
  % Rotate coordinates.
  Nuclei.Coordinates = Nuclei.Coordinates*Rotation';
  
  % Save rotation matrix.
  Nuclei.RandomRotationMatrix = Rotation;
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function ... % Nuclei = 
  defineSpinOperators()%Nuclei,System, Method)

nuT = System.Methyl.tunnel_splitting;
maxSize = 10;
maxClusterSize = [0,6,6];
if System.Methyl.methylMethylCoupling
  methylFactor = 2*System.Methyl.include + 1;
  maxClusterSize(2) = min(maxSize,methylFactor*Method.order);
else
  maxClusterSize(2)= min(maxSize,Method.order + 3);
end
Nuclei.maxClusterSize = maxClusterSize;

Nuclei.SpinOperators = cell(1,4);

multiplicity = 1;
Nuclei.SpinOperators{multiplicity} = 1;

for multiplicity = 2:3
  S = (multiplicity-1)/2;
  Nuclei.SpinOperators{multiplicity} = generateSpinOperators(S,maxClusterSize(multiplicity),Method.gpu);
%   rot = generateRotationMatrices(spinDim,numberMethyl)
end
if Method.allowHDcoupling % allowMixedSpins
  [Nuclei.SpinOperators_HD,Nuclei.SpinXiXjOperators_HD,...
    Nuclei.SpinXiXjOperators_DH] ...
    = assembleMixedSpinOperator(...
    Nuclei.SpinOperators{2}{1},Nuclei.SpinOperators{3}{1},true);
  
  Nuclei.SpinOperators_DH = assembleMixedSpinOperator(...
    Nuclei.SpinOperators{3}{1},Nuclei.SpinOperators{2}{1}, false);
end
  
if System.Methyl.method==0
  % Nuclei.rotationalMatrix{clusterSize,numberOfMethyls}
  % generateRotationMatrices(spinMultiplicity,numberOfMethyls)
  Nuclei.rotationalMatrix{1,1} = -nuT/3*generateRotationMatrices(2^3,1);
  Nuclei.rotationalMatrix{2,1} = -nuT/3*generateRotationMatrices(2*2^3,1);
  if System.Methyl.methylMethylCoupling
    Nuclei.rotationalMatrix{2,2} = -nuT/3*generateRotationMatrices(2^3*2^3,2);
  end
%   
%   Nuclei.rotationalMatrix{2,1} = -nuT/3*generateRotationMatrices(2,1);
%   Nuclei.rotationalMatrix{4,1} = -nuT/3*generateRotationMatrices(4,1);
%   Nuclei.rotationalMatrix{4,2} = -nuT/3*generateRotationMatrices(4,2);
%   Nuclei.rotationalMatrix{3,1} = -nuT/3*generateRotationMatrices(3,1);
%   Nuclei.rotationalMatrix{9,1} = -nuT/3*generateRotationMatrices(9,1);
%   Nuclei.rotationalMatrix{9,1} = -nuT/3*generateRotationMatrices(9,2);
else
  Nuclei.rotationalMatrix = [];
end
Nuclei.SpinXiXjOperators = generateXiXjSpinOperators(1,maxClusterSize(3));
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function ... %[Nuclei,iNuc] = 
  setNucleus()%Nuclei,type, inucleus,iNuc, Delta_R,...
%   System,Method,isSolvent,Exchangeable,Type,ElectronCenteredCoordinates, ...
%   NuclearCoordinates,pdbCoordinates,Conect,MoleculeID,isWater,spinCenter)
% H =============================================================
isProtium_ = (strcmp(type,'H') && System.protium);

if System.newIsotopologuePerOrientation && ~Method.reparseNuclei
  doParseAs1H_ = isProtium_ && ~isSolvent(inucleus);
else
  
  isDeuteriumTurnedProtium_ = ( strcmp(type,'D') && isSolvent(inucleus) ...
    && (rand() > System.deuteriumFraction) );
  
  doParseAs1H_ = ...
    (isProtium_ || isDeuteriumTurnedProtium_);
  
end

if  doParseAs1H_
  iNuc = iNuc +1;
  Nuclei.Index(iNuc) = iNuc;
  Nuclei.Type{iNuc} = '1H';
  Nuclei.Element{iNuc} = type;
  Nuclei.Connected{iNuc} = Conect;
  Nuclei.Spin(iNuc) = 0.5; % hbar
  Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
  Nuclei.Nuclear_g(iNuc) = 5.58569;
  Nuclei.Coordinates((iNuc),:) = NuclearCoordinates;
  Nuclei.PDBCoordinates((iNuc),:)= pdbCoordinates(inucleus,:) + Delta_R;
  Nuclei.MoleculeID(iNuc) = MoleculeID(inucleus);
  Nuclei.Exchangeable(iNuc) = Exchangeable(inucleus);
  Nuclei.NumberStates(iNuc) = int8(2);
  Nuclei.MethylID(iNuc) = 0;
  Nuclei.valid(iNuc)= true;
  Nuclei.Abundance(iNuc) = 1;
  Nuclei.isSolvent(iNuc) = isSolvent(inucleus);
  Nuclei.isWater(iNuc) = isWater(inucleus);
  if Nuclei.Exchangeable
    Nuclei.number_1H_exchangeable + Nuclei.number_1H_exchangeable + 1;
  else
    Nuclei.number_1H_nonExchangeable + Nuclei.number_1H_nonExchangeable + 1;
  end
  % CH3_A =========================================================
elseif strcmp(type,'CH3')
  iNuc = iNuc +1;
  Nuclei.Index(iNuc) = iNuc;
  Nuclei.Type{iNuc} = 'CH3';
  Nuclei.Element{iNuc} = type;
  Nuclei.Connected{iNuc} = Conect;
  Nuclei.Spin(iNuc) = 1/2; % hbar
  Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
  Nuclei.Nuclear_g(iNuc) = 5.58569;
  Nuclei.Coordinates((iNuc),:) = NuclearCoordinates;
  Nuclei.PDBCoordinates((iNuc),:)= pdbCoordinates(inucleus,:) + Delta_R;
  Nuclei.MoleculeID(iNuc) = MoleculeID(inucleus);
  Nuclei.Exchangeable(iNuc) = false;
  Nuclei.NumberStates(iNuc) = int8(8);
  Nuclei.valid(iNuc)= System.Methyl.method~=2; %
  
  ID_ref_ = find(Methyl_Data.ID(:,1)==inucleus);
  methylID_ = Methyl_Data.ID(ID_ref_,2) + max(Methyl_Data.ID(:,2))*(uc - 1);
  Nuclei.Group_ID{iNuc} = methylID_;
  Nuclei.MethylID(iNuc) = -methylID_;
  Nuclei.Auxiliary_ID(iNuc,:) = Methyl_Data.Hydron_ID{inucleus};
  
  
  [Nuclei.State{iNuc}, Nuclei.Abundance(iNuc)] = getMethylState(System);
  
  
  iNuc = iNuc +1;
  Nuclei.Index(iNuc) = iNuc;
  Nuclei.Type{iNuc} = 'CH3_1H';
  Nuclei.Element{iNuc} = type;
  Nuclei.Connected{iNuc} = Conect;
  Nuclei.Spin(iNuc) = 0.5; % hbar
  Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
  Nuclei.Nuclear_g(iNuc) = 5.58569;
  Nuclei.Coordinates((iNuc),:) = ...
    scaleFactor*(...
    Methyl_Data.Hydron_Coordinates{inucleus}(1,:) ...
    + Delta_R - Nuclei.Electron_pdbCoordinates);
  Nuclei.PDBCoordinates((iNuc),:) = ...
    Methyl_Data.Hydron_Coordinates{inucleus}(1,:) + Delta_R;
  Nuclei.NumberStates(iNuc) = int8(2);
  Nuclei.MethylID(iNuc) = methylID_;
  Nuclei.valid(iNuc) = System.Methyl.method==2;
  Nuclei.Abundance(iNuc) = 1;
  
  iNuc = iNuc +1;
  Nuclei.Index(iNuc) = iNuc;
  Nuclei.Type{iNuc} = 'CH3_1H';
  Nuclei.Element{iNuc} = type;
  Nuclei.Connected{iNuc} = Conect;
  Nuclei.Spin(iNuc) = 0.5; % hbar
  Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
  Nuclei.Nuclear_g(iNuc) = 5.58569;
  Nuclei.Coordinates((iNuc),:) = ...
    scaleFactor*(...
    Methyl_Data.Hydron_Coordinates{inucleus}(2,:) ...
    + Delta_R - Nuclei.Electron_pdbCoordinates);
  Nuclei.PDBCoordinates((iNuc),:) = ...
    Methyl_Data.Hydron_Coordinates{inucleus}(2,:) + Delta_R;
  Nuclei.NumberStates(iNuc) = int8(2);
  Nuclei.MethylID(iNuc) = methylID_;
  Nuclei.valid(iNuc)= System.Methyl.method==2;
  Nuclei.Abundance(iNuc) = 1;
  
  iNuc = iNuc +1;
  Nuclei.Index(iNuc) = iNuc;
  Nuclei.Type{iNuc} = 'CH3_1H';
  Nuclei.Element{iNuc} = type;
  Nuclei.Connected{iNuc} = Conect;
  Nuclei.Spin(iNuc) = 0.5; % hbar
  Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
  Nuclei.Nuclear_g(iNuc) = 5.58569;
  Nuclei.Coordinates((iNuc),:) = ...
    scaleFactor*(...
    Methyl_Data.Hydron_Coordinates{inucleus}(3,:) ...
    + Delta_R - Nuclei.Electron_pdbCoordinates);
  Nuclei.PDBCoordinates((iNuc),:) = ...
    Methyl_Data.Hydron_Coordinates{inucleus}(3,:) + Delta_R;
  Nuclei.NumberStates(iNuc) = int8(2);
  Nuclei.MethylID(iNuc) = methylID_;
  Nuclei.valid(iNuc)= System.Methyl.method==2;
  Nuclei.isWater(iNuc) = false;
  Nuclei.Abundance(iNuc) = 1;
  
  % D =============================================================
elseif strcmp(type,'D') && System.deuterium
  if Exchangeable(inucleus)
    Nuclei.number_2H_exchangeable + Nuclei.number_2H_exchangeable + 1;
  else
    Nuclei.number_2H_nonExchangeable + Nuclei.number_2H_nonExchangeable + 1;
  end
  if System.limitToSpinHalf
    return;
  end
  iNuc = iNuc +1;
  Nuclei.Index(iNuc) = iNuc;
  Nuclei.Type{iNuc} = '2H';
  Nuclei.Element{iNuc} = type;
  Nuclei.Connected{iNuc} = Conect;
  Nuclei.Spin(iNuc) = 1; % hbar
  Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
  Nuclei.Nuclear_g(iNuc) = 0.857438;
  Nuclei.Coordinates((iNuc),:) = NuclearCoordinates;
  Nuclei.PDBCoordinates((iNuc),:)= pdbCoordinates(inucleus,:) + Delta_R;
  %Nuclei.pdbID(iNuc) = pdbID(inucleus);
  Nuclei.MoleculeID(iNuc) = MoleculeID(inucleus);
  Nuclei.Exchangeable(iNuc) = Exchangeable(inucleus);
  Nuclei.NumberStates(iNuc) = int8(3);
  Nuclei.MethylID(iNuc) = 0;
  Nuclei.valid(iNuc)= true;
  Nuclei.isWater(iNuc) = isWater(inucleus);
  Nuclei.Abundance(iNuc) = 1;
  Nuclei.isSolvent(iNuc) = isSolvent(inucleus);
  
  
  % Set up quadrupole tensors for water deuterons
  if System.nuclear_quadrupole && ~System.spinHalfOnly
    
    
    if isempty(Conect)
      error(['Nucleus %d is not connected to anything - ',...
        'cannot build NQ tensor.'],inucleus);
    end
    for iconnect = Conect
      switch Type{iconnect}
        case {'O','C'}
          zQ = ElectronCenteredCoordinates(iconnect,:) - NuclearCoordinates;
        case {'M','D'}
          xQ = ElectronCenteredCoordinates(iconnect,:) - NuclearCoordinates;
      end
    end
    
    if isWater(inucleus)
      
      % Water Quadrupole Values
      % Edmonds, D. T.; Mackay, A. L.
      % The Pure Quadrupole Resonance of the Deuteron in Ice.
      % Journal of Magnetic Resonance (1969) 1975, 20 (3), 515–519.
      % https://doi.org/10.1016/0022-2364(75)90008-6.
      eta_ = 0.112;
      e2qQh_ = 213.4e3; % Hz
    else
      % ORCA
      eta_ = 0; % from eta_ = 0.0161;
      e2qQh_ = 0.1945e6; % Hz
      xQ = [0,0,0];
    end
    Nuclei = setQuadrupoleTensor(e2qQh_,eta_,zQ,xQ,iNuc,Nuclei,System);
    
  end
  
  % C ============================================================
elseif strcmp(type,'C') && System.carbon
  iNuc = iNuc +1;
  Nuclei.Index(iNuc) = iNuc;
  Nuclei.Type{iNuc} = '13C';
  Nuclei.Element{iNuc} = type;
  Nuclei.Connected{iNuc} = Conect;
  Nuclei.Spin(iNuc) = 0.5; % hbar
  Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
  Nuclei.Nuclear_g(iNuc) = 1.4048;
  Nuclei.Coordinates((iNuc),:) = NuclearCoordinates;
  Nuclei.PDBCoordinates((iNuc),:)= pdbCoordinates(inucleus,:);
  %Nuclei.pdbID(iNuc) = pdbID(inucleus);
  Nuclei.MoleculeID(iNuc) = MoleculeID(inucleus);
  Nuclei.Exchangeable(iNuc) = Exchangeable(inucleus);
  Nuclei.NumberStates(iNuc) = int8(2);
  Nuclei.MethylID(iNuc) = 0;
  Nuclei.valid(iNuc)= true;
  Nuclei.isWater(iNuc) = isWater(inucleus);
  Nuclei.isSolvent(iNuc) = isSolvent(inucleus);
  
  Nuclei.Abundance(iNuc) = 0.0107;
  
  
  
  % N =============================================================
elseif strcmp(type,'N') && System.nitrogen  && ~System.limitToSpinHalf
  iNuc = iNuc +1;
  Nuclei.Index(iNuc) = iNuc;
  Nuclei.Type{iNuc} = '14N';
  Nuclei.Element{iNuc} = type;
  Nuclei.Connected{iNuc} = Conect;
  Nuclei.Spin(iNuc) = 1; % hbar
  Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
  Nuclei.Nuclear_g(iNuc) = 0.403761;
  Nuclei.Coordinates((iNuc),:) = NuclearCoordinates;
  Nuclei.PDBCoordinates((iNuc),:)= pdbCoordinates(inucleus,:);
  %Nuclei.pdbID(iNuc) = pdbID(inucleus);
  Nuclei.MoleculeID(iNuc) = MoleculeID(inucleus);
  Nuclei.Exchangeable(iNuc) = Exchangeable(inucleus);
  Nuclei.NumberStates(iNuc) = int8(3);
  Nuclei.MethylID(iNuc) = 0;
  Nuclei.valid(iNuc)= true;
  Nuclei.isWater(iNuc) = isWater(inucleus);
  Nuclei.isSolvent(iNuc) = isSolvent(inucleus);
  Nuclei.Abundance(iNuc) = 0.99632;
  
  
  switch spinCenter
    case 'TEMPO'
      if norm(NuclearCoordinates) < System.angstrom
        % Aurich, H. G.; Hahn, K.; Stork, K.; Weiss, W. Aminyloxide
        % (nitroxide)—XXIV.
        %Tetrahedron 1977, 33 (9), 969–975.
        % https://doi.org/10.1016/0040-4020(77)80210-X.
        % Nuclei.FermiContact(iNuc) = 42.7*1e6; % Hz
        
        % Owenius, R.; Engström, M.; Lindgren, M.; Huber, M.
        % Influence of Solvent Polarity and Hydrogen Bonding on the EPR
        % Parameters of a Nitroxide Spin Label Studied by 9-GHz and
        % 95-GHz EPR Spectroscopy and DFT Calculations.
        %J. Phys. Chem. A 2001, 105 (49), 10967–10977.
        % https://doi.org/10.1021/jp0116914.
        Nuclei.FermiContact(iNuc) = 31.528e+06; %Hz
        Nuclei.Azz(iNuc) = 90.801e+06; % Hz
        
        
        if isempty(Conect)
          error(['Nucleus %d is not connected to anything',...
            '- cannot build NQ tensor.'],inucleus);
        end
        for iconnect = Conect
          % Marsh, D.
          % Bonding in Nitroxide Spin Labels from 14 N
          % Electric–Quadrupole Interactions.
          % J. Phys. Chem. A 2015, 119 (5), 919–921.
          % https://doi.org/10.1021/jp512764w.
          
          switch Type{iconnect}
            case 'O'
              xQ = ElectronCenteredCoordinates(iconnect,:) ...
                - NuclearCoordinates;
            case 'C'
              yQ = ElectronCenteredCoordinates(iconnect,:) ...
                - NuclearCoordinates;
          end
        end
        zQ = cross(xQ,yQ);
        Atensor_L = setHyperfineTensor(Nuclei.Azz(iNuc),...
          Nuclei.FermiContact(iNuc),zQ,xQ,iNuc,Nuclei);
        Nuclei.Atensor(iNuc,:) = Atensor_L(:)';

        % Jeong, J.; Briere, T.; Sahoo, N.; Das, T. P.;
        % Ohira, S.; Nishiyama, O.
        % Theory of Nuclear Quadrupole Interactions of 14 N, 17O,
        % and 35 CI Nuclei in p-Cl-Ph-CH-N=TEMPO.
        % Z. Naturforsch 2002.
        
        if System.nuclear_quadrupole
          % e2qQh_ = 4.807*1e6; % Hz
          % eta_ = 0.408;
          
          % de Oliveira, M.; Knitsch, R.; Sajid, M.; Stute, A.;
          % Elmer, L.-M.; Kehr, G.; Erker, G.; Magon, C. J.;
          % Jeschke, G.; Eckert, H.
          % Aminoxyl Radicals of B/P Frustrated Lewis Pairs:
          % Refinement of the Spin-Hamiltonian Parameters by Field- and
          % Temperature-Dependent Pulsed EPR Spectroscopy.
          % PLoS ONE 2016, 11 (6), e0157944.
          % https://doi.org/10.1371/journal.pone.0157944.
          
          e2qQh_ = 3.5*1e6; % Hz
          eta_ = 0.68;
          
          Nuclei = setQuadrupoleTensor(...
            e2qQh_,eta_,zQ,xQ,iNuc,Nuclei,System);
        end
        
      end
      
    otherwise
      Nuclei.FermiContact(iNuc) = 0;
  end
  
  % Si ============================================================
elseif strcmp(type,'Si') && System.silicon
  iNuc = iNuc +1;
  Nuclei.Index(iNuc) = iNuc;
  Nuclei.Type{iNuc} = '29Si';
  Nuclei.Element{iNuc} = type;
  Nuclei.Connected{iNuc} = Conect;
  Nuclei.Spin(iNuc) = 0.5; % hbar
  Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
  Nuclei.Nuclear_g(iNuc) = -1.11058;
  Nuclei.Coordinates((iNuc),:) = NuclearCoordinates;
  Nuclei.PDBCoordinates((iNuc),:)= pdbCoordinates(inucleus,:);
  %Nuclei.pdbID(iNuc) = pdbID(inucleus);
  Nuclei.MoleculeID(iNuc) = MoleculeID(inucleus);
  Nuclei.Exchangeable(iNuc) = Exchangeable(inucleus);
  Nuclei.NumberStates(iNuc) = int8(2);
  Nuclei.MethylID(iNuc) = 0;
  Nuclei.valid(iNuc)= true;
  Nuclei.isWater(iNuc) = isWater(inucleus);
  Nuclei.isSolvent(iNuc) = isSolvent(inucleus);
  
  Nuclei.Abundance(iNuc) = 0.046832;
  
  
  % electron ======================================================
elseif strcmp(type,'e')
  
  % electron, not a nucleus
  iNuc = iNuc +1;
  Nuclei.Index(iNuc) = iNuc;
  Nuclei.Type{iNuc} = 'e';
  Nuclei.Element{iNuc} = type;
  Nuclei.Connected{iNuc} = Conect;
  Nuclei.Spin(iNuc) = 0.5; % hbar
  Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
  Nuclei.Nuclear_g(iNuc) = 2.0023;
  Nuclei.Coordinates((iNuc),:) = NuclearCoordinates;
  Nuclei.PDBCoordinates((iNuc),:)= pdbCoordinates(inucleus,:);
  %Nuclei.pdbID(iNuc) = pdbID(inucleus);
  Nuclei.MoleculeID(iNuc) = MoleculeID(inucleus);
  Nuclei.Exchangeable(iNuc) = Exchangeable(inucleus);
  Nuclei.NumberStates(iNuc) = int8(2);
  Nuclei.MethylID(iNuc) = 0;
  Nuclei.valid(iNuc)= true;
  Nuclei.isWater(iNuc) = isWater(inucleus);
  Nuclei.isSolvent(iNuc) = isSolvent(inucleus);
  Nuclei.Abundance(iNuc) = 1;
  
elseif System.allAtoms
  iNuc = iNuc +1;
  Nuclei.Index(iNuc) = iNuc;
  Nuclei.Type{iNuc} = type;
  Nuclei.Element{iNuc} = type;
  Nuclei.Connected{iNuc} = Conect;
  Nuclei.Spin(iNuc) = 0; % hbar
  Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
  Nuclei.Nuclear_g(iNuc) = 0;
  Nuclei.Coordinates((iNuc),:) = NuclearCoordinates;
  Nuclei.PDBCoordinates((iNuc),:)= pdbCoordinates(inucleus,:);
  Nuclei.MoleculeID(iNuc) = MoleculeID(inucleus);
  Nuclei.Exchangeable(iNuc) = Exchangeable(inucleus);
  Nuclei.NumberStates(iNuc) = int8(2);
  Nuclei.valid(iNuc)= true;
  Nuclei.isWater(iNuc) = isWater(inucleus);
  Nuclei.isSolvent(iNuc) = isSolvent(inucleus);
  Nuclei.Abundance(iNuc) = 1;
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function ... %Nuclei = 
  computeNuclearInteractions()%Nuclei,System, Method,scaleFactor)
  
Nuclei.Statistics = getPairwiseStatistics(System, Nuclei);
Nuclei.DistanceMatrix = Nuclei.Statistics.DistanceMatrix;
if(Nuclei.Statistics.Distance > System.radius*scaleFactor)
  error(['Error in parseNuclei(): ','Nuclei beyond the distance cutoff ', ...
    'remain in the system.'])
end


if Method.lock_bAmax
  Nuclei.bAmax_lim = lock_bAmax(Nuclei.Statistics, Method);
  Method.cutoff.bAmax(:) = Nuclei.bAmax_lim;
  doAddCriterion = true;
  for icriterion = Method.Criteria
    if strcmp(icriterion,'bAmax')
      doAddCriterion = false;
    end
  end
  if doAddCriterion
    Method.Criteria{end+1} = 'bAmax';
  end
end
% Get the highest spin value. 
Nuclei.maxSpin = max(Nuclei.Spin);

Nuclei.Adjacency = getAdjacencyMatrix(System,Nuclei, Method);
Nuclei.AntiAdjacency = getAntiAdjacencyMatrix(System, Nuclei, Method); 
% Set the starting spin index and ending spin index.
Nuclei.startSpin = max(1, floor(Method.startSpin));
Nuclei.endSpin = min(Nuclei.number, floor(Method.endSpin));

% Check for consistancy.
if Nuclei.startSpin > Nuclei.endSpin
  disp(['Starting cluster spin cannot be greater than ending spin.  ', ...
    'Swapping assignment.']);
  Nuclei.startSpin = max(0, floor(Method.endSpin));
  Nuclei.endSpin = min(Nuclei.number, floor(Method.startSpin));
end

Nuclei.numberStartSpins = ...
  min(Nuclei.number, Nuclei.endSpin - Nuclei.startSpin + 1);


% set thermal energy
Nuclei.kT = System.kT;

% set thermal equilibrium state
[Nuclei.State, ~]= setThermalEnsembleState(System,Nuclei);

Nuclei.ZeemanStates = setRandomZeemanState(Nuclei);
[Nuclei.RandomDenityMatrices,Nuclei.RandomSpinVector] = ...
  setRandomDensityMatrix(Nuclei);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [Coordinates,Type,Connected,Indices_nonSolvent,pdbID,MoleculeID,...
  numberH,isSolvent,isWater,Exchangeable,VanDerWaalsRadii] ...
  = addRandomSpins(System, Electron_pdbCoordinates, ...
  Coordinates,Type,Connected,Indices_nonSolvent,pdbID,...
  MoleculeID,numberH,isSolvent,isWater,Exchangeable,VanDerWaalsRadii)

% radius to add spheres within
R_ = System.RandomEnsemble.radius;

% radius within which to avoid adding spheres
R0_ = System.RandomEnsemble.innerRadius;

% hard sphere radius
r0_ = System.RandomEnsemble.sphereRadius;

% Get random packing of hard spheres, centered on the origin.
[inCoor, N_]  = generateRandomConcentrationEnsemble(...
  System.RandomEnsemble.concentration,R_,R0_,2*r0_);

% Translate random ensemble to be centered on the elecron spin.
inCoor = inCoor + Electron_pdbCoordinates;

% Remove spheres that overlap with pdb input.
if ~isempty(Coordinates)
  [outCoor,  N_] = removeOverlap(inCoor,Coordinates,VanDerWaalsRadii);
else
  outCoor = inCoor;
end
% number of pdb coordinates
M_ = size(Coordinates,2);

% Add random spins to pdb spin.
Coordinates = [Coordinates; outCoor];

% Loop through the new spins.
for ii = M_+N_:-1: M_+1
  
  % Update spin info.
  Type{ii} = System.RandomEnsemble.Type;
  if strcmp(Type{ii},'H')
    numberH(1) = numberH(1) + 1;
  elseif strcmp(Type{ii},'D')
    numberH(2) = numberH(2) + 1;
  end
  Exchangeable(ii) = System.RandomEnsemble.Exchangeable;
  isSolvent(ii) = System.RandomEnsemble.isSolvent;
  isWater(ii) = System.RandomEnsemble.isWater;
  MoleculeID(ii) = ii;
  pdbID(ii) = - ii + M_;
  Indices_nonSolvent(ii) = ~isSolvent(ii);
  Connected{ii} = [];
  
end
numberH(3) = numberH(1) + numberH(3);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function R = alignCoordinates(MolecularX,MolecularZ)

NormZ = MolecularZ(:)/norm(MolecularZ);
NormX = MolecularX(:)/norm(MolecularX);
if abs(NormX'*NormZ) > 1e-12
  projectionXonZ = NormX'*NormZ;
  NormX = NormX - projectionXonZ*NormZ;
  NormX = NormX/norm(NormX);
end
Z = NormZ;
alpha = acos( Z(1)/sqrt( Z(1)^2 + Z(2)^2 ) );
if Z(2) < 0
  alpha = -alpha;
end
Rz1 = rotateZ(-alpha);
Z = Rz1*Z;

theta = acos(Z(3));
Ry2 = rotateY(-theta);

X = Ry2*Rz1*NormX;
phi = acos(X(1)/sqrt( X(1)^2 + X(2)^2 ) );
if X(2) < 0
  phi = -phi;
end

Rz3 = rotateZ(-phi);

R = Rz3*Ry2*Rz1;
Z = R*NormZ;
Z = Z/norm(Z);
X = R*NormX;
X = X/norm(X);
X = X - (X'*Z)*Z;
X = X/norm(X);
Z = abs(Z - [0;0;1]);
X = abs(X - [1;0;0]);
if (sum(X) + sum(Z)) > 1e-12
  disp(sum(X) + sum(Z));
  warning('Coordinates not aligned.');
  
  if (sum(X) + sum(Z)) > 1e-9
    disp(sum(X) + sum(Z));
    error('Coordinates not aligned.');
  end
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Sets the initial bath state with a Boltzmann distribution.
% Both output variables contain the same information, but are formated
% differently.
function [State,ZeemanStates] = setThermalEnsembleState(System,Nuclei)

ZeemanStates = zeros(Nuclei.number, 2*Nuclei.maxSpin+1);

% Loop through all nuclei.
for iinucleus = Nuclei.number:-1:1
  
  spin_multiplicity = Nuclei.NumberStates(iinucleus);
  
  if isfield(Nuclei,'State') && (length(Nuclei.State) >= iinucleus) ...
      && ~isempty(Nuclei.State{iinucleus})
    
    % Assign state.
    State{iinucleus} = Nuclei.State{iinucleus};
    
    ZeemanStates(iinucleus,1:length(Nuclei.State{iinucleus})) = ...
      Nuclei.State{iinucleus};
    continue
    
  end
  
  % Initialize state.
  State{iinucleus} = zeros(spin_multiplicity,1);
  
  
  % Set spin antiparallel to the magnetic field.
  mI = -double(Nuclei.Spin(iinucleus));
  
  % Set reference energy.
  en0 = -mI*System.magneticField*System.muN*Nuclei.Nuclear_g(iinucleus);
  
  % Loop through spin states.
  for ithresh = 1:spin_multiplicity
    % Increment z-projection.
    mI = double(ithresh - Nuclei.Spin(iinucleus) - 1);
    
    % Calculate the difference in energy from the reference energy.
    deltaE = -mI*System.magneticField*System.muN*Nuclei.Nuclear_g(iinucleus)...
      - en0;
    
    
    % Set population according to the Boltzmann distribution.
    State{iinucleus}(ithresh) = ...
      State{iinucleus}(ithresh) + exp(-deltaE/Nuclei.kT/2);
    
    ZeemanStates(iinucleus,ithresh)  = ...
      ZeemanStates(iinucleus,ithresh) + exp(-deltaE/Nuclei.kT/2);
  end
  
  % Normalize states.
  normalization = sqrt(State{iinucleus}'*State{iinucleus});
  State{iinucleus} = State{iinucleus}/normalization;
  ZeemanStates(iinucleus,:) = ZeemanStates(iinucleus,:)./normalization;
  
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


function [randDen,randSpin] = setRandomDensityMatrix(Nuclei)

N = Nuclei.number;
nStates = Nuclei.nStates; 
maxnStates = max(nStates);
maxMult = max(Nuclei.StateMultiplicity);
randDen = zeros(maxMult,maxMult, maxnStates,N);
randSpin = zeros(3, maxnStates,N);
S = cell(4,maxMult);
for iMult=2:maxMult
  spin_ = (iMult-1)/2;
  S{1,iMult} = spinX(spin_);
  S{2,iMult} = spinY(spin_);
  S{3,iMult} = spinZ(spin_);
end
% Loop through all nuclei.
for iSpin = 1:N
  mult_ = Nuclei.StateMultiplicity(iSpin);
  psi_ = rand(mult_,1,maxnStates)+rand(mult_,1,maxnStates)*1i;
  
  for iave =1:maxnStates
      rho_ = psi_(:,:,iave)*psi_(:,:,iave)';
      rho_ = rho_/trace(rho_);
      randDen(1:mult_,1:mult_,iave,iSpin) = rho_;
      
      for ix = 1:3
        randSpin(ix,iave,iSpin) = trace(rho_*S{ix,mult_});
      end

  end
  
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function System = setIsotopeDefaults(System, Method)

defaultValue = ~any(strcmp(Method.Criteria,'methyl only'));
if ~isfield(System,'protium')
  System.protium = defaultValue;
end
if ~isfield(System,'protiumFractionSolvent')
  System.protiumFractionSolvent = 1;
end
if ~isfield(System,'protiumFractionProtein')
  System.protiumFractionProtein = 1;
end
if ~isfield(System,'deuterium')
  System.deuterium = defaultValue;
end

if isfield(System,'hydrogen')
  System.protium = System.hydrogen;
  System.deuterium = System.hydrogen;
end
if ~isfield(System,'carbon')
  System.carbon = defaultValue;
end
if ~isfield(System,'nitrogen')
  System.nitrogen = defaultValue;
end
if ~isfield(System,'oxygen')
  System.oxygen = false;
end
if ~isfield(System,'silicon')
  System.silicon = defaultValue;
end
if ~isfield(System,'allAtoms')
  System.allAtoms = false;
end
if ~isfield(System,'scale')
  System.scale = 1;
end

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function  [Methyl_Data,Coordinates,Type_out,UnitCell,Connected,...
  Indices_nonWater] ...
  = findMethyls(System, Coordinates,Type,UnitCell,Connected,Indices_nonWater)

Type_out  = Type;
Number_Nuclei = length(Type);
Methyl_Data.number_methyls = 0;
Methyl_Data.Hydron_Coordinates = cell(1);
Methyl_Data.ID = [];
Methyl_Data.Hydron_radius  = cell(1);
Methyl_Data.Moment_of_Inertia = [];

for ispin = 1:Number_Nuclei
  
  % Check for a carbon.
  if ~strcmp(Type{ispin},'C')
    continue;
  end
  
  % Get connected nuclei.
  local_group = [Type{Connected{ispin}} ];
  
  % Skip lone carbons.
  if isempty(local_group)
    continue;
  end
  
  % Check if the carbon is a methyl carbon.
  indexH = [];
  for ii=0:length(local_group)
    
    % Remove the forth carbon bond.
    methyl = local_group([1:ii-1,(ii+1):end]);
    
    if strcmp(methyl,'CHHH')
      % Record indices.
      indexH = Connected{ispin}([ 2,3,4 ]);
      % Adjust for the forth carbon bond.
      indexH = indexH + (indexH <=ii);
      break;
    elseif strcmp(methyl,'HCHH')
      indexH = Connected{ispin}([ 1,3,4 ]);
      indexH = indexH + (indexH <=ii);
      break;
    elseif strcmp(methyl,'HHCH')
      indexH = Connected{ispin}([ 1,2,4]);
      indexH = indexH + (indexH <=ii);
      break;
    elseif strcmp(methyl,'HHHC')
      indexH = Connected{ispin}([ 1,2,3 ]);
      indexH = indexH + (indexH <=ii);
      break;
    end
    
  end
  
  % Check for no detected methyl groups.
  if isempty(indexH)
    continue;
  end
  
  new_index = Number_Nuclei + Methyl_Data.number_methyls + 1;
  %   Methyl_Data
  %   Coordinates
  %   Type,
  %   UnitCell
  Connected{new_index} =  Connected{ispin};
  %   Indices_nonWater
  Methyl_Data.number_methyls = Methyl_Data.number_methyls + 1;
  Methyl_Data.ID = [Methyl_Data.ID; new_index,Methyl_Data.number_methyls];
  
  Methyl_Data.Hydron_Coordinates{new_index} = zeros(3,3);
  Methyl_Data.Hydron_ID{new_index} = indexH;
  
  Center_of_Mass_Coor = sum(Coordinates(indexH,:),1)/3;
  Methyl_Data.Center_of_Mass_Coor{new_index} = Center_of_Mass_Coor;
  Coordinates(new_index,:) = Center_of_Mass_Coor;
  
  Methyl_Data.Hydron_radius{Methyl_Data.number_methyls} = [0,0,0];
  for iH = 1:3
    % Save hydron coordinates for Hamiltonian construction.
    Methyl_Data.Hydron_Coordinates{new_index}(iH,:) = Coordinates(indexH(iH),:);
    Methyl_Data.Hydron_radius{Methyl_Data.number_methyls}(iH) = norm(Coordinates(indexH(iH),:) - Center_of_Mass_Coor);
  end
  Methyl_Data.Moment_of_Inertia(Methyl_Data.number_methyls) = System.m1H*sum(Methyl_Data.Hydron_radius{Methyl_Data.number_methyls}.^2); 
  % Change hydrogens into methy pseudo-particles.
  Type_out{new_index} = 'CH3';
  Type_out{indexH(1)} = 'CH3_H_source';
  Type_out{indexH(2)} = 'CH3_H_source';
  Type_out{indexH(3)} = 'CH3_H_source';
  
  
end
ep = exp(2*1i*pi/3);
sq3 = sqrt(3);
%                       [aaa, aab, aba, abb, baa, bab, bba, bbb]
Methyl_Data.Transform = [0  , 1  , ep , 0  , ep', 0  , 0  , 0  ; ... % |+1/2 Ea>
  0  , 0  , 0  , ep', 0  , ep , 1  , 0  ; ... % |-1/2 Ea>
  0  , 1  , ep', 0  , ep , 0  , 0  , 0  ; ... % |+1/2 Eb>
  0  , 0  , 0  , ep , 0  , ep', 1  , 0  ; ... % |-1/2 Eb>
  sq3, 0  , 0  , 0  , 0  , 0  , 0  , 0  ; ... % |+3/2 A >
  0  , 1  , 1  , 0  , 1  , 0  , 0  , 0  ; ... % |+1/2 A >
  0  , 0  , 0  , 1  , 0  , 1  , 1  , 0  ; ... % |-1/2 A >
  0  , 0  , 0  , 0  , 0  , 0  , 0  , sq3]/sq3;% |-3/2 A >

Methyl_Data.Projection_A = zeros(4,8);
Methyl_Data.Projection_A(5,5) = 1;  Methyl_Data.Projection_A(6,6) = 1;
Methyl_Data.Projection_A(7,7) = 1;  Methyl_Data.Projection_A(8,8) = 1;

Methyl_Data.Projection_Ea = zeros(2,8);
Methyl_Data.Projection_Ea(1,1) = 1;  Methyl_Data.Projection_Ea(2,2) = 1;

Methyl_Data.Projection_Eb = zeros(2,8);
Methyl_Data.Projection_Eb(1,3) = 1;  Methyl_Data.Projection_Eb(2,4) = 1;


PA = zeros(8);
for ii=5:8
  for jj=5:8
    PA=PA +Methyl_Data.Transform(ii,:)'*Methyl_Data.Transform(jj,:);
  end
end
Pap = Methyl_Data.Transform(1,:)'*Methyl_Data.Transform(1,:);
Pam = Methyl_Data.Transform(2,:)'*Methyl_Data.Transform(2,:);

Pbp = Methyl_Data.Transform(3,:)'*Methyl_Data.Transform(3,:);
Pbm = Methyl_Data.Transform(4,:)'*Methyl_Data.Transform(4,:);

P=zeros(8,8,5);
P(:,:,1) = PA;
P(:,:,2) = Pap;
P(:,:,3) = Pam;
P(:,:,4) = Pbp;
P(:,:,5) = Pbm;
E = eye(2);
[Methyl_Data.Projection_3,Methyl_Data.Projection_4,...
  Methyl_Data.Projection_5,Methyl_Data.Projection_6] = ...
  getMethylProjections(P,E);
  
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function cellshifts = getCellShifts(pdbCoordinates, ElectronCoor, System)

if ~System.isUnitCell
  cellshifts = zeros(1,3);
  return;
end

Coordinates = pdbCoordinates - ElectronCoor;

pdbMinX = min(Coordinates,[],1);
pdbMaxX = max(Coordinates,[],1);
pdbDeltaX = pdbMaxX - pdbMinX;


% Find the minimum number of additional cells needed.
numUCplus = ceil( (System.radius - pdbMaxX)./pdbDeltaX );
numUCminus = ceil( ( System.radius - abs(pdbMinX) )./pdbDeltaX );


% consistency check
if any(numUCplus < 0) || any(numUCminus < 0)
  error('The electron spin is not located within the defined system.');
end
  
% Multiply the current number of cells by the number in this dimension.
numUnitCells = prod( 1 +  numUCplus + numUCminus);
if numUnitCells == 1
  cellshifts = zeros(1,3);
  return;
end

if  ~System.UnitCell.isUnitCell 
  error(['System radius extends beyond the volume spanned by the pdb, ', ...
    'and no unit cell information is available.  ',...
    'If this is intentional, please set System.isUnitCell = false; ',...
    'otherwise, please add unit cell information to the pdb.']);  
end

Angles = System.UnitCell.Angles;

ABC(1,:) = System.UnitCell.ABC(1)*[1,0,0];
ABC(2,:) = System.UnitCell.ABC(2)*[cos(Angles(3)),sin(Angles(3)),0];
cx = cos(Angles(2));
cy = ( cos(Angles(1)) - cos(Angles(2))*cos(Angles(3)) )/sin(Angles(3));
cz = sqrt(1-cx^2 - cy^2);
ABC(3,:) = System.UnitCell.ABC(3)*[cx,cy,cz];


idx = 0;
cellshifts = zeros(numUnitCells,3);

for a = -numUCminus(1):numUCplus(1)
  for b = -numUCminus(2):numUCplus(2)
    for c = -numUCminus(3):numUCplus(3)
      idx = idx+1;
      cellshifts(idx,:) = ABC(1,:)*a + ABC(2,:)*b + ABC(3,:)*c;
    end
  end
end

if idx ~= numUnitCells                                                         
  error('Error: created %d out of %d unit cells.',idx, numUnitCells)
end

R = vecnorm(cellshifts,2,2);
[~,idx] = sort(R);
cellshifts = cellshifts(idx,:);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function Electron_Coordinates = getElectronCoordinates(System,Coordinates,pdbID)

Electron_Coordinates = System.Electron.Coordinates;

if iscell(Electron_Coordinates)
  % find nuclei to average over
  replaceNuclei = [Electron_Coordinates{:}];
  ReplaceNuclei = zeros(size(replaceNuclei));
  for irep = 1:length(replaceNuclei)
    ReplaceNuclei(irep) = find(pdbID==replaceNuclei(irep));
  end
  % place the electron at the mean coordinates
  System.Electron.Coordinates = mean( Coordinates(ReplaceNuclei,:),1);
  
  % set initial electron coordinates to a 3-vector
  Electron_Coordinates = System.Electron.Coordinates;
  
  
  
elseif length(System.Electron.Coordinates) == 1
  
  replaceNucleus = System.Electron.Coordinates;
  System.Electron.Coordinates = Coordinates(replaceNucleus,:);
  Electron_Coordinates = Coordinates(replaceNucleus,:);
  
elseif length(System.Electron.Coordinates) == 2
  
  replaceNucleus1 = System.Electron.Coordinates(1);
  replaceNucleus2 = System.Electron.Coordinates(2);
  System.Electron.Coordinates = ...
    0.5*( Coordinates(replaceNucleus1,:) +Coordinates(replaceNucleus2,:));
  Electron_Coordinates = System.Electron.Coordinates;
  
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [State, gA]= getMethylState(System)

Nucleus_.kT=System.kT;
Nucleus_.Nuclear_g(1) = 5.58569;
Nucleus_.Spin(1) = 3/2;
Nucleus_.maxSpin = 3/2;
Nucleus_.NumberStates(1) = int8(4);
Nucleus_.number = 1;
State_A = setThermalEnsembleState(System,Nucleus_);

Nucleus_.Spin(1) = 1/2;
Nucleus_.NumberStates(1) = int8(2);
State_E = setThermalEnsembleState(System,Nucleus_);

Spin_Density = [State_E{1}; State_E{1}; State_A{1}].^2;
e_EkT= exp(-System.hbar*System.Methyl.tunnel_splitting/System.Methyl.kT);
Density = Spin_Density;%.*[e_EkT;e_EkT;e_EkT;e_EkT; 1;1;1;1];
Density = Density/sum(Density);
State = sqrt(Density);
gA = 1/(1+e_EkT);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function allPass = checkNuclei(Nuclei)

METHYLID = 1;
pass(METHYLID) = true;

for ii=1:Nuclei.number
  if Nuclei.MethylID(ii) > 0
  
    pass(METHYLID) = pass(METHYLID)*( ...
      sum(Nuclei.MethylID==Nuclei.MethylID(ii)) ==3);

  elseif Nuclei.MethylID(ii) < 0
    
    pass(METHYLID) = pass(METHYLID)*( ...
      sum(Nuclei.MethylID==Nuclei.MethylID(ii)) ==1);
  
  end
end
if ~pass(METHYLID)
  error(['Error in checkNuclei(): ', 'methylID is not correctly set up.']);
end

  allPass = all(pass(METHYLID)); 
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
