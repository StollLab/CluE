% Load pdb file into structure
%
%  Nuclei = parseNuclei(System,Method,datafile)
%
% Input:
%   System       structure with fields for the spin system
%   Method       structure with fields for the method
%   pdbFileName  name of PDB file

function [Nuclei, System]= centralSpinSystem(System,Method,Data,pdb_)
 
System = setCentralSpinSystemDefaults(System,Method);

ienum = 1;
PARTICLE_CENTRALSPIN = ienum; ienum = ienum + 1;
PARTICLE_ELECTRON = ienum; ienum = ienum + 1;

PARTICLE_HYDROGEN = ienum; ienum = ienum + 1;

PARTICLE_PROTIUM = ienum; ienum = ienum + 1;
PARTICLE_1H_EXCHANGEABLE = ienum; ienum = ienum + 1;
PARTICLE_1H_NONEXCHANGEABLE = ienum; ienum = ienum + 1;
PARTICLE_1H_METHYL = ienum; ienum = ienum + 1;

PARTICLE_DEUTERIUM = ienum; ienum = ienum + 1;
PARTICLE_2H_EXCHANGEABLE = ienum; ienum = ienum + 1;
PARTICLE_2H_NONEXCHANGEABLE = ienum; ienum = ienum + 1;
PARTICLE_2H_METHYL = ienum; ienum = ienum + 1;

PARTICLE_CARBON = ienum; ienum = ienum + 1;
PARTICLE_12C = ienum; ienum = ienum + 1;
PARTICLE_13C = ienum; ienum = ienum + 1;
PARTICLE_C_METHYL = ienum; ienum = ienum + 1;
PARTICLE_13C_METHYL = ienum; ienum = ienum + 1;

PARTICLE_NITROGEN = ienum; ienum = ienum + 1;
PARTICLE_14N = ienum; ienum = ienum + 1;
PARTICLE_15N = ienum; ienum = ienum + 1;

PARTICLE_OXYGEN = ienum; ienum = ienum + 1;
PARTICLE_16O = ienum; ienum = ienum + 1;
PARTICLE_17O = ienum; ienum = ienum + 1;
PARTICLE_SILICON = ienum; ienum = ienum + 1;
PARTICLE_28SI = ienum; ienum = ienum + 1;
PARTICLE_29SI = ienum; ienum = ienum + 1;

PARTICLE_UNSET = ienum; ienum = ienum + 1;
PARTICLE_NONE = ienum; ienum = ienum + 1;
PARTICLE_VOID = ienum; ienum = ienum + 1;
PARTICLE_ALL = ienum; ienum = ienum + 1;

PARTICLE_ENUM = ienum; ienum = ienum + 1;
if(PARTICLE_ENUM+1-ienum ~= 0)
  error(['Error in centralSpinSystem(): ', 'enums are incorrect.']);
end
associatedParticleMatrix_ = sparse(0,0); % indexed by CluE indices.
associatedPDBMatrix_ = sparse(pdb_.number,pdb_.number); % indexed by pdb ID.
cellShifts_              = [];
% methylIDs_               = [];
moleculeID_              = [];
number_                  = 0; 
numberMethyls_           = 0;
numberParticleClasses_   = 0;
particleClassID_         = [];
particles_               = cell(0,0);
%pdb_                     = parsePDBfile(Data.InputData, System.angstrom);
originVec_               = getElectronCoordinates(...
                            System,[pdb_.x,pdb_.y,pdb_.z],pdb_.serial);
pdbID_                   = []; 
residueList_             = cell(0,0); 
uniqueResidueList_       = {}; 
zeroIndex_               = {};


buildSystem();

% Finish setup.
coordinates_ = [];
refreshSystem();

% Print total number of particles.
fprintf('\nFound %d particles.\n', number_);

% Loop through all initialized ParticleClasses and print details.
for iitype = 1:numberParticleClasses_
  printParticleInfo( particles_{iitype} );
end


Nuclei = struct();
assembleNuclei();

return;

% Start of nested functions

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function uniqueID = addParticle(particleEnum, resName, ...
    associatedParticles, iparticle, ucpdb)
 
% Initialize return value.
uniqueID  = -1;
  
% Get PDB number.     
ipdb = mod(ucpdb, pdb_.number);
if ipdb==0, ipdb = pdb_.number; end
uc = ceil(ucpdb/pdb_.number);


% Check if particle type with resName has been initialized.
isUnique = true;

for itype=1:numberParticleClasses_
  if (particleEnum==particles_{itype}.particleEnum) && ...
    strcmp( resName, particles_{itype}.resName)
    isUnique = false;
  end
  
  if ~isUnique

    uniqueID = particles_{itype}.ID;
   
    % Record which particle class this particle is in;
    particleClassID_(iparticle) = itype;
     
     break;
  end
end

% Initialize if it does not already exist.
if isUnique && isParticleValid(particleEnum, pdb_.resName{ipdb}) 

  addParticle_new(particleEnum, resName);
  
  % Create new particle type. 
  uniqueID = particles_{numberParticleClasses_}.ID;

  % Record which particle class this particle is in;
  particleClassID_(iparticle) = particles_{numberParticleClasses_}.ID;

end

% Check if the particle is already a member of this particle class.
if any(pdbID_(particles_{uniqueID}.members)==ucpdb)
  return;
end

% Add parttcle to member set.
particles_{uniqueID}.members  = unique(...
  [particles_{uniqueID}.members,iparticle]);

particles_{uniqueID}.number = length(particles_{uniqueID}.members);

% TO DO: REPLACE associatedParticlesCollection WITH associatedPDBMatrix.
particles_{uniqueID}.associatedParticlesCollection{...
  particles_{uniqueID}.number} = associatedParticles;

if uc==1
  for iAsso = 1:numel(associatedParticles)
    iaParticle = associatedParticles(iAsso);
    associatedPDBMatrix_(iaParticle,ipdb) = iaParticle;
  end
end

% Record molecule ID number.
moleculeID_(iparticle) = pdb_.resSeq(ipdb) + (uc-1)*pdb_.number;

% Record PDB number.
if numel(pdbID_) >= iparticle && pdbID_(iparticle) ~= ucpdb
  error(['Error in addParticle(): ',...
    'particle ', num2str(iparticle), 'is already asigned pdb ID ',...
    num2str(pdbID_(iparticle)), ' and cannot be re-asigned pdb ID ', ...
    num2str(ucpdb), '.']);
end
  
pdbID_(iparticle) = ucpdb;


if iparticle == number_+1
  % Count this particle.
  number_ = number_ + 1;
elseif iparticle > number_ + 1
  error('Error in addParticle(): particle count is incorrect.');
end

% Update residue list.
residueList_{iparticle} = resName;

if ~particles_{uniqueID}.exchangeable && ~any(moleculeID_(iparticle)...
    ==moleculeID_(particles_{uniqueID}.membersOnePerMolecule))
  
  particles_{uniqueID}.membersOnePerMolecule = ...
    [particles_{uniqueID}.membersOnePerMolecule,iparticle];
  
  particles_{uniqueID}.numberMolecules = particles_{uniqueID}.numberMolecules+1;

end

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function addParticle_new(particleEnum, resName)     
 
% Update particle type number.  
numberParticleClasses_ = numberParticleClasses_ + 1;

% Set identification values.  
particles_{numberParticleClasses_}.particleName = ...
  getParticleString(particleEnum);
particles_{numberParticleClasses_}.resName = resName;
particles_{numberParticleClasses_}.particleEnum = particleEnum;
particles_{numberParticleClasses_}.ID = numberParticleClasses_;

% Initialize other values.
particles_{numberParticleClasses_}.members = [];
particles_{numberParticleClasses_}.membersOnePerMolecule = [];
particles_{numberParticleClasses_}.number = 0;
particles_{numberParticleClasses_}.numberMolecules = 0;

% Initialize option values.
particles_{numberParticleClasses_}.abundance = 1;
particles_{numberParticleClasses_}.active = true;
particles_{numberParticleClasses_}.atomType = PARTICLE_UNSET;
particles_{numberParticleClasses_}.associatedParticlesCollection = {};
particles_{numberParticleClasses_}.barrierPotential = 0;
particles_{numberParticleClasses_}.doRandomIsotopes = [];
particles_{numberParticleClasses_}.exchangeable = true;
particles_{numberParticleClasses_}.gFactor = 0;
particles_{numberParticleClasses_}.hf_Azz = 0;
particles_{numberParticleClasses_}.hf_FermiContact = 0;
particles_{numberParticleClasses_}.hf_x = {};
particles_{numberParticleClasses_}.hf_y = {};
particles_{numberParticleClasses_}.hf_z = {};
particles_{numberParticleClasses_}.isNucleus = false; 
particles_{numberParticleClasses_}.NQ_e2qQh = 0;
particles_{numberParticleClasses_}.NQ_eta = 0;
particles_{numberParticleClasses_}.NQ_x = {};
particles_{numberParticleClasses_}.NQ_y = {};
particles_{numberParticleClasses_}.NQ_z = {};
particles_{numberParticleClasses_}.spinMultiplicity = 1;
particles_{numberParticleClasses_}.switchParticle = PARTICLE_UNSET;
particles_{numberParticleClasses_}.extraCellSwitchParticle = PARTICLE_UNSET;
particles_{numberParticleClasses_}.tunnelSplitting = 0;

% Set default options.
setParticle(numberParticleClasses_, particleEnum,resName);

% Set user specified options.
setParticleOptions(numberParticleClasses_, System.particleOptions);

particleClassIndex = findParticleClasses(...
  particles_{numberParticleClasses_}.switchParticle, resName);
 
if numel(particleClassIndex)==1 ...
    && particles_{particleClassIndex}.doRandomIsotopes
  particles_{numberParticleClasses_}.abundance = ...
    1 - particles_{particleClassIndex}.abundance;
elseif numel(particleClassIndex) > 1
  error(['Error in addParticle_new():',...
    ' multiple switch particles found.']);
end

 
if isempty(particles_{numberParticleClasses_}.doRandomIsotopes)
  
  particles_{numberParticleClasses_}.doRandomIsotopes = ...
    (particles_{numberParticleClasses_}.abundance<1) && ...
    particles_{numberParticleClasses_}.switchParticle ~= PARTICLE_UNSET;
  
end

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function preliminaryParticle = analyzeParticle(index_pdb)
  
pdbParticle =   getParticleClass(pdb_.element{index_pdb});
preliminaryParticle.particleEnum = PARTICLE_NONE;

preliminaryParticle.associatedParticles = [];

switch pdbParticle
  % Hydrogen
  case {PARTICLE_HYDROGEN, PARTICLE_PROTIUM, PARTICLE_DEUTERIUM}
    
    if ~System.include_hydrogen
      return;
    end
    
    % Assume protium unless specified otherwise.
    isH1 = pdbParticle ~= PARTICLE_DEUTERIUM;
    
    if (isH1 && ~System.include_1H)  || (~isH1 && ~System.include_2H)
      return;
    end
    
    isExchangeable = isHydronExchangeable(index_pdb);
    
    % Set hydron type.
    if isH1 && isExchangeable
      preliminaryParticle.particleEnum = PARTICLE_1H_EXCHANGEABLE;
      
    elseif ~isH1 && isExchangeable
      preliminaryParticle.particleEnum = PARTICLE_2H_EXCHANGEABLE;
      
    elseif isH1 && ~isExchangeable
      preliminaryParticle.particleEnum = PARTICLE_1H_NONEXCHANGEABLE;
      
    else
      preliminaryParticle.particleEnum = PARTICLE_2H_NONEXCHANGEABLE;
    end
    
    % Carbon
  case PARTICLE_CARBON
    if System.Methyl.include
      
      % Check if carbon is a methyl carbon.
      methylHydrons =  isMethyl(index_pdb);
      
      % All indices will be unique if it is a methyl.
      if ~(methylHydrons(1)==methylHydrons(2) ...
          || methylHydrons(1) == methylHydrons(3) ...
          || methylHydrons(2) == methylHydrons(3))
        
        preliminaryParticle.particleEnum = PARTICLE_C_METHYL;
        preliminaryParticle.associatedParticles = methylHydrons;
        
        numberMethyls_ = numberMethyls_ + 1;
        
%         methylIDs_(numberMethyls_) = index_pdb;
      end
    elseif System.include_13C
      preliminaryParticle.particleEnum = PARTICLE_12C;
    end
    
    % Nitrogen
  case PARTICLE_NITROGEN
    if System.include_14N
      preliminaryParticle.particleEnum = PARTICLE_14N;
    elseif System.include_15N
      preliminaryParticle.particleEnum = PARTICLE_15N;
    end
    
    
    % Silicon
  case PARTICLE_SILICON
    if System.include_29Si
      preliminaryParticle.particleEnum = PARTICLE_29SI;
    end
    
    % Electron
  case PARTICLE_ELECTRON
    if System.include_electron
      preliminaryParticle.particleEnum = PARTICLE_ELECTRON;
    end
    
  otherwise
    preliminaryParticle.particleEnum = PARTICLE_NONE;
end

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% TO DO: THIS FUNCTION IS TEMPORARY.  CluE SHOULD BE RECONFIGURED TO WORK WITH
% THE DATA STRUCTURE OF centralSpinSystem(), WITHOUT A CONVERSION STEP.
function assembleNuclei()
  

% Initialize
Nuclei.Index = 1:(number_-1); % set here.
Nuclei.Type = cell(1,number_-1);
Nuclei.Element = cell(number_-1,1);
% Nuclei.Connected = cell(number_,1); % not needed.
Nuclei.Spin = zeros(1,number_-1);
Nuclei.StateMultiplicity = zeros(1,number_-1);
Nuclei.Nuclear_g = zeros(1,number_-1);
Nuclei.Coordinates = coordinates_(:,2:number_)'; % set here.
Nuclei.PDBCoordinates = zeros(number_ - 1,3);
Nuclei.pdbID = zeros(1,number_-1);
Nuclei.ucpdbID = zeros(1,number_-1);
Nuclei.MoleculeID = moleculeID_; % set here.
Nuclei.methylTunnelingSplitting = sparse(number_-1,number_-1);
Nuclei.Exchangeable = zeros(1,number_-1);
Nuclei.NumberStates = int8(zeros(1,number_-1));
Nuclei.MethylID = zeros(1,number_-1);
methylCounter = 0;
% Nuclei.valid = true(1,number_-1); % not needed.
% Nuclei.isWater = true(1,number_-1); % not needed.
Nuclei.isSolvent = true(1,number_-1);
Nuclei.valid = true(1,number_-1);

% Nuclei.hyperfine2lab = sparse(3*3,num_);
Nuclei.Atensor = sparse(number_-1,9);
Nuclei.FermiContact = sparse(number_-1,1);
Nuclei.Azz = sparse(number_-1,1);


Nuclei.number_1H_exchangeable = ...
  findParticles(PARTICLE_1H_EXCHANGEABLE, 'ALL', true);

Nuclei.number_1H_nonExchangeable = ...
  findParticles(PARTICLE_1H_NONEXCHANGEABLE, 'ALL', true)...
  + findParticles(PARTICLE_1H_METHYL, 'ALL', true);

Nuclei.number_2H_exchangeable = ...
  findParticles(PARTICLE_2H_EXCHANGEABLE, 'ALL', true);

Nuclei.number_2H_nonExchangeable = ...
  findParticles(PARTICLE_2H_NONEXCHANGEABLE, 'ALL', true)...
  + findParticles(PARTICLE_2H_METHYL, 'ALL', true);

if ~System.limitToSpinHalf
  Nuclei.Qtensor = zeros(3,3,number_-1);
end

for iparticle = 2:number_
  iNuc = iparticle - 1;
  
  itype = particleClassID_(iparticle);
  particleEnum = particles_{itype}.particleEnum;
  
  ipdb = mod(pdbID_(iparticle),pdb_.number);
  if ipdb==0, ipdb = pdb_.number; end
  
  uc = ceil(pdbID_(iparticle)/pdb_.number);
  
  Nuclei.Type{iNuc} = getParticleTypeString(particleEnum);
  Nuclei.Element{iNuc} = pdb_.element{ipdb};
  Nuclei.PDBCoordinates(iNuc,:) = [pdb_.x(ipdb),pdb_.y(ipdb),pdb_.z(ipdb)]...
    + cellShifts_(:,uc)';
  Nuclei.pdbID(iNuc) = pdbID_(iparticle);
  Nuclei.Spin(iNuc) = (particles_{itype}.spinMultiplicity - 1)/2;
  Nuclei.StateMultiplicity(iNuc) = particles_{itype}.spinMultiplicity;
  Nuclei.Nuclear_g(iNuc) = particles_{itype}.gFactor;
  Nuclei.Exchangeable(iNuc) = particles_{itype}.exchangeable;
  Nuclei.NumberStates(iNuc) = int8(particles_{itype}.spinMultiplicity);
  Nuclei.isSolvent(iNuc) = isMoleculeSolvent(particles_{itype}.resName);
  Nuclei.valid(iNuc) = particles_{itype}.active;
  
  if particles_{itype}.NQ_e2qQh ~= 0
    idx = particles_{itype}.members==iparticle;
    Nuclei.Qtensor(:,:,iNuc) = particles_{itype}.Qtensor{idx};
  end
  
  if particles_{itype}.hf_Azz ~= 0 || particles_{itype}.hf_FermiContact ~= 0
    idx = particles_{itype}.members==iparticle;
    Nuclei.Atensor(iNuc,:) = particles_{itype}.Atensor{idx}(:);
  end
  
  if particleEnum == PARTICLE_C_METHYL
    methylCounter = methylCounter + 1;

    %     hydrons = getAssociatedParticles(itype,iparticle);
    hydrons = nonzeros(associatedParticleMatrix_(:,iparticle)) - 1;
    
    nuT = particles_{itype}.tunnelSplitting;
    
    for ih = 1:2
      iH = hydrons(ih);
      for jh = ih+1:3 
        jH = hydrons(jh);
        Nuclei.methylTunnelingSplitting(iH,jH) = nuT;
      end
    end

    Nuclei.MethylID(iNuc) = -methylCounter;
    for ih = 1:3
      iH = hydrons(ih);
      Nuclei.MethylID(iH) = methylCounter;
    end
%     error(['Error in assembleNuclei(): ', ...
%       'Nuclei.MethylID not set.  Please update source code.'])
  end
  
end
  
defineSpinOperators();

Nuclei.nStates = System.nStates; 
Nuclei.number = 0;

% Define spin operators.
%Nuclei = 
 %Nuclei,System, Method);

% Copy graphCriterion to Nuclei.
Nuclei.graphCriterion = Method.graphCriterion;


% Set electron coordinates in the pdb frame.
Nuclei.Electron_pdbCoordinates = originVec_';
  
if System.Methyl.include %&& System.Methyl.method ~= 2
  Nuclei.Methyl_Data.number_methyls = methylCounter;
%     error('Error in centralSpinSystem(): methyl code not yet implemented.');
else
  Nuclei.Methyl_Data = [];
end

% Decide whether or not to add randomly distributed hard spheres.
if System.RandomEnsemble.include
  error('Error in centralSpinSystem(): addRandomSpins() not yet implemented.');
end  


% get number of nuclei
Nuclei.number = number_-1; 


% Rotate coordinates if requested by user via System.X/Y/Z
% [Nuclei,System] = 
setOrientation(); %Nuclei,System, pdbCoordinates);

% Clear excess entries.
% Nuclei = 
% cleanUpNuclei(); %Nuclei,System,Method,Npdb);

% Get coupling statistics.
% Nuclei = 
computeNuclearInteractions(); %Nuclei,System, Method,scaleFactor);

if Data.writeSpinPDB
  try
    writeSpinPDB(Nuclei,ones(1,Nuclei.number),...
      [Data.OutputData, '_spinSystem.pdb']);
  catch
    warning('Warning: could not write spin pdb.')
  end
end

% Clean
% Nuclei.State = [];

if ~Method.getNuclearStatistics
  Statistics = struct();
  Statistics.Nuclear_Dipole = Nuclei.Statistics.Nuclear_Dipole;
  Statistics.Nuclear_Dipole_x_iy_Z = Nuclei.Statistics.Nuclear_Dipole_x_iy_Z; 

  Nuclei.Statistics = Statistics;
  Nuclei.DistanceMatrix = [];
end

if ~Method.getNuclearContributions
  Nuclei.PDBCoordinates = [];
  Nuclei.Element = [];
end

% if System.newIsotopologuePerOrientation
  Nuclei.MoleculeIDunique = unique(Nuclei.MoleculeID);
% else
%   Nuclei.MoleculeID = [];
%   Nuclei.Connected = [];
%   Nuclei.isWater = [];
% end    

checkNuclei();  
return;
end 
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function buildAssociatedParticlematrix()

% Resize matrix.
associatedParticleMatrix_(number_,number_) = 0;

% Loop through all particles in system.
for iParticle = 1:number_

  % Get pdb information.
  ucpdb = pdbID_(iParticle);
  
  uc = ceil(ucpdb/pdb_.number);
  
  ipdb = mod(ucpdb,pdb_.number);
  if ipdb==0, ipdb = pdb_.number; end

  % Get non-zero column elements, the particles associated with ipdb. 
  associatedPDBcol = nonzeros(associatedPDBMatrix_(:,ipdb));
  for iAsso = 1:numel(associatedPDBcol)

    % Find associated particles if they exist.
    assoiParticle = find(pdbID_ ==  ...
      (associatedPDBcol(iAsso) + (uc-1)*pdb_.number));
    
    % Add missing associated particles.
    if isempty(assoiParticle)
      itype = particleClassID_(iParticle);

      assoiParticle = number_;

      associatedPart = analyzeParticle(associatedPDBcol(iAsso));

      addParticle(associatedPart.particleEnum, particles_{itype}.resName, ...
        [], number_+1, (uc-1)*pdb_.number + associatedPDBcol(iAsso));

      if assoiParticle==number_
        error(['Error in buildAssociatedParticlematrix(): ', ...
        'particle could not be added.']);
      else
        assoiParticle = number_;
      end

    elseif numel(assoiParticle)>1
      error(['Error in buildAssociatedParticlematrix(): ', ...
        'particle has been added multiple times.']);
    end

    associatedParticleMatrix_(assoiParticle,iParticle) = assoiParticle;
  end
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function buildExtraCells()

cellShifts_ = getCellShifts()';

generateUniqueResidueList();

% Print matrix.
disp("cells shifts (meter) = ");
disp(cellShifts_);

% Copy pdbID before erasing.
pdbID = pdbID_;

% Remove all particles but leave settings.
clearBath();

numUnitCells = size( cellShifts_,2);
if numUnitCells == 0
  return;
end

for itype = 1:numberParticleClasses_
  
  % Get particle label.
  particleEnum = particles_{itype}.particleEnum;
  
  if particleEnum==PARTICLE_CENTRALSPIN
    continue;
  end
  
  % Get number of particles in one cell.
  if particles_{itype}.exchangeable
    % Copy particle indices.
    indices = particles_{itype}.members;
    N = particles_{itype}.number;
  else
    % Copy particle indices.
    indices = particles_{itype}.membersOnePerMolecule;
    N = particles_{itype}.numberMolecules;
  end
  
  ucN =  numUnitCells*N ;

  clearMembers(itype);

  useRandomParticles = false;

  if particles_{itype}.switchParticle == PARTICLE_VOID ...
      && particles_{itype}.abundance < 1 

    useRandomParticles = true;

    % Get total number of particles in all unit cells.
    ucNtotal = numUnitCells*N;

    % Get number to use.
    pd = makedist('Binomial','N',ucNtotal,'p',particles_{itype}.abundance);
    ucN = random(pd);

    % Draw a random set of indices. 
    selectedIndices = randperm(ucNtotal,ucN);
  elseif particles_{itype}.extraCellSwitchParticle == PARTICLE_VOID ...
      && particles_{itype}.abundance < 1 

    useRandomParticles = true;

    % Get total number of particles in all unit cells outside the primary cell.
    ucNtotal = (numUnitCells-1)*N;

    % Get number to use.
    pd = makedist('Binomial','N',ucNtotal,'p',particles_{itype}.abundance);
    ucN = random(pd);

    % Use all indices from the primary cell 
    % and draw a random set of indices for the extra cells. 
    selectedIndices =[1:N, N + randperm(ucNtotal,ucN)];
    ucN = ucN + N;

  end

  % Loop through particles.
  for ii=1:ucN
    
    particleIndex =  ii;
    
    if useRandomParticles
      particleIndex = selectedIndices(ii); 
    end
    
    % Get unit cell.  Note that N = particles_{itype}.number.  The devision by
    % N rather than pdb_.number is because particleIndex runs over particles
    % of particles_{itype}.
    uc = ceil( particleIndex/ N );
    
    % Define index to store unit cell and ipdb info together.
    ucIndexOffset = (uc - 1)*pdb_.number;
        
    % Get index of particle copy in the primary cell.
    particleIndex =  mod( particleIndex, N);
    if particleIndex==0, particleIndex = N; end
    
    baseIndex = indices( particleIndex );
    
    % Get PDB number.
    ipdb = pdbID( baseIndex );
    
    
    if particles_{itype}.exchangeable
      sameMoleculeList = pdbID(baseIndex);
    else
      sameMoleculeList = pdb_.serial( pdb_.resSeq==pdb_.resSeq(ipdb) );
    end
    

    
    % Add particle copy of particle in this unit cell.
    resName = particles_{itype}.resName;
    for jpdb=sameMoleculeList'
      
      if ~isParticleWithinSystem(jpdb, cellShifts_(:,uc),true), continue; end
      
      % Get particleEnum and associatedParticles.
      preliminaryParticle = analyzeParticle(jpdb);
     
      if preliminaryParticle.particleEnum ~= particleEnum, continue; end
      
      % Define index to store unit cell and ipdb info together.
      ucpdb = jpdb + ucIndexOffset;
      
      % Add unit cell offset to associatedParticles. 
      associatedParticles = preliminaryParticle.associatedParticles ...
        + ucIndexOffset;
      
      % Since methyl carbons add their hydrons, some particle may already be
      % added.
      alreadyAdded = find(pdbID_==ucpdb);
      if ~isempty(alreadyAdded)

        if numel(alreadyAdded) > 1
          error(['Error in buidExtraCellss(): ', ...
            'pdb ID ',num2str(ucpbd),' added multiple times.']);
        end

        alreadyAddedID = particleClassID_(alreadyAdded);
        if particles_{alreadyAddedID}.particleEnum ...
            == preliminaryParticle.particleEnum
          error(['Error in buildExtraCells: ',...
            'cannot overwrite particle index ', num2str(alreadyAdded), ...
            ', ', getParticleString(particles_{alreadyAddedID}.particleEnum),...
            ' with ',getParticleString(preliminaryParticle.particleEnum), '.']);
        end

        continue;    
      end
      
      
      addParticle(particleEnum,resName, associatedParticles, number_+1, ucpdb);

    end
  end
end


end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function buildPrimaryCell()

numberMethyls_ = 0;

% Loop through all nuclei.
for ipdb = 1:pdb_.number
    
  % Ignore particles outside the system radius.
  if ~isParticleWithinSystem(ipdb, 0, false), continue; end
  
  preliminaryParticle = analyzeParticle(ipdb);
  
  if preliminaryParticle.particleEnum ~= PARTICLE_NONE
    % Add particle.
    resName = pdb_.resName{ipdb};
    addParticle(preliminaryParticle.particleEnum,resName,...
      preliminaryParticle.associatedParticles, number_+1, ipdb);

  end
end

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function buildSystem()
  
  % Initialize the observed spin.
  addParticle_new(PARTICLE_CENTRALSPIN, 'detect');
  particles_{numberParticleClasses_}.members = 0;
  particles_{numberParticleClasses_}.number = 1;
  residueList_{numberParticleClasses_} = 'detect';
  particleClassID_(numberParticleClasses_) = numberParticleClasses_;
  moleculeID_(numberParticleClasses_) = 0;
  
  number_ = number_ + 1;
  
  zeroIndex_ = number_;
  
  generateUniqueResidueList();
  
  buildPrimaryCell();
  
  buildExtraCells();

  buildAssociatedParticlematrix()

  % Adjust indices for dropped particles.
  %adjustAssociatedParticles();

  % Set methyls.
  setMethyls();

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function changeIsotope(iParticle)

% Get ParticleClass index number.
type_idx = particleClassID_(iParticle);  

% Get the new particle type.
newParticle = particles_{type_idx}.switchParticle;   

% Get the unit cell pdb number.
ucpdb = pdbID_(iParticle);

% Get unit cell index.
uc = ceil(ucpdb/pdb_.number);

% Check if the particle is outside the primary cell and has a different 
% switch particle defined for the extra cell.
% Note, void particle should have accounted for already. 
if uc > 1 ...
    && isParticleValid(particles_{type_idx}.extraCellSwitchParticle, ...
    particles_{type_idx}.resName) ...
    && particles_{type_idx}.extraCellSwitchParticle ~= PARTICLE_VOID
  
  newParticle = particles_{type_idx}.extraCellSwitchParticle;   
end

% Get the indices for any associated particles.
associatedParticles = ...
  particles_{type_idx}.associatedParticlesCollection{...
  particles_{type_idx}.members==iParticle};


% Remove particle from current ParticleClass.                               
if ~removeMember(type_idx, iParticle) 
  error(['Error in changeIsotope(): could not remove particle ', ...         
    num2str(ucpdb), ' from particle class ', ...                    
  getParticleString(particles_{type_idx}.particleEnum), ' ', ...
  particles_{type_idx}.resName, '.']);
end
resName = particles_{type_idx}.resName;      
addParticle(newParticle, resName, associatedParticles , iParticle, ucpdb);

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function allPass = checkNuclei()

METHYLID = 1;
pass(METHYLID) = true;

if System.Methyl.include
  for ii=1:Nuclei.number
    if Nuclei.MethylID(ii) > 0
      
      pass(METHYLID) = pass(METHYLID)*( ...
        sum(Nuclei.MethylID==Nuclei.MethylID(ii)) ==3);
      
    elseif Nuclei.MethylID(ii) < 0
      
      pass(METHYLID) = pass(METHYLID)*( ...
        sum(Nuclei.MethylID==Nuclei.MethylID(ii)) ==1);
      
    end

    if ~pass(METHYLID)
      error(['Error in checkNuclei(): ', 'methylID is not correctly set up.']);
    end
  end
end



  allPass = all(pass(METHYLID)); 
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% This function checks the system for errors.  More checks should be added as
% as bugs are are found and/or the code grows. 
function checkSystem()

  
n = 0;
nCmethyl = 0;
nHmethyl = 0;

allMembers = zeros(number_,1) - 1;

% Check that the number of particles matches the count from all particle
% classes.
for itype=1:numberParticleClasses_
  n0 = particles_{itype}.number;
  n1 = length(particles_{itype}.members);
  
  if n0 ~= n1
    error([ ...
      'Error in checkSystem(): ', 'particle class ', ...
      getParticleString(particles_{itype}.particleEnum), ' ',...
      particles_{itype}.resName, ' has' num2str(n1),' members, but reports ',...
      num2str(n0), ' members'...
      ])  ;
  end
  allMembers(n+1:n+n0) = particles_{itype}.members;
  n = n + n0;

  if particles_{itype}.particleEnum == PARTICLE_1H_METHYL
    nHmethyl = nHmethyl + particles_{itype}.number;
  elseif particles_{itype}.particleEnum == PARTICLE_C_METHYL
    nCmethyl = nCmethyl + particles_{itype}.number; 
  end 
end

% Check that the system has the correct number of particles.
if n ~= number_
    error([ ...
      'Error in checkSystem(): ', 'system has ', num2str(n),...
      ' members, but reports ', num2str(number_), ' members'...
      ])  ;
end

if 3*nCmethyl ~= nHmethyl
  error(['Error in checkSystem(): ', 'there ', num2str(nCmethyl), ...
    'methyl carbons, but ', num2str(nHmethyl), ' methy hydrodrons.']);  
end

% Check for multi-counted particles.
if length(unique(allMembers)) ~= number_
    error([ ...
      'Error in checkSystem(): ', 'system has ', ...
      num2str(length(unique(allMembers))),...
      ' unique members, but reports ', num2str(number_), ' members'...
      ])  ;
end

  
return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{
function ... % Nuclei = 
  cleanUpNuclei()%Nuclei,System,Method,Npdb)

if ~System.limitToSpinHalf
  Nuclei.Qtensor(:,:,Nuclei.number+1:end) = [];
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
%}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function clearBath()
number_ = zeroIndex_;                                                                     
particleClassID_ = particleClassID_(1:number_);                                                
moleculeID_ = moleculeID_(1:number_);                                                     
pdbID_ = pdbID_(1:number_);                                                          
residueList_ = residueList_(1:number_);  
 
numberMethyls_ = 0;
% methylIDs_ = [];
return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function clearMembers(index) 
  
particles_{index}.associatedParticlesCollection_ = {}; 
particles_{index}.members = [];
particles_{index}.membersOnePerMolecule = [];
particles_{index}.number = 0;
particles_{index}.numberMolecules = 0;

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{
function numberInPDB = countValidPDB(radius)

numberInPDB = 0;

% Loop through all nuclei.
for ipdb = 1:pdb_.number
  
  % Get pdb coordinate vector;
  rVec = [pdb_.x(ipdb); pdb_.y(ipdb); pdb_.z(ipdb) ];

  % Translate to system coordinates.
  rVec = rVec - originVec_;

  r = norm(rVec,2);
    
  % Ignore particles outside the system radius.
  if r > radius, continue; end
  numberInPDB = numberInPDB + 1;
end

return;
end
%}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function particleIndices = findParticles(particleEnum, resName, countOnly)
  
indices = findParticleClasses(particleEnum, resName);

if isempty(indices)
  if countOnly
    particleIndices  = 0;
  else
    particleIndices =  [];
  end
  return;
end

num = 0;

for itype = indices'
  num = num + particles_{itype}.number;  
end

if ~islogical(countOnly), countOnly = false; end

if countOnly
  particleIndices = num;
  return;
end

particleIndices = zeros(num,0);
num = 0;

for itype = indices
  particleIndices(num+1:num+particles_{itype}.number) = ...
    particles_{itype}.members;
  num = num + particles_{itype}.number;  
end

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function indices = findParticleClasses(particleEnum, resName)

indices = zeros(numberParticleClasses_,1) - 1;
for itype = 1:numberParticleClasses_
  doesMatch = true;
  if particleEnum ~= PARTICLE_ALL ...
      && particleEnum ~= particles_{itype}.particleEnum 
    doesMatch = false;
  end
  
  if ~strcmp(resName,'ALL') && ~strcmp(resName,particles_{itype}.resName)
    doesMatch = false;
  end
  
  if doesMatch
    indices(itype) = itype;
  end
  
end

indices(indices<0) = [];

return;  
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function generateUniqueResidueList()

increment = 10;   
currentSize = increment;
uniqueResidueList_= cell(increment,1);
resCounter = 0;

% Loop through all PDB atoms.  
for ipdb = 1:pdb_.number

  % Initialize test variable.
  isUnique = true;

  % Check all current molecule types for matches.
  for imol = 1:resCounter
    if uniqueResidueList_{imol} == pdb_.resName{ipdb}
      isUnique = false;
      break;
    end
  end
 
  % If the current type is unique, add it to the list.
  if isUnique
    resCounter = resCounter + 1;
    if resCounter>currentSize
      uniqueResidueList_{currentSize+increment} = '';
    end
    
    uniqueResidueList_{resCounter} = pdb_.resName{ipdb};
  end
end
uniqueResidueList_ = uniqueResidueList_(1:resCounter);

% Print list.
fprintf( ['Found ', num2str(length(uniqueResidueList_)),' unique types: { ']);
for imol = 1:length(uniqueResidueList_)-1                       
  fprintf([uniqueResidueList_{imol}, ', ']);                            
end   
fprintf( [uniqueResidueList_{end}, ' }.\n']);   
return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function atomEnum  = getAtomType(particleEnum)
switch particleEnum 
  case PARTICLE_UNSET
    atomEnum = particleEnum;
    
  case PARTICLE_ELECTRON
    atomEnum = particleEnum;
    
  case PARTICLE_HYDROGEN
    atomEnum = particleEnum;
    
  case PARTICLE_PROTIUM
    atomEnum = PARTICLE_HYDROGEN;
    
  case PARTICLE_1H_EXCHANGEABLE
    atomEnum = PARTICLE_HYDROGEN;
    
  case PARTICLE_1H_NONEXCHANGEABLE
    atomEnum = PARTICLE_HYDROGEN;
    
  case PARTICLE_1H_METHYL
    atomEnum = PARTICLE_HYDROGEN;
    
  case PARTICLE_DEUTERIUM
    atomEnum = PARTICLE_HYDROGEN;
    
  case PARTICLE_2H_EXCHANGEABLE
    atomEnum = PARTICLE_HYDROGEN;
    
  case PARTICLE_2H_NONEXCHANGEABLE
    atomEnum = PARTICLE_HYDROGEN;
    
  case PARTICLE_2H_METHYL
    atomEnum = PARTICLE_HYDROGEN;
    
  case PARTICLE_CARBON
    atomEnum = PARTICLE_CARBON;
    
  case PARTICLE_13C
    atomEnum = PARTICLE_CARBON;
    
  case PARTICLE_C_METHYL
    atomEnum = PARTICLE_CARBON;
    
  case PARTICLE_NITROGEN
    atomEnum = PARTICLE_NITROGEN;
    
  case PARTICLE_14N
    atomEnum = PARTICLE_NITROGEN;
    
  case PARTICLE_15N
    atomEnum = PARTICLE_NITROGEN;
    
  case PARTICLE_OXYGEN
    atomEnum = PARTICLE_OXYGEN;
    
  case PARTICLE_16O
    atomEnum = PARTICLE_OXYGEN;
    
  case PARTICLE_17O
    atomEnum = PARTICLE_OXYGEN;
    
  case PARTICLE_SILICON
    atomEnum = PARTICLE_SILICON;
    
  case PARTICLE_29SI
    atomEnum = PARTICLE_SILICON;
    
  case PARTICLE_NONE
    atomEnum = particleEnum;
    
  case PARTICLE_VOID
    atomEnum = particleEnum;
    
  case PARTICLE_CENTRALSPIN
    atomEnum = PARTICLE_ELECTRON;
    
  otherwise
    error('Error in getAtomType: could not identify the particle type.');
end 
return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function cellshifts = getCellShifts()

if ~System.isUnitCell
  cellshifts = zeros(1,3);
  return;
end

Coordinates = [pdb_.x,pdb_.y,pdb_.z] - originVec_';

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

if  ~System.isUnitCell 
  error(['System radius extends beyond the volume spanned by the pdb, ', ...
    'and no unit cell information is available.  ',...
    'If this is intentional, please set System.isUnitCell = false; ',...
    'otherwise, please add unit cell information to the pdb.']);  
end


if pdb_.a<0 || pdb_.b<0 || pdb_.c<0
  error(['Error in getCellShifts(): ', 'no unit cell is set.  ', ...
    'Please ensure the input PDB contains a properly formated CRYST1 ',...
    'header or set System.isUnitCell to false, to stop PBC.']);
end

ABC(1,:) = pdb_.a*[1,0,0];
ABC(2,:) = pdb_.b*[cos(pdb_.gamma),sin(pdb_.gamma),0];
cx = cos(pdb_.beta);
cy = ( cos(pdb_.alpha) - cos(pdb_.beta)*cos(pdb_.gamma) )/sin(pdb_.gamma);
cz = sqrt(1-cx^2 - cy^2);
ABC(3,:) = pdb_.c*[cx,cy,cz];


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
[~,indices] = sort(R);
cellshifts = cellshifts(indices,:);

if any(isnan(cellshifts(:)))
  error(['Error in getCellShifts(): ', 'cell shifts contain NaN values.']);
end

ucTest = (cellshifts(:,1)==cellshifts(:,1)').* ...
  (cellshifts(:,2)==cellshifts(:,2)').* ...
  (cellshifts(:,3)==cellshifts(:,3)') - eye(numUnitCells);

if max(abs(ucTest(:)))>Method.errorTolerance
  error(['Error in getCellShifts(): ', 'non-unique cell shifts created.']);
end

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function particleEnum = getParticleClass(particleStr)
    
if strcmp(particleStr, 'unset')
  particleEnum = PARTICLE_UNSET;
  
elseif strcmp(particleStr, 'e')
  particleEnum = PARTICLE_ELECTRON;
  
elseif strcmp(particleStr, 'hydrogen') || strcmp(particleStr, 'H')
  particleEnum = PARTICLE_HYDROGEN;
  
elseif strcmp(particleStr, '1H') || strcmp(particleStr, 'protium')
  particleEnum = PARTICLE_PROTIUM;
  
elseif strcmp(particleStr, '1H_exchangeable')
  particleEnum = PARTICLE_1H_EXCHANGEABLE;
  
elseif strcmp(particleStr, '1H_nonExchangeable')
  particleEnum = PARTICLE_1H_NONEXCHANGEABLE; 
  
elseif strcmp(particleStr, '1H_METHYL')
  particleEnum = PARTICLE_1H_METHYL;
  
elseif strcmp(particleStr, '2H') || strcmp(particleStr, 'deuterium')
  particleEnum = PARTICLE_DEUTERIUM;
  
elseif strcmp(particleStr, '2H_exchangeable')
  particleEnum = PARTICLE_2H_EXCHANGEABLE;
  
elseif strcmp(particleStr, '2H_nonExchangeable')
  particleEnum = PARTICLE_2H_NONEXCHANGEABLE;
  
elseif strcmp(particleStr, '2H_methyl') 
  particleEnum = PARTICLE_2H_METHYL;
  
elseif strcmp(particleStr, 'carbon') || strcmp(particleStr, 'C')
  particleEnum = PARTICLE_CARBON;
  
elseif strcmp(particleStr, '13C')
  particleEnum = PARTICLE_13C;
  
elseif strcmp(particleStr, 'C_methyl')
  particleEnum = PARTICLE_C_METHYL;
  
elseif strcmp(particleStr, 'nitrogen') || strcmp(particleStr, 'N')
  particleEnum = PARTICLE_NITROGEN;
  
elseif strcmp(particleStr, '14N')
  particleEnum = PARTICLE_14N;
  
elseif strcmp(particleStr, '15N')
  particleEnum = PARTICLE_15N;

elseif strcmp(particleStr, 'oxygen') || strcmp(particleStr, 'O')
  particleEnum = PARTICLE_OXYGEN;
  
elseif strcmp(particleStr, '17O')
  particleEnum = PARTICLE_17O;
  
elseif strcmp(particleStr, '29Si')
  particleEnum = PARTICLE_29SI;
  
elseif strcmp(particleStr, 'none')
  particleEnum = PARTICLE_NONE;
  
elseif strcmp(particleStr, 'void')
  particleEnum = PARTICLE_VOID;

elseif strcmp(particleStr, 'M')
  particleEnum = PARTICLE_VOID;

elseif strcmp(particleStr, 'S')
  particleEnum = PARTICLE_VOID;

elseif strcmp(particleStr, 'P')
  particleEnum = PARTICLE_VOID;  

elseif strcmp(particleStr, 'centralSpin')
  particleEnum = PARTICLE_CENTRALSPIN;
else
  error(['Error in getParticleClass: could not identify particle type "',...
    particleStr,'".']);
end
    
return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function nameStr = getParticleString(particleEnum)

switch particleEnum 
  case PARTICLE_UNSET
    nameStr = 'unset';

  case PARTICLE_ELECTRON
   nameStr = 'electron';

  case PARTICLE_HYDROGEN
    nameStr = 'hydrogen';

  case PARTICLE_PROTIUM
    nameStr = '1H';

  case PARTICLE_1H_EXCHANGEABLE
    nameStr = '1H_exchangeable';

  case PARTICLE_1H_NONEXCHANGEABLE
    nameStr = '1H_nonExchangeable';

  case PARTICLE_1H_METHYL
    nameStr = '1H_METHYL';
 
  case PARTICLE_DEUTERIUM
    nameStr = '2H';

  case PARTICLE_2H_EXCHANGEABLE
    nameStr = '2H_exchangeable';

  case PARTICLE_2H_NONEXCHANGEABLE
    nameStr = '2H_nonExchangeable';

  case PARTICLE_2H_METHYL
    nameStr = '2H_methyl';

  case PARTICLE_CARBON
    nameStr = '13C';

  case PARTICLE_13C
    nameStr = '13C';

  case PARTICLE_C_METHYL
    nameStr = 'C_methyl';

  case PARTICLE_NITROGEN
    nameStr = 'nitrogen';

  case PARTICLE_14N
    nameStr = '14N';

  case PARTICLE_15N
    nameStr = '15N';

  case PARTICLE_OXYGEN
    nameStr = 'oxygen';
    
  case PARTICLE_17O
    nameStr = '17O';

  case PARTICLE_SILICON
    nameStr = 'silicon';
    
  case PARTICLE_29SI
    nameStr = '29Si';

  case PARTICLE_NONE
    nameStr = 'none';

  case PARTICLE_VOID
    nameStr = 'void';

  case PARTICLE_CENTRALSPIN
    nameStr = 'centralSpin';
    
  otherwise
    error('Error in getParticleType: could not identify particle type.');
end

return;
end  
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% TO DO: THIS FUNCTION IS TEMPORARY.  CluE SHOULD BE RECONFIGURED TO WORK WITH
% THE DATA STRUCTURE OF centralSpinSystem(), WITHOUT A CONVERSION STEP.
function nameStr = getParticleTypeString(particleEnum)
  

switch particleEnum 
  case PARTICLE_UNSET
    nameStr = 'unset';

  case PARTICLE_ELECTRON
   nameStr = 'electron';

  case {PARTICLE_HYDROGEN,PARTICLE_PROTIUM,...
      PARTICLE_1H_EXCHANGEABLE,PARTICLE_1H_NONEXCHANGEABLE}
    nameStr = '1H';

  case PARTICLE_1H_METHYL
    nameStr = 'CH3_1H';
 
  case {PARTICLE_DEUTERIUM, ...
      PARTICLE_2H_EXCHANGEABLE,PARTICLE_2H_NONEXCHANGEABLE}
    nameStr = '2H';

  case PARTICLE_CARBON
    nameStr = '13C';

  case PARTICLE_13C
    nameStr = '13C';

  case PARTICLE_C_METHYL
    nameStr = 'CH3';

  case {PARTICLE_NITROGEN,PARTICLE_14N}
    nameStr = '14N';

  case {PARTICLE_SILICON,PARTICLE_29SI}
    nameStr = 'silicon';
     
  otherwise
    error('Error in getParticleTypeString: not a valid particle type.');
end

return;
end  
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function sameMoleculeList = getSameMoleculeList(iParticle, particleType, ...
    atomType)

if particleType == getAtomType(particleType)                            
  atomType = true;                                                               
elseif atomType
  particleType = getAtomType(particleType);
end

% Initialize output.
sameMoleculeList = zeros(number_,0);
counter = 0;

% Get molecule ID of input particle.
id = moleculeID_(iParticle);
if particleType==PARTICLE_ALL
  sameMoleculeList = moleculeID_(moleculeID_==id);
  return
end

% Loop through all particles.
for jParticle = 1:number_
  
  % Get molecule ID of the jth particle.
  testID = moleculeID_(jParticle);
  
  % Check if the particle is on the same molecule.
  if testID ~= id, continue; end
  
  testEnum = particles_{particleClassID_(jParticle)}.particleEnum;
  
  if atomType
    testEnum = getAtomType(testEnum);
  end
  
  if particleType == testEnum
    
    % Add the index to the list.
    counter = counter + 1;
    sameMoleculeList(counter) = jParticle;
  end
  
end
sameMoleculeList(counter+1:end) = [];

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function xyz = getTensorAxes(X,Y,Z,particle_index,particleEnum)
xyz = zeros(3);

xyz(:,1) = getTensorAxis(X,particle_index,particleEnum);
xyz(:,2) = getTensorAxis(Y,particle_index,particleEnum);
xyz(:,3) = getTensorAxis(Z,particle_index,particleEnum);

isSet = vecnorm(xyz) > 0;

if(sum(isSet)<2)
  error('Error in setTensorAxes(): at least two directions must be specified.');
end

if isSet(3) && isSet(1)
  xyz = orthonormalize(xyz, 3, 1);
 
elseif isSet(3) && isSet(2)
  xyz = orthonormalize(xyz, 3, 2);

elseif isSet(1) && isSet(2)
  xyz = orthonormalize(xyz, 1, 2);
  
end


err = max( max( abs(xyz'*xyz - eye(3)) ) );
if err > Method.errorTolerance
  error('Error in setTensorAxes(): could not set axes.')
end

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function x = getTensorAxis(X,particle_index,particleEnum)

x = zeros(3,1);  
sameMoleculeList = [];

if isempty(X)
  return;
  
elseif ~iscell(X) && ~ischar(X)
  x(1,1) = X(1);
  x(2,1) = X(2);
  x(3,1) = X(3);
  
elseif numel(X)>=2 && ~ischar(X)
  x(1,1) = pdb_.x(X{2}) - pdb_.x(X{1});
  x(2,1) = pdb_.y(X{2}) - pdb_.y(X{1});
  x(3,1) = pdb_.z(X{2}) - pdb_.z(X{1});
  
elseif ischar(X)
  
  if strcmp(X,'random')
  x(:,1) = rand(3,1);
  
  elseif length(X)>6 && strcmp(X(1:6),'bonded') 
    
    idx_pdb = mod(pdbID_(particle_index),pdb_.number);
    indices = find(pdb_.connections(:,idx_pdb));
    
    dashIdx = find(X=='-');
    
    selection = [str2double(X(7:dashIdx-1)), str2double(X(dashIdx+1:end))];
    selection(selection~=0) = indices(selection(selection~=0));
    selection(selection==0) = idx_pdb;
    
    x(1,1) = pdb_.x( selection(2) ) - pdb_.x( selection(1) );
    x(2,1) = pdb_.y( selection(2) ) - pdb_.y( selection(1) );
    x(3,1) = pdb_.z( selection(2) ) - pdb_.z( selection(1) );
    
    
  elseif length(X)>=4 && strcmp(X(1:4),'bond')
    
    idx_pdb = mod(pdbID_(particle_index),pdb_.number);
    pdb_connect_idx = find(pdb_.connections(:,idx_pdb));
    
     atomType = PARTICLE_ALL;
     spcIdx = find(X==' ');
     if length(X)>=5 &&  X(5)=='-'
       atomType = getAtomType(getParticleClass(X(6:spcIdx-1)));
     end

     if atomType ~= PARTICLE_ALL
       N = length(pdb_connect_idx);
       for iidx = 1:N
         pdb_connect_idx(iidx) = pdb_connect_idx(iidx)*...
           (atomType==getAtomType(getParticleClass(...
           pdb_.element(pdb_connect_idx(iidx)))));
       end
       pdb_connect_idx(pdb_connect_idx==0) = [];
     end

     N = length(pdb_connect_idx);
     if isempty(spcIdx)
       k = 1;
     else
       k = mod(str2double(X(spcIdx+1:end)),N);
       if k==0, k=N; end
     end
     if(numel(pdb_connect_idx) == 0)
      error(['Error in getTensorAxes(): ',...
        'could not identify a unique bond axis.']);
    end
    
    other_index = pdb_connect_idx(k);
    x(1,1) = pdb_.x(other_index) - pdb_.x(idx_pdb);
    x(2,1) = pdb_.y(other_index) - pdb_.y(idx_pdb);
    x(3,1) = pdb_.z(other_index) - pdb_.z(idx_pdb);
    
  elseif length(X)>4 && strcmp(X(1:4),'same') 
    if isempty(sameMoleculeList)

      atomType = PARTICLE_ALL;
      spcIdx = find(X==' ');
      if(X(5)=='-')
        atomType = getAtomType(getParticleClass(X(6:spcIdx-1)));
      end
      sameMoleculeList = getSameMoleculeList(particle_index,atomType,false);
      N = length(sameMoleculeList);
      n = find(sameMoleculeList == particle_index);
      k = mod(n + str2double(X(spcIdx+1:end)),N);
      if k==0, k=N; end
      x(1,1) = pdb_.x(k) - pdb_.x(pdbID_(particle_index));
      x(2,1) = pdb_.y(k) - pdb_.y(pdbID_(particle_index));
      x(3,1) = pdb_.z(k) - pdb_.z(pdbID_(particle_index));
    end
    
  end
  
end

xNorm = vecnorm(x);
if xNorm>0
  x = x/xNorm;
end

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function isExchangeable = isHydronExchangeable(ipdb)

% Initialize the output.
isExchangeable = false;

% Get rows.
rows = find(pdb_.connections(:,ipdb));


% Loop over connections.
for irow = 1:length(rows)
  
  % Get row index.
  row = rows(irow);
  
  % Check if the bond is to a hydrogen.
  if getParticleClass(pdb_.element{row}) == PARTICLE_OXYGEN ...
    || getParticleClass(pdb_.element{row}) == PARTICLE_NITROGEN
    isExchangeable = true;
  end
end

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function isSol = isMoleculeSolvent(resName)
  
isSol = false;  
N = length(System.solventMolecules);
  
for imol = 1:N
  isSol = isSol || strcmp(resName, System.solventMolecules{imol});
  if isSol
    break;
  end
end

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  function inSys = isParticleWithinSystem(ipdb, cellShift,useInnerCutoff)
  
% Get pdb coordinate vector;
rVec = [pdb_.x(ipdb); pdb_.y(ipdb); pdb_.z(ipdb) ];
  
% Translate to system coordinates.
rVec = rVec - originVec_ + cellShift;
  
r = norm(rVec,2);
  
% Ignore particles outside the system radius.
if r > System.load_radius 
  inSys  = false; 
else
  inSys = true;
end

if inSys && useInnerCutoff && r < System.inner_radius
      inSys  = false;
end

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function methylHydrons = isMethyl(ipdb)

% Initialize output for the not-a-methyl case.
methylHydrons = [0;0;0]; 

% Initialize the output for methylcarbons.
hydrogenIndexVec = [0;0;0];

% Get rows.
rows = find(pdb_.connections(:,ipdb));

% Initialize counters.
numBonded = 0;
numHydrogen = 0;

% Loop over connections.
for irow = 1:length(rows)
  
  % Count this bound.
  numBonded = numBonded + 1;
  
  % Methyls have 4 bonds, so exit if there are more than 4 bonds.
  if numBonded > 4, return; end
  
  % Get row index.
  row = rows(irow);
  
  % Check if the bond is to a hydrogen.
  particleEnum = getParticleClass(pdb_.element{row});
  if particleEnum == PARTICLE_HYDROGEN || ...
      particleEnum == PARTICLE_PROTIUM || ...
      particleEnum == PARTICLE_DEUTERIUM
    
    % Count this hydrogen.
    numHydrogen = numHydrogen + 1;
    
    % Record the hydrogen's index
    hydrogenIndexVec(numHydrogen) = row;
    
    % Methyls have 3 hydrogens, so exit if this carbon does not.
    if numHydrogen > 3, return; end
    
  end
end

% Methyl carbons  have 4 bonds and 3 hydrogens,
% so exit if this carbon does not.
if numHydrogen ~= 3 || numBonded ~=4 , return; end

% Print methyl information.
fprintf([ pdb_.element{ipdb}, ' pdb# ', num2str(ipdb), ...
  ' identified as a methyl carbon, with hydrons {', ...
  pdb_.element{hydrogenIndexVec(1) }, ' pdb# ', num2str(hydrogenIndexVec(1)), ...
  ', ', pdb_.element{hydrogenIndexVec(2)}, ...
  ' pdb# ', num2str(hydrogenIndexVec(2)), ...
  ', ',pdb_.element{hydrogenIndexVec(3)}, ...
  ' pdb# ', num2str(hydrogenIndexVec(3)), ... 
  '}.\n']);

% Return the PDB indices of the methyl hydrogens.  
methylHydrons = hydrogenIndexVec;
return;
end 
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function valid = isParticleValid(particleEnum, resName)
  
valid = false;

switch particleEnum 
  case PARTICLE_UNSET
    valid = false;

  case PARTICLE_ELECTRON
    valid = System.includeElectron;

  case PARTICLE_HYDROGEN
    valid = System.includeHydrogen;

  case PARTICLE_PROTIUM
    valid = System.include_1H;

  case PARTICLE_1H_EXCHANGEABLE
    valid = System.include_1H_exchangeable;

  case PARTICLE_1H_NONEXCHANGEABLE
    valid = System.include_1H_nonExchangeable;

  case PARTICLE_1H_METHYL
    valid = System.include_1H_METHYL;
 
  case PARTICLE_DEUTERIUM
    valid = System.include_2H;

  case PARTICLE_2H_EXCHANGEABLE
    valid = System.include_2H_exchangeable;

  case PARTICLE_2H_NONEXCHANGEABLE
    valid = System.include_2H_nonExchangeable;

  case PARTICLE_2H_METHYL
    valid = System.include_2H_methyl;

  case PARTICLE_CARBON
    valid = System.include_13C || System.include_C_methyl;

  case PARTICLE_13C
    valid = System.include_13C;

  case PARTICLE_C_METHYL
    valid = System.include_C_methyl;

  case PARTICLE_NITROGEN
    valid = System.include_14N || Nuclei.include_15N;

  case PARTICLE_14N
    valid = System.include_14N;

  case PARTICLE_15N
    valid = System.include_15N;

  case {PARTICLE_OXYGEN,PARTICLE_17O}
    valid = System.include_17O;

  case {PARTICLE_SILICON, PARTICLE_29SI}
    valid = System.include_29Si;

  case PARTICLE_NONE
    valid = false;

  case PARTICLE_VOID
    valid = true;

  case PARTICLE_CENTRALSPIN
    valid = true;
end



if strcmp('detect', resName) && particleEnum ~=PARTICLE_CENTRALSPIN
  valid = false;
end

if any(strcmp( resName, System.excludedResidues ))
  valid = false;
end

return;
end  
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% This function orthonormalizes a 3 by 3 matrix of column vectors.
% The column vector specified by the index 'primary' is taken directly as the 
% direction of that axis.  The secondary vector is modified by removing the
% cits component along the primary direction.  The final vector is the cross
% product of the primary and secondary vectors.
function xyz = orthonormalize(xyz, primary, secondary)
  
xyz(:,primary) = xyz(:,primary)/norm(xyz(:,primary));
  
xyz(:,secondary) = xyz(:,secondary) ...
  - xyz(:,primary)*(xyz(:,secondary)'*xyz(:,primary));

xyz(:,secondary) = xyz(:,secondary)/norm(xyz(:,secondary));


tertiary = 6/primary/secondary;
if mod(primary+1,3) == mod(secondary,3)
  sign = +1;
elseif  mod(primary-1,3) == mod(secondary,3)
  sign = -1;
else
  error(['Error in orthogonalize(): ',...
    'axis indices must be two unique numbers in {1,2,3}.']);
end

xyz(:,tertiary) = sign*cross(xyz(:,primary),xyz(:,secondary));

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function printParticleInfo( particle )
str = [ getParticleString(particle.particleEnum),', ', particle.resName,': '];

if particle.active  
  str  = [str,'active'];
else
  str  = [str,'inactive'];
end
disp(str);
if strcmp(particle.resName,'detect')
  printVector(particle.members, 'members');
else
  printVector(pdbID_(particle.members), 'members');
end

  
return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function success  = removeMember(address, index)
  
toRemove = find(particles_{address}.members==index);

if isempty(toRemove)
  err = ['Error in removeMember(): ', 'particle ' num2str(index),...
    ' is not a member of partilce class ', ...
    getParticleString(particles_{address}.particleEnum), ' ',...
    particles_{address}.resName, '.'];
  error(err)
elseif numel(toRemove)>1
  err = ['Error in removeMember(): ', 'particle ' num2str(index),...
    ' occurs ', num2str(numel(toRemove)), ' times in partilce class ', ...
    getParticleString(particles_{address}.particleEnum), ' ',...
    particles_{address}.resName, '.'];
  error(err)
end

particles_{address}.members(toRemove) = [];
particles_{address}.associatedParticlesCollection(toRemove) = [];
particles_{address}.number = length(particles_{address}.members);


particleClassID_(index) = -1;     

success = true;
return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{
function setAssociatedParticles(itype,iparticle,associatedParticles)

idx = 0;  
for imem = 1:particles_{itype}.number
  if iparticle == particles_{itype}.members(imem)
    idx=imem;
    break;
  end
end

if idx>0
  particles_{itype}.associatedParticlesCollection{idx} = associatedParticles;
else
  error(['Error in setAssociatedParticles(): ', 'particle ',...
    num2str(iparticle), ' is not a member of particle class ', ...
    particles_{itype}.name, ' ', particles_{itype}.resName, '.']);
end

return;
end
%}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% This function loops through all the nuclei and randomly generates a new
% isotopologue.
function setIsotopologue()

% Initialize vector of vectors of bools.
doSwitchIsotope = false(number_,1);

% Read loop, go through all initialized ParticleClasses.                
for itype = 1:numberParticleClasses_   

  % Skip particles set to be ignored.
  if ~particles_{itype}.doRandomIsotopes, continue; end

  % Void particles are accounted for when building the system.
  if particles_{itype}.particleEnum==PARTICLE_VOID, continue; end
  if particles_{itype}.switchParticle==PARTICLE_VOID, continue; end
  
  % Get number of particles of this type.
  nMember = particles_{itype}.number();
  
  % Loop through all particles of this type.
  for imem = 1:nMember

    if particles_{itype}.extraCellSwitchParticle==PARTICLE_VOID
      iParticle = particles_{itype}.members(imem);
      % Get unit cell index.
      uc = ceil(pdbID_(iParticle)/pdb_.number);
      if uc>1, continue; end
    end
    if rand() > particles_{itype}.abundance
    
      iParticle = particles_{itype}.members(imem);
      
      % Flag to switch.
      doSwitchIsotope(iParticle) = true;

    end

  end


end

% Write loop.
for iParticle=1:number_
  if doSwitchIsotope(iParticle)

    % Get ParticleClass index number.
    type_idx0 =  particleClassID_(iParticle);

    % If the particle is exchangeable, 
    if particles_{type_idx0}.exchangeable
      % then it can change isotopes independently.
      changeIsotope(iParticle);
    else

      % Get particle enum for error checking.
      particleType = particles_{type_idx0}.particleEnum;

      % Get list of all non-exchangeable hydrons on the same molecule.
      sameMolList = getSameMoleculeList(iParticle, particleType, false);
     
      % Only roll the die once per molecule.
      if(iParticle ~= sameMolList(1)), continue; end
      
      % Loop through all like-atoms of the molecule.
      for n = sameMolList
        % Get particle class index number.
        type_idx1 =  particleClassID_(n);
        
        % Change all non-exchangeable hydrons on this molecule.
        if ~particles_{type_idx1}.exchangeable
          % Get particle enum for error checking.
          testParticleType = particles_{type_idx1}.particleEnum;
          
          if particles_{type_idx0}.atomType == PARTICLE_NONE || ...
              particles_{type_idx0}.atomType == PARTICLE_UNSET
            error(['Error in setIsotopologue(): ', ...
              'atom type has not been set.']);
            
          elseif particleType == testParticleType
            changeIsotope(n);
            
          elseif particles_{type_idx0}.atomType ...
              == particles_{type_idx1}.atomType
            error( ['Error in setIsotopologue(): ', ...
              'in molecule ', num2str(moleculeID_(iParticle)), ...
              ',  particle ', num2str(iParticle), ', pdb# ', ...
              num2str(pdbID_(iParticle)),' is a ', ...
              getParticleString(particles_{type_idx0}), ...
              ', but particle ', num2str(n), ', pdb# ', num2str(pdbID_(n)),...
              ' is a ', getParticleString(particles_{type_idx1}), '.'] ) ;
            
          end
        end
      end
    end
  end
end
return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function setOptions()

deactivateIndices = [];

if System.Methyl.include
  if System.Methyl.method == 2
    deactivateIndices = findParticleClasses(PARTICLE_C_METHYL, 'ALL');
  else
    error(['Error in setOptions(): ', ...
      'methyl method ',  num2str(System.Methyl.method), ' not supported.']);
  end
end

for itype = deactivateIndices'
  particles_{itype}.active = false;
end

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function setParticle(particleIndex, particleEnum,resName)

switch particleEnum 
  case PARTICLE_UNSET

  case PARTICLE_ELECTRON
    particles_{particleIndex}.atomType = PARTICLE_ELECTRON;
    % https://www.physics.nist.gov/cgi-bin/cuu/Value?gem
    particles_{particleIndex}.gFactor = -2.00231930436256;
    particles_{particleIndex}.isNucleus = false;
    particles_{particleIndex}.spinMultiplicity = 2;
    particles_{particleIndex}.switchParticle = PARTICLE_VOID;

  case {PARTICLE_HYDROGEN, PARTICLE_PROTIUM}
    particles_{particleIndex}.atomType = PARTICLE_HYDROGEN;
    % https://physics.nist.gov/cgi-bin/cuu/Value?gp  
    particles_{particleIndex}.gFactor = 5.5856946893;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 2;
    particles_{particleIndex}.switchParticle = PARTICLE_DEUTERIUM;

  case PARTICLE_1H_EXCHANGEABLE
    particles_{particleIndex}.atomType = PARTICLE_HYDROGEN;
    % https://physics.nist.gov/cgi-bin/cuu/Value?gp  
    particles_{particleIndex}.gFactor = 5.5856946893;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 2;
    particles_{particleIndex}.switchParticle = PARTICLE_2H_EXCHANGEABLE;

  case PARTICLE_1H_NONEXCHANGEABLE
    particles_{particleIndex}.atomType = PARTICLE_HYDROGEN;
    particles_{particleIndex}.exchangeable = false;
    % https://physics.nist.gov/cgi-bin/cuu/Value?gp  
    particles_{numberParticleClasses_}.gFactor = 5.5856946893;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 2;
    particles_{particleIndex}.switchParticle = PARTICLE_2H_NONEXCHANGEABLE;

  case PARTICLE_1H_METHYL
    particles_{particleIndex}.atomType = PARTICLE_HYDROGEN;
    particles_{particleIndex}.exchangeable = false;
    % https://physics.nist.gov/cgi-bin/cuu/Value?gp  
    particles_{numberParticleClasses_}.gFactor = 5.5856946893;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 2;
    particles_{particleIndex}.switchParticle = PARTICLE_2H_METHYL;
 
  case PARTICLE_DEUTERIUM
    % https://physics.nist.gov/cgi-bin/cuu/Value?gdn 
    particles_{particleIndex}.gFactor = 0.8574382338;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 3;
    particles_{particleIndex}.switchParticle = PARTICLE_PROTIUM;

  case PARTICLE_2H_EXCHANGEABLE
    particles_{particleIndex}.atomType = PARTICLE_HYDROGEN;
    % https://physics.nist.gov/cgi-bin/cuu/Value?gdn 
    particles_{particleIndex}.gFactor = 0.8574382338;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 3;
    particles_{particleIndex}.switchParticle = PARTICLE_1H_EXCHANGEABLE;
    
    % Water Quadrupole Values
    % Edmonds, D. T.; Mackay, A. L.
    % The Pure Quadrupole Resonance of the Deuteron in Ice.
    % Journal of Magnetic Resonance (1969) 1975, 20 (3), 515–519.
    % https:%doi.org/10.1016/0022-2364(75)90008-6.
    particles_{particleIndex}.NQ_eta = 0.112;
    particles_{particleIndex}.NQ_e2qQh = 213.4e3; % Hz
    particles_{particleIndex}.NQ_z = 'bond-oxygen +1';
    if strcmp(resName,'WAT') || strcmp(resName,'SOL')
      particles_{particleIndex}.NQ_x = 'same-hydrogen +1';
    else
      particles_{particleIndex}.NQ_x = 'random';
    end
      

  case {PARTICLE_2H_NONEXCHANGEABLE,PARTICLE_2H_METHYL}
    particles_{particleIndex}.atomType = PARTICLE_HYDROGEN;
    particles_{particleIndex}.exchangeable = false;
    % https://physics.nist.gov/cgi-bin/cuu/Value?gdn 
    particles_{particleIndex}.gFactor = 0.8574382338;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 3;
    
    if particleEnum  == PARTICLE_2H_METHYL
      particles_{particleIndex}.switchParticle = PARTICLE_1H_METHYL;
    else
      particles_{particleIndex}.switchParticle = PARTICLE_1H_NONEXCHANGEABLE;
    end
    
    % ORCA
    particles_{particleIndex}.NQ_eta = 0; % from eta_ = 0.0161;
    particles_{particleIndex}.NQ_e2qQh = 0.1945e6; % Hz
    particles_{particleIndex}.NQ_x = 'random';
    particles_{particleIndex}.NQ_z = 'bond';

  case {PARTICLE_CARBON, PARTICLE_12C}
    particles_{particleIndex}.atomType = PARTICLE_CARBON;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 1;
    particles_{particleIndex}.switchParticle = PARTICLE_13C;
    
  case PARTICLE_13C  
    particles_{particleIndex}.atomType = PARTICLE_CARBON;
    % Weast, R. Handbook of Chemistry and Physics, 54th ed.; CRC Press, 1973.
    particles_{particleIndex}.gFactor = 1.4048; 
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 2;
    particles_{particleIndex}.switchParticle = PARTICLE_12C;

  case PARTICLE_C_METHYL
    particles_{particleIndex}.atomType = PARTICLE_CARBON;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 1;
    particles_{particleIndex}.switchParticle = PARTICLE_13C_METHYL;
    
  case PARTICLE_13C_METHYL
    particles_{particleIndex}.atomType = PARTICLE_CARBON;
    % Weast, R. Handbook of Chemistry and Physics, 54th ed.; CRC Press, 1973.
    particles_{particleIndex}.gFactor = 1.4048; 
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 2;
    particles_{particleIndex}.switchParticle = PARTICLE_C_METHYL;  

  case {PARTICLE_NITROGEN, PARTICLE_14N}
    particles_{particleIndex}.atomType = PARTICLE_NITROGEN;
    % Weast, R. Handbook of Chemistry and Physics, 54th ed.; CRC Press, 1973.
    particles_{particleIndex}.gFactor = 0.4036;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 3;
    particles_{particleIndex}.switchParticle = PARTICLE_15N;
    if strcmp(resName,'TEM')
      % Owenius, R.; Engström, M.; Lindgren, M.; Huber, M.
      % Influence of Solvent Polarity and Hydrogen Bonding on the EPR
      % Parameters of a Nitroxide Spin Label Studied by 9-GHz and
      % 95-GHz EPR Spectroscopy and DFT Calculations.
      %J. Phys. Chem. A 2001, 105 (49), 10967–10977.
      % https:%doi.org/10.1021/jp0116914.
      
      particles_{particleIndex}.hf_FermiContact = 31.528e+06; %Hz
      particles_{particleIndex}.hf_Azz = 90.801e+06; % Hz
      
      particles_{particleIndex}.hf_x = 'bonded 1-2'; % C--C
      particles_{particleIndex}.hf_y = 'bonded 0-3'; % N-O bond
      
      
      % de Oliveira, M.; Knitsch, R.; Sajid, M.; Stute, A.;
      % Elmer, L.-M.; Kehr, G.; Erker, G.; Magon, C. J.;
      % Jeschke, G.; Eckert, H.
      % Aminoxyl Radicals of B/P Frustrated Lewis Pairs:
      % Refinement of the Spin-Hamiltonian Parameters by Field- and
      % Temperature-Dependent Pulsed EPR Spectroscopy.
      % PLoS ONE 2016, 11 (6), e0157944.
      % https:%doi.org/10.1371/journal.pone.0157944.
      
      particles_{particleIndex}.NQ_e2qQh = 3.5*1e6; % Hz
      particles_{particleIndex}.NQ_eta = 0.68;
      
      particles_{particleIndex}.NQ_x = 'bonded 1-2'; % C--C
      particles_{particleIndex}.NQ_y = 'bonded 0-3'; % N-O bond
    end

  case PARTICLE_15N
    particles_{particleIndex}.atomType = PARTICLE_NITROGEN;
    % Weast, R. Handbook of Chemistry and Physics, 54th ed.; CRC Press, 1973.
    particles_{particleIndex}.gFactor = -0.5662;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 2;
    particles_{particleIndex}.switchParticle = PARTICLE_14N;

  case {PARTICLE_OXYGEN, PARTICLE_16O}
    particles_{particleIndex}.atomType = PARTICLE_OXYGEN;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 1;
    particles_{particleIndex}.switchParticle = PARTICLE_17O;
    
  case PARTICLE_17O
    particles_{particleIndex}.atomType = PARTICLE_OXYGEN;
    % Weast, R. Handbook of Chemistry and Physics, 54th ed.; CRC Press, 1973.
    particles_{particleIndex}.gFactor = -0.7575;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 6;
    particles_{particleIndex}.switchParticle = PARTICLE_16O;

  case {PARTICLE_SILICON, PARTICLE_28SI}
    particles_{particleIndex}.atomType = PARTICLE_SILICON;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 1;
    particles_{particleIndex}.switchParticle = PARTICLE_29SI;
    
  case PARTICLE_29SI
    particles_{particleIndex}.atomType = PARTICLE_SILICON;
    % Weast, R. Handbook of Chemistry and Physics, 54th ed.; CRC Press, 1973.
    particles_{particleIndex}.gFactor = 1.4447;
    particles_{particleIndex}.isNucleus = true;
    particles_{particleIndex}.spinMultiplicity = 2;
    particles_{particleIndex}.switchParticle = PARTICLE_28SI;

  case PARTICLE_NONE
    particles_{particleIndex}.atomType = PARTICLE_NONE;

  case PARTICLE_VOID
    particles_{particleIndex}.atomType = PARTICLE_VOID;

  case PARTICLE_CENTRALSPIN
    particles_{particleIndex}.atomType = PARTICLE_ELECTRON;
    
  otherwise
    error('Error in setParticle(): could not identify particle type.');
end

return;
end  
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function setParticleOptions(particleIndex,particleOptions)
  
particleEnum = particles_{particleIndex}.particleEnum;
atomEnum = particles_{particleIndex}.atomType;
resName = particles_{particleIndex}.resName;
particleStr = getParticleString(particleEnum);
atomStr = getParticleString(atomEnum);

Nopt = length(particleOptions);
for iopt=1:4:Nopt
  if strcmp(particleOptions{iopt}, particleStr) || ...
      strcmp(particleOptions{iopt}, atomStr)

    resBool = ~strcmp(particleOptions{iopt+1}(1),'!');
    if resBool
      optRes = particleOptions{iopt+1};
    else
      % Remove '!' from start.
      optRes = particleOptions{iopt+1}(2:end);
    end

    % Try to set option if resName is specified in particleOptions{iopt+1},
    % if resName not specified as avoided using '!', 
    % or if particleOptions{iopt+1}=='all'.
    if (strcmp(optRes, resName) == resBool) || ...
        (strcmp(optRes, 'all') && resBool)
      
      if strcmp(particleOptions{iopt+2},'abundance')
        particles_{particleIndex}.abundance = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'active')
        particles_{particleIndex}.active = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'associatedParticlesCollection')
        particles_{particleIndex}.associatedParticesCollection = ...
          particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'barrierPotential')
        particles_{particleIndex}.barrierPotential = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'doRandomIsotopes')
        particles_{particleIndex}.doRandomIsotopes = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'exchangeable')
        particles_{particleIndex}.exchangeable = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'hf_Azz')
        particles_{particleIndex}.hf_Azz = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'hf_FermiContact')
        particles_{particleIndex}.hf_FermiContact = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'hf_x')
        particles_{particleIndex}.hf_x = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'hf_y')
        particles_{particleIndex}.hf_y = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'hf_z')
        particles_{particleIndex}.hf_z = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'isNucleus')
        particles_{particleIndex}.isNucleus = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'NQ_e2qQh')
        particles_{particleIndex}.NQ_e2qQh = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'NQ_eta')
        particles_{particleIndex}.NQ_eta = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'NQ_x')
        particles_{particleIndex}.NQ_x = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'NQ_y')
        particles_{particleIndex}.NQ_y = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'NQ_z')
        particles_{particleIndex}.NQ_z = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'spinMultiplicity')
        particles_{particleIndex}.spinMultiplicity = particleOptions{iopt+3};
      elseif strcmp(particleOptions{iopt+2},'switchParticle')
        particles_{particleIndex}.switchParticle = ...
          getParticleClass(particleOptions{iopt+3});
      elseif strcmp(particleOptions{iopt+2},'extraCellSwitchParticle')
        particles_{particleIndex}.extraCellSwitchParticle = ...
          getParticleClass(particleOptions{iopt+3});
      elseif strcmp(particleOptions{iopt+2},'tunnelSplitting')
        particles_{particleIndex}.tunnelSplitting = particleOptions{iopt+3};
      end
      
    end
  end
end

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{
function  optionIDs = searchParticles(particleEnum, resName )

N = length(System.particleOptions);
optionIDs = zeros(N,1);
for ii = 1:N

  if strcmp(System.particleOptions{1}{1},getParticleString(particleEnum)) && ...
      strcmp(System.particleOptions{1}{2},resName)
    optionIDs(ii) = ii;
  end
end
optionIDs(optionIDs ==0) = [];

return;
end 
%}   
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%{
  function index = ucpdbID2Index(ucpdb)
  
  % Loop through all indices, starting in the correct unit cell.
  index = find(pdbID_ == ucpdb);
  
  if numel(index)==1
    return;
  end
  
  % Get unit cell index.
  uc = ceil(ucpdb/pdb_.number);

  % Get PDB number.
  ipdb = mod(ucpdb,pdb_.number);

  % The list should be complete.  Report error.
  if isempty(index)
    error(['Error in ucpdbID2Index(): could not find index for', ...
      num2str(ipdb), ' unit cell ', num2str(uc), '.']);
  else
    disp('index = ');
    disp(index);

    error(['Error in ucpdbID2Index(): found multiple indices for', ...
      num2str(ipdb), ' unit cell ', num2str(uc), '.']);
  end
end
%}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function refreshSystem()

% Set coordinates.  
setCoordinates();

% Set isotopologue
setIsotopologue();

% Set methyls.
% setMethyls();

% Set nuclear quadrupoles.
setQuadrupoles()

% Set nuclear quadrupoles.
setHyperfine()

% Set post setup options.
setOptions();

% Check for errors.
checkSystem();  
return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function setCoordinates()

% Set Coordinates to have 3 rows and a column for each particle.
coordinates_ = zeros(3,number_);

% Loop through all particle other than the central spin.
for iSpin = zeroIndex_+1:number_
  
  % Initialize coordinate vector.
  xyz = [0;0;0];
  
  % Get pdb ID.
  ucpdb = pdbID_(iSpin);

  % Get PDB number.
  ipdb = mod(ucpdb,pdb_.number);
  if ipdb==0, ipdb = pdb_.number; end

  % Get unit cell index.
  uc = ceil(ucpdb/pdb_.number);
  
  % Get x, y, and z coordinates.
  xyz(1) = pdb_.x(ipdb);
  xyz(2) = pdb_.y(ipdb);
  xyz(3) = pdb_.z(ipdb);

  % Translate system.
  coordinates_(:,iSpin) = xyz - originVec_ + cellShifts_(:,uc);
end  

coorTest = ...
  (coordinates_(1,zeroIndex_+1:number_)==...
  coordinates_(1,zeroIndex_+1:number_)').* ...
  (coordinates_(2,zeroIndex_+1:number_)==...
  coordinates_(2,zeroIndex_+1:number_)').* ...
  (coordinates_(3,zeroIndex_+1:number_)==...
  coordinates_(3,zeroIndex_+1:number_)')...
  - eye(number_-zeroIndex_);

if max(abs(coorTest(:))) > Method.errorTolerance
  error(['Error in setCoordinates(): ', ...
    'multiple particles have the same coordinates']);
end

return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function setHyperfine()
  
  for itype=1:numberParticleClasses_
    
    if particles_{itype}.hf_Azz==0 && particles_{itype}.hf_FermiContact == 0
      continue;
    end
    
    for imem=1:particles_{itype}.number
      
      R_A2L = getTensorAxes(...
        particles_{itype}.hf_x,particles_{itype}.hf_y,particles_{itype}.hf_z,...
        particles_{itype}.members(imem),particles_{itype}.particleEnum);
      
      
      Azz = particles_{itype}.hf_Azz;
      fc = particles_{itype}.hf_FermiContact;
      
      Atensor_A = eye(3)*fc + diag(Azz/2*[-1,-1,2]);
      Atensor_L = R_A2L'*Atensor_A*R_A2L;
      
      particles_{itype}.Atensor{imem} = Atensor_L;
    end
  end
  return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Rotate coordinates if requested by user via System.X/Y/Z
%-------------------------------------------------------------------------------
function setOrientation()

pdbCoordinates = [pdb_.x,pdb_.y,pdb_.z];
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
function setMethyls()

if ~System.Methyl.include
  return;
end

% Loop through all initialized ParticleClasses.
for itype = 1:numberParticleClasses_
   
  % Check if particle is a methyl.
  if particles_{itype}.particleEnum ~= PARTICLE_C_METHYL
    continue;
  end

  % Initialize hydron index to large value.
  methylHID = -1;
  
  % Determine if methyl hydrons are unique.
  isUnique = true;
  % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  %
  % ISSUE NOTES (DELETE THIS COMMENT ONCE ADDRESSED)
  %
  % DOUBLE CHECK THIS CODE ON A SYSTEM WITH MULTIPLE METHYL GROUPS,
  % AND WHEN THIS FUNCTION IS CALLED MORE THAN ONCE PER SESSION.
  % 
  % THERE SHOULD BE A WAY TO SIMPLIFY THIS SECTION: REPLACE THE isUnique
  % SECTIONS WITH A CALL TO A MORE GENERAL FUNCTION LIKE addParticle(),
  % OR A NEW FUNCTION FOR THAT PURPOSE.
  %
  % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  for jtype=1:numberParticleClasses_
     
    if itype==jtype, continue; end

    isUnique = isUnique && ...
      ~(particles_{jtype}.particleEnum == PARTICLE_1H_METHYL && ...
         strcmp(particles_{itype}.resName, particles_{jtype}.resName) ...
       );
       
    if ~isUnique
      % Record ID number for this particle.
      methylHID = jtype;
      break;
    end
    
  end

  % If this kind of methyl is new, assign a new ID number.
  if isUnique
    methylHID = numberParticleClasses_ + 1;
  end
  
  % Set hydron tunnel splitting to the methyl tunnel splitting.   
  particles_{methylHID}.tunnelSplitting = particles_{itype}.tunnelSplitting;
  
  % Get hydrons associated with this kind of metthyl.
  %   methylAssociatedParticlesCollection = ...
  %     particles_{itype}.associatedParticlesCollection;

  % Loop over methyls of this kind.
  %   for imethyl=1:length(methylAssociatedParticlesCollection)
  for imethyl=1:particles_{itype}.number

    % Get methyl ID within Particle Class.
    methylCID = particles_{itype}.members(imethyl);

    % Get the indices for the 3 methyl hydrons.
    %     methylHydrons = ...
    %       methylAssociatedParticlesCollection{imethyl};
    methylHydrons = nonzeros(associatedParticleMatrix_(:,methylCID));
    if numel(methylHydrons) ~= 3
      error(['Error in setMethyls(): ', ...
        'methyl carbon ', num2str(methylCID), ' has ', ...
        num2str(numel(methylHydrons)), ' instead of 3 hydrons.' ]);
    end

    % Loop over methyl hydrons.  
    for iHydron = 1:3

      % Get hydron ID number.
      hydronID = methylHydrons(iHydron);
 
      % Determine which ParticleClass the hydron is currently in.
      hydronAddress = particleClassID_(hydronID);
        
      % Get the current particle type.
      particleH = particles_{hydronAddress}.particleEnum;
      
      % Check the current type.
      if particleH==PARTICLE_1H_NONEXCHANGEABLE ...
          || particleH == PARTICLE_1H_METHYL
        methylEnum = PARTICLE_1H_METHYL;
      elseif particleH == PARTICLE_2H_NONEXCHANGEABLE
        disp(['There is a 2H on methyl group ', num2str(methylCID), ...
          '.  Methyl tunneling will be ignored.']);
        particles_{itype}.active = false;
        methylEnum = PARTICLE_2H_METHYL;
      else
        
        printParticleInfo( particles_{hydronAddress} );
                
        error(['Error in setMethyls(): methyl hydron, ', ...
          num2str(hydronAddress), ' catagorized as ',  ...
          particles_{hydronAddress}.resName, ' ',...
          getParticleString(hydronAddress), '.']);
      end


      % Remove hydron from current ParticleClass.
      if ~removeMember(hydronAddress, hydronID)
        
        printParticleInfo(particles_{hydronAddress});
        
        error(['Error in setMethyls(): could not remove hydron ', ...
          num2str(hydronID), ' from the following ParticleClass.']);       
        
      end

      % Add hydron to new ParticleClass.
      addParticle(methylEnum, particles_{itype}.resName,...
        methylCID, hydronID, pdbID_(hydronID));

      % Update address.
      particleClassID_(hydronID) = methylHID;

    end
  end
end
  
return;
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function setQuadrupoles()
  
  for itype=1:numberParticleClasses_
    
    if particles_{itype}.NQ_e2qQh ==0
      % ~particles_{itype}.isNucleus || particles_{itype}.spinMultiplicity < 3
      continue;
    end
    
    for imem=1:particles_{itype}.number
      
      R_Q2L = getTensorAxes(...
        particles_{itype}.NQ_x,particles_{itype}.NQ_y,particles_{itype}.NQ_z,...
        particles_{itype}.members(imem),particles_{itype}.particleEnum);
      
      I = (particles_{itype}.spinMultiplicity - 1)/2;
      
      e2qQh = particles_{itype}.NQ_e2qQh*System.nuclear_quadrupole_scale_e2qQh;
      eta = particles_{itype}.NQ_eta*System.nuclear_quadrupole_scale_eta;
      
      Qtensor_Q = e2qQh/4/I/(2*I-1)  *  diag([-1+eta, -1-eta, 2]);
      Qtensor_L = R_Q2L'*Qtensor_Q*R_Q2L;
      
      particles_{itype}.Qtensor{imem} = Qtensor_L;
    end
  end
  return;
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
function ... %Nuclei = 
  computeNuclearInteractions()%Nuclei,System, Method,scaleFactor)
  
Nuclei.Statistics = getPairwiseStatistics(System, Method, Nuclei);
Nuclei.DistanceMatrix = Nuclei.Statistics.DistanceMatrix;
if any(Nuclei.Statistics.Distance > System.radius*System.scale)
  if System.Methyl.include
    % TO DO: ADD CHECK.
  else
    error(['Error in parseNuclei(): ','Nuclei beyond the distance cutoff ', ...
      'remain in the system.'])
  end
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
% [Nuclei.State, ~]= setThermalEnsembleState(System,Nuclei);

Nuclei.ZeemanStates = setRandomZeemanState(Nuclei);
[Nuclei.RandomDenityMatrices,Nuclei.RandomSpinVector] = ...
  setRandomDensityMatrix(Nuclei);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end



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
%{
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
%}
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
function System = setCentralSpinSystemDefaults(System, Method)

if ~isfield(System,'isUnitCell')
  System.isUnitCell = true;
end

if ~isfield(System,'particleOptions')
  System.particleOptions = {};
end

if ~isfield(System,'excludedResidues')
  System.excludedResidues = {''};
end
if ~isfield(System,'includeElectron')
  System.include_electron = true;
end
if ~isfield(System,'includeHydrogen')
  System.include_hydrogen = true;
end
if ~isfield(System,'include_1H')
  System.include_1H = System.include_hydrogen;
end
if ~isfield(System,'include_1H_exchangeable')
  System.include_1H_exchangeable = System.include_1H;
end
if ~isfield(System,'include_1H_nonExchangeable')
  System.include_1H_nonExchangeable = System.include_1H;
end
if ~isfield(System,'include_1H_METHYL')
  System.include_1H_METHYL = System.include_1H;
end
if ~isfield(System,'include_2H')
  System.include_2H = System.include_hydrogen;
end
if ~isfield(System,'include_2H_exchangeable')
  System.include_2H_exchangeable = System.include_2H;
end
if ~isfield(System,'include_2H_nonExchangeable')
  System.include_2H_nonExchangeable = System.include_2H;
end
if ~isfield(System,'include_2H_methyl')
  System.include_2H_methyl = System.include_2H;
end
if ~isfield(System,'include_13C')
  System.include_13C = false;
end
if ~isfield(System,'include_C_methyl')
  System.include_C_methyl = System.Methyl.include;
end
if ~isfield(System,'include_14N')
  System.include_14N = true;
end
if ~isfield(System,'include_15N')
  System.include_15N = false;
end
if ~isfield(System,'include_17O')
  System.include_17O = false;
end
if ~isfield(System,'include_29Si')
  System.include_29Si = true;
end
    
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

if ~isfield(System, 'solventMolecules')
  System.solventMolecules = {'SOL','WAT','MGL'};
end

if ~isfield(System,'particleOptions')
  System.particleOptions = {};
end
for iopt=1:4:numel(System.particleOptions)
  if strcmp(System.particleOptions{iopt}, 'methyl')
    System.particleOptions{iopt} = 'C_methyl';
  end
end

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
if all( size(Electron_Coordinates)==[3,1]) 
  return;
elseif all( size(Electron_Coordinates)==[1,3])
  Electron_Coordinates = Electron_Coordinates';
else 
  error(['Error in getElectronCoordinates():', ...
    'could not get electron coordinates.']);
end

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
