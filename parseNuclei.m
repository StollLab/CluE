function Nuclei = parseNuclei(System, Method,datafile)

% set values to unspecified fields
System = setDefault(System);


% Define spin operators.
Nuclei.SpinOperators = cell(1,4);

multiplicity = 1;
Nuclei.SpinOperators{multiplicity}=1;
[Spin2Op_1,Spin2Op_2,Spin2Op_3,Spin2Op_4,Spin2Op_5,Spin2Op_6] = generateSpinOperators(1/2);
[Spin3Op_1,Spin3Op_2,Spin3Op_3,Spin3Op_4,Spin3Op_5,Spin3Op_6] = generateSpinOperators(1);
[Spin4Op_1,Spin4Op_2,Spin4Op_3,Spin4Op_4,Spin4Op_5,Spin4Op_6] = generateSpinOperators(3/2);

multiplicity = 2;
Nuclei.SpinOperators{multiplicity} = {Spin2Op_1,Spin2Op_2,Spin2Op_3,Spin2Op_4,Spin2Op_5,Spin2Op_6};
multiplicity = 3;
Nuclei.SpinOperators{multiplicity} = {Spin3Op_1,Spin3Op_2,Spin3Op_3,Spin3Op_4,Spin3Op_5,Spin3Op_6};
multiplicity = 4;
Nuclei.SpinOperators{multiplicity} = {Spin4Op_1,Spin4Op_2,Spin4Op_3,Spin4Op_4,Spin4Op_5,Spin4Op_6};

multiplicity = 3;
[SpinXiXjOp_1,SpinXiXjOp_2,SpinXiXjOp_3,SpinXiXjOp_4,SpinXiXjOp_5,SpinXiXjOp_6]= generateXiXjSpinOperators(1);
Nuclei.SpinXiXjOperators{multiplicity} = {SpinXiXjOp_1,SpinXiXjOp_2,SpinXiXjOp_3,SpinXiXjOp_4,SpinXiXjOp_5,SpinXiXjOp_6};

% ENUM
E = 1;  Z = 2; RAISE = 3; LOWER = 4;

% Copy graphCriterion to Nuclei.
Nuclei.graphCriterion = Method.graphCriterion;

% scale volume
scaleFactor = System.scale;

if strcmp(datafile,'user')
  Nuclei = System.Nuclei - System.Electron.Coordinates;
  Nuclei.kT = System.kT;
else
  
  Nuclei.dataSource = datafile;
  Nuclei.number = 0;
  % open data file
  [Coordinates,Type,UnitCell,Connected,Indices_nonWater,pdbID,numberH] = parsePDB(datafile,System);
  % Connected = formConnection(Connected_,Indices_nonWater);
  
  % number of PDB entries used
  Npdb = length(Type); 
  Nuclei.quadrupole2lab = zeros(3,3,numberH(2));
  Nuclei.Qtensor = zeros(3,3,numberH(2));
  Nuclei.quadrupoleXaxis = zeros(numberH(2),3);
  Nuclei.quadrupoleYaxis = zeros(numberH(2),3);
  Nuclei.quadrupoleZaxis = zeros(numberH(2),3);
  
if System.Methyl.include
  [Methyl_Data,Coordinates,Type,UnitCell,Connected,Indices_nonWater] = findMethyls(Coordinates,Type,UnitCell,Connected,Indices_nonWater);
  Nuclei.Methyl_Data = Methyl_Data;
else
  Nuclei.Methyl_Data = [];
end

% Get nuclear quadrupole parameters.
if System.nuclear_quadrupole
  
    % TEMPORARY VALUES FOR TESTING
    quadrupoleMoment = 0.00286*ones(1,Npdb)*System.barn;
    quadrupoleAxis = cell(1,Npdb);
    for ii =1:Npdb
      quadrupoleAxis{ii} = rand(3,1);
      quadrupoleAxis{ii} = quadrupoleAxis{ii}/norm(quadrupoleAxis{ii});
    end
    Vzz = -1.8*ones(1,Npdb)*0.97175e22;
    eta = 0.2*ones(1,Npdb);
    e2qQh = 220e3; % Hz
end

  % Initialize the number of unit cells o include along each direction.
  a_limit = 1;
  b_limit = 1;
  c_limit = 1;
  ABC = eye(3); % unit cell edge vectors.
  
  if ~isfield(System,'UnitCell') && UnitCell.isUnitCell
    System.UnitCell = UnitCell;
  end
  
  if (isfield(System,'UnitCell') && isfield(System,'radius'))
    Angles = System.UnitCell.Angles;
    
    ABC(1,:) = System.UnitCell.ABC(1)*[1,0,0];
    ABC(2,:) = System.UnitCell.ABC(2)*[cos(Angles(3)),sin(Angles(3)),0];
    cx = cos(Angles(2));
    cy = ( cos(Angles(1)) - cos(Angles(2))*cos(Angles(3)) )/sin(Angles(3));
    cz = sqrt(1-cx^2 - cy^2);
    ABC(3,:) = System.UnitCell.ABC(3)*[cx,cy,cz];
    
    a_limit = 1 + 2*ceil(System.radius/ABC(1,1)/scaleFactor);
    b_limit = 1 + 2*ceil(System.radius/ABC(2,2)/scaleFactor);
    c_limit = 1 + 2*ceil(System.radius/ABC(3,3)/scaleFactor);
    
  end
  if ~isfield(System,'radius')
    System.radius = inf;
  end
  

  Initial_Electron_Coordinates = System.Electron.Coordinates;

  if iscell(Initial_Electron_Coordinates)
    % find nuclei to average over
    replaceNuclei = [Initial_Electron_Coordinates{:}];

    % place the electron at the mean coordinates
    System.Electron.Coordinates = mean( Coordinates(replaceNuclei,:),1);
    
    % set initial electron coordinates to a 3-vector
    Initial_Electron_Coordinates = System.Electron.Coordinates;
    

    
  elseif length(System.Electron.Coordinates) == 1
 
    replaceNucleus = System.Electron.Coordinates;
    System.Electron.Coordinates = Coordinates(replaceNucleus,:);
    Initial_Electron_Coordinates = Coordinates(replaceNucleus,:);
    
  elseif length(System.Electron.Coordinates) == 2
 
    replaceNucleus1 = System.Electron.Coordinates(1);
    replaceNucleus2 = System.Electron.Coordinates(2);
    System.Electron.Coordinates = 0.5*( Coordinates(replaceNucleus1,:) +Coordinates(replaceNucleus2,:));
    Initial_Electron_Coordinates = System.Electron.Coordinates;
    
  end
  nucleiCounter = uint32(0);
  
  % initialize A edge vector of unit cell
  Delta_A = 0*ABC(1,:);
  
  % loop over x unit cell spacings
  for xx =1:a_limit
    
    % initialize B edge vector of unit cell
    Delta_B=0*ABC(2,:);
    
    % loop over y unit cell spacings
    for yy=1:b_limit
      
      % initialize C edge vector of unit cell
      Delta_C = 0*ABC(3,:);
      
      % loop over z unit cell spacings
      for zz =1:c_limit
        
        % 3-vector offset to put each nucleus in the correct unit cell
        Delta_R = Delta_A+Delta_B+Delta_C;
        
        ElectronCenteredCorrdinates = scaleFactor*(Coordinates + Delta_R - System.Electron.Coordinates);
        
        % loop over all nuclei
        for inucleus = 1:size(Type,2)
          type = Type{inucleus};
          
          % get nuclear connetion data
          try
            Conect = Connected{inucleus};
          catch
            Conect = {};
          end
          
          % set nuclear coordinates relative to the electron
          NuclearCorrdinates = ElectronCenteredCorrdinates(inucleus,:);
          
          % skip if the electron-nuclear separation is over the set cutoff
          if norm(NuclearCorrdinates)>System.radius
            continue;
          end
          if norm(NuclearCorrdinates)<System.inner_radius
            continue;
          end
          % switch nuclear type
 
          % H =============================================================
          if strcmp(type,'H') && System.protium
            nucleiCounter = nucleiCounter +1;
            Nuclei.Index(nucleiCounter) = nucleiCounter;
            Nuclei.Tile{nucleiCounter} = [sign(Delta_A(1))*norm(Delta_A)/norm(ABC(1,:) ),sign(Delta_B(2))*norm(Delta_B)/norm(ABC(2,:) ),sign(Delta_C(3))*norm(Delta_C)/norm(ABC(3,:) )];
            Nuclei.Type{nucleiCounter} = '1H';
            Nuclei.Element{nucleiCounter} = type;
            Nuclei.Connected{nucleiCounter} = Conect;
            Nuclei.Spin(nucleiCounter) = 0.5; % hbar
            Nuclei.StateMultiplicity(nucleiCounter) = 2*Nuclei.Spin(nucleiCounter) +1;
            Nuclei.Nuclear_g(nucleiCounter) = 5.58569;
            Nuclei.Coordinates((nucleiCounter),:) = NuclearCorrdinates;
            Nuclei.PDBCoordinates((nucleiCounter),:)= Coordinates(inucleus,:);
            Nuclei.pdbID(nucleiCounter) = pdbID(inucleus);
            Nuclei.NumberStates(nucleiCounter) = int8(2);
            Nuclei.valid(nucleiCounter)= true;
            if any(Indices_nonWater==inucleus)
              Nuclei.Abundance(nucleiCounter) = System.protiumFractionProtein;
            else
              Nuclei.Abundance(nucleiCounter) = System.protiumFractionSolvent;
            end
            
            if isfield(System.Electron,'P')
              x = Nuclei.Coordinates(nucleiCounter,1);
              y = Nuclei.Coordinates(nucleiCounter,2);
              z = Nuclei.Coordinates(nucleiCounter,3);
              c = (2*System.mu0/3)*System.Electron.g*System.muB*Nuclei.Nuclear_g(nucleiCounter)*System.muN/(2*pi*System.hbar);
              Nuclei.Hyperfine(nucleiCounter) = System.Electron.P(x,y,z)*c;
            elseif norm( Nuclei.Coordinates((nucleiCounter) ))< System.angstrom*1e-2
              Nuclei.Hyperfine(nucleiCounter) = 1420e6; % Hz.
            else
              Nuclei.Hyperfine(nucleiCounter) = 0;
            end
          % CH3_A =========================================================
          elseif strcmp(type,'CH3')
            nucleiCounter = nucleiCounter +1;
            Nuclei.Index(nucleiCounter) = nucleiCounter;
            Nuclei.Tile{nucleiCounter} = [sign(Delta_A(1))*norm(Delta_A)/norm(ABC(1,:) ),sign(Delta_B(2))*norm(Delta_B)/norm(ABC(2,:) ),sign(Delta_C(3))*norm(Delta_C)/norm(ABC(3,:) )];
            Nuclei.Type{nucleiCounter} = 'CH3';
            Nuclei.Element{nucleiCounter} = type;
            Nuclei.Connected{nucleiCounter} = Conect;
            Nuclei.Spin(nucleiCounter) = 1/2; % hbar
            Nuclei.StateMultiplicity(nucleiCounter) = 2*Nuclei.Spin(nucleiCounter) +1;
            Nuclei.Nuclear_g(nucleiCounter) = 5.58569;
            Nuclei.Coordinates((nucleiCounter),:) = NuclearCorrdinates;
            Nuclei.PDBCoordinates((nucleiCounter),:)= Coordinates(inucleus,:);
            Nuclei.pdbID(nucleiCounter) = pdbID(inucleus);
            Nuclei.NumberStates(nucleiCounter) = int8(8);
            Nuclei.valid(nucleiCounter)= true;
            
            ID_ref_ = find(Methyl_Data.ID(:,1)==inucleus);
            Nuclei.Group_ID{nucleiCounter} = Methyl_Data.ID(ID_ref_,2);
            Nuclei.Auxiliary_ID(nucleiCounter,:) = Methyl_Data.Hydron_ID{inucleus};
%             Nuclei.Auxiliary_Coordinates{nucleiCounter} = Methyl_Data.Hydron_Coordinates{inucleus};

            [Nuclei.State{nucleiCounter}, Nuclei.Abundance(nucleiCounter)] = getMethylState(System); 
            
            nucleiCounter = nucleiCounter +1;
            Nuclei.Index(nucleiCounter) = nucleiCounter;
            Nuclei.Tile{nucleiCounter} = [sign(Delta_A(1))*norm(Delta_A)/norm(ABC(1,:) ),sign(Delta_B(2))*norm(Delta_B)/norm(ABC(2,:) ),sign(Delta_C(3))*norm(Delta_C)/norm(ABC(3,:) )];
            Nuclei.Type{nucleiCounter} = 'CH3_1H';
            Nuclei.Element{nucleiCounter} = type;
            Nuclei.Connected{nucleiCounter} = Conect;
            Nuclei.Spin(nucleiCounter) = 0.5; % hbar
            Nuclei.StateMultiplicity(nucleiCounter) = 2*Nuclei.Spin(nucleiCounter) +1;
            Nuclei.Nuclear_g(nucleiCounter) = 5.58569;
            Nuclei.Coordinates((nucleiCounter),:) = ...
              scaleFactor*(Methyl_Data.Hydron_Coordinates{inucleus}(1,:) + Delta_R - System.Electron.Coordinates);
            Nuclei.PDBCoordinates((nucleiCounter),:)= Methyl_Data.Hydron_Coordinates{inucleus}(1,:);
            Nuclei.NumberStates(nucleiCounter) = int8(2);
            Nuclei.valid(nucleiCounter)= false;
            Nuclei.Abundance(nucleiCounter) = 1;
            
            nucleiCounter = nucleiCounter +1;
            Nuclei.Index(nucleiCounter) = nucleiCounter;
            Nuclei.Tile{nucleiCounter} = [sign(Delta_A(1))*norm(Delta_A)/norm(ABC(1,:) ),sign(Delta_B(2))*norm(Delta_B)/norm(ABC(2,:) ),sign(Delta_C(3))*norm(Delta_C)/norm(ABC(3,:) )];
            Nuclei.Type{nucleiCounter} = 'CH3_1H';
            Nuclei.Element{nucleiCounter} = type;
            Nuclei.Connected{nucleiCounter} = Conect;
            Nuclei.Spin(nucleiCounter) = 0.5; % hbar
            Nuclei.StateMultiplicity(nucleiCounter) = 2*Nuclei.Spin(nucleiCounter) +1;
            Nuclei.Nuclear_g(nucleiCounter) = 5.58569;
            Nuclei.Coordinates((nucleiCounter),:) = ...
              scaleFactor*(Methyl_Data.Hydron_Coordinates{inucleus}(2,:) + Delta_R - System.Electron.Coordinates);
            Nuclei.PDBCoordinates((nucleiCounter),:)= Methyl_Data.Hydron_Coordinates{inucleus}(2,:);
            Nuclei.NumberStates(nucleiCounter) = int8(2);
            Nuclei.valid(nucleiCounter)= false;
            Nuclei.Abundance(nucleiCounter) = 1;
            
            nucleiCounter = nucleiCounter +1;
            Nuclei.Index(nucleiCounter) = nucleiCounter;
            Nuclei.Tile{nucleiCounter} = [sign(Delta_A(1))*norm(Delta_A)/norm(ABC(1,:) ),sign(Delta_B(2))*norm(Delta_B)/norm(ABC(2,:) ),sign(Delta_C(3))*norm(Delta_C)/norm(ABC(3,:) )];
            Nuclei.Type{nucleiCounter} = 'CH3_1H';
            Nuclei.Element{nucleiCounter} = type;
            Nuclei.Connected{nucleiCounter} = Conect;
            Nuclei.Spin(nucleiCounter) = 0.5; % hbar
            Nuclei.StateMultiplicity(nucleiCounter) = 2*Nuclei.Spin(nucleiCounter) +1;
            Nuclei.Nuclear_g(nucleiCounter) = 5.58569;
            Nuclei.Coordinates((nucleiCounter),:) = ...
              scaleFactor*(Methyl_Data.Hydron_Coordinates{inucleus}(3,:) + Delta_R - System.Electron.Coordinates);
            Nuclei.PDBCoordinates((nucleiCounter),:)= Methyl_Data.Hydron_Coordinates{inucleus}(3,:);
            Nuclei.NumberStates(nucleiCounter) = int8(2);
            Nuclei.valid(nucleiCounter)= false;
            Nuclei.Abundance(nucleiCounter) = 1;

%             
%             nucleiCounter = nucleiCounter + 1;
%             Nuclei.Type{nucleiCounter} = 'CH3_Ea';;
%             Nuclei.Element{nucleiCounter} = 'CH3_Ea';
%             Nuclei.Coordinates((nucleiCounter),:) = NuclearCorrdinates*inf;
%             
%             nucleiCounter = nucleiCounter + 1;
%             Nuclei.Type{nucleiCounter}= 'CH3_Eb';
%             Nuclei.Element{nucleiCounter} = 'CH3_Eb';
%             Nuclei.Coordinates((nucleiCounter),:) = NuclearCorrdinates*inf;
            
%             Nuclei.SelectionRules{nucleiCounter} = eye(4);
%             Nuclei.Transform{nucleiCounter} = MethylTransform;
%             Nuclei.Projection{nucleiCounter} = zeros(4,8);
%             Nuclei.Projection{nucleiCounter} = Methyl_Data.Projection_A;
          % CH3_E =========================================================
          elseif false %strcmp(type,'CH3_Ea') || strcmp(type,'CH3_Eb')
            nucleiCounter = nucleiCounter +1;
            Nuclei.Index(nucleiCounter) = nucleiCounter;
            Nuclei.Tile{nucleiCounter} = [sign(Delta_A(1))*norm(Delta_A)/norm(ABC(1,:) ),sign(Delta_B(2))*norm(Delta_B)/norm(ABC(2,:) ),sign(Delta_C(3))*norm(Delta_C)/norm(ABC(3,:) )];
            Nuclei.Type{nucleiCounter} = type;
            Nuclei.Element{nucleiCounter} = type;
            Nuclei.Connected{nucleiCounter} = Conect;
            Nuclei.Spin(nucleiCounter) = 1/2; % hbar
            Nuclei.StateMultiplicity(nucleiCounter) = 2*Nuclei.Spin(nucleiCounter) +1;
            Nuclei.Nuclear_g(nucleiCounter) = 5.58569;
            Nuclei.Coordinates((nucleiCounter),:) = NuclearCorrdinates*inf;
            Nuclei.PDBCoordinates((nucleiCounter),:)= Coordinates(inucleus,:);
            Nuclei.NumberStates(nucleiCounter) = int8(2);
            ID_ref_ = find(Methyl_Data.ID(:,1)==inucleus);
            Nuclei.MethylID(nucleiCounter) = Methyl_Data.ID(ID_ref_,2);
            Nuclei.Hydron_Coordinates{nucleiCounter} = Methyl_Data.Hydron_Coordinates{Nuclei.MethylID};
            Nuclei.SelectionRules{nucleiCounter} = eye(2);
            Nuclei.Transform{nucleiCounter} = Methyl_Data.Transform;
            if strcmp(type,'CH3_Ea')
              Nuclei.Projection{nucleiCounter} = Methyl_Data.Projection_Ea;
            else
              Nuclei.Projection{nucleiCounter} = Methyl_Data.Projection_Eb;
            end
          % D =============================================================
          elseif strcmp(type,'D') && System.deuterium && ~System.limitToSpinHalf
            nucleiCounter = nucleiCounter +1;
            Nuclei.Index(nucleiCounter) = nucleiCounter;
            Nuclei.Tile{nucleiCounter} = [sign(Delta_A(1))*norm(Delta_A)/norm(ABC(1,:) ),sign(Delta_B(2))*norm(Delta_B)/norm(ABC(2,:) ),sign(Delta_C(3))*norm(Delta_C)/norm(ABC(3,:) )];
            Nuclei.Type{nucleiCounter} = '2H';
            Nuclei.Element{nucleiCounter} = type;
            Nuclei.Connected{nucleiCounter} = Conect;
            Nuclei.Spin(nucleiCounter) = 1; % hbar
            Nuclei.StateMultiplicity(nucleiCounter) = 2*Nuclei.Spin(nucleiCounter) +1;
            Nuclei.Nuclear_g(nucleiCounter) = 0.857438;
            Nuclei.Coordinates((nucleiCounter),:) = NuclearCorrdinates;
            Nuclei.PDBCoordinates((nucleiCounter),:)= Coordinates(inucleus,:);
            Nuclei.pdbID(nucleiCounter) = pdbID(inucleus);
            Nuclei.NumberStates(nucleiCounter) = int8(3);
            Nuclei.valid(nucleiCounter)= true;

            Nuclei.Abundance(nucleiCounter) = 1;
            if isfield(System.Electron,'P')
              x = Nuclei.Coordinates(nucleiCounter,1);
              y = Nuclei.Coordinates(nucleiCounter,2);
              z = Nuclei.Coordinates(nucleiCounter,3);
              c = (2*System.mu0/3)*System.Electron.g*System.muB*Nuclei.Nuclear_g(nucleiCounter)*System.muN/(2*pi*System.hbar);
              Nuclei.Hyperfine(nucleiCounter) = System.Electron.P(x,y,z)*c;
            elseif norm(Nuclei.Coordinates((nucleiCounter) )) < System.angstrom*1e-2
              
              Nuclei.Hyperfine(nucleiCounter) = 220e6; % Hz.
            else
              Nuclei.Hyperfine(nucleiCounter) = 0;
            end
            
            if System.nuclear_quadrupole
              for iconnect = Conect
                switch Type{iconnect}
                  case 'O'
                    Nuclei.quadrupoleZaxis(nucleiCounter,:) = ElectronCenteredCorrdinates(iconnect,:) - NuclearCorrdinates;
                    Nuclei.quadrupoleZaxis(nucleiCounter,:) = Nuclei.quadrupoleZaxis(nucleiCounter,:)/norm(Nuclei.quadrupoleZaxis(nucleiCounter,:));
                  case 'M'
                    Nuclei.quadrupoleXaxis(nucleiCounter,:) = ElectronCenteredCorrdinates(iconnect,:) - NuclearCorrdinates;
                  otherwise
                end
              end
              if norm(Nuclei.quadrupoleZaxis(nucleiCounter,:)) == 0 || norm(Nuclei.quadrupoleXaxis(nucleiCounter,:)) < 0
                error('Failed to set quadrupole tensor orientation.')
              end
              Nuclei.quadrupoleXaxis(nucleiCounter,:) = Nuclei.quadrupoleXaxis(nucleiCounter,:) - Nuclei.quadrupoleXaxis(nucleiCounter,:)*Nuclei.quadrupoleZaxis(nucleiCounter,:)'*Nuclei.quadrupoleZaxis(nucleiCounter,:);
              Nuclei.quadrupoleXaxis(nucleiCounter,:) = Nuclei.quadrupoleXaxis(nucleiCounter,:)/norm(Nuclei.quadrupoleXaxis(nucleiCounter,:));

              Nuclei.quadrupoleYaxis(nucleiCounter,:) = cross(Nuclei.quadrupoleZaxis(nucleiCounter,:),Nuclei.quadrupoleXaxis(nucleiCounter,:));
              Nuclei.quadrupoleYaxis(nucleiCounter,:) = Nuclei.quadrupoleYaxis(nucleiCounter,:)/norm(Nuclei.quadrupoleYaxis(nucleiCounter,:));
              
              
              Nuclei.quadrupole2lab(:,:,nucleiCounter) = [Nuclei.quadrupoleXaxis(nucleiCounter,:); ...
                                       Nuclei.quadrupoleYaxis(nucleiCounter,:); ...
                                       Nuclei.quadrupoleZaxis(nucleiCounter,:)];
              
                                     
              Nuclei.Qtensor(:,:,nucleiCounter) = ... %quadrupoleMoment(inucleus)*Vzz(inucleus)*System.e/4/Nuclei.Spin(nucleiCounter)/(2*Nuclei.Spin(nucleiCounter) - 1) ...
                e2qQh/4/Nuclei.Spin(nucleiCounter)/(2*Nuclei.Spin(nucleiCounter) - 1) ...    
                *[-1 + eta(nucleiCounter), 0 ,0; 0, -1 - eta(nucleiCounter), 0; 0, 0, 2];
              
              % Rotate from qudrupole egienframe to lab frame:
              %  P_{lab} = R*P_{q}*R'.
              
              Nuclei.Qtensor(:,:,nucleiCounter) = Nuclei.quadrupole2lab(:,:,nucleiCounter)*Nuclei.Qtensor(:,:,nucleiCounter)*Nuclei.quadrupole2lab(:,:,nucleiCounter)';
              
              % Nuclei.HNQ will soon be obsolete.  Please use Nuclei.Qtensor instead.   
%               HNQ = ... % quadrupoleMoment(inucleus)*Vzz(inucleus)*System.e/4/Nuclei.Spin(nucleiCounter)/(2*Nuclei.Spin(nucleiCounter) - 1) ...
%                     e2qQh/4/Nuclei.Spin(nucleiCounter)/(2*Nuclei.Spin(nucleiCounter) - 1) ...
%                     *( ... 
%                        2*Nuclei.SpinOperators{3}{1}(:,:,Z)^2 ...
%                        - eta(inucleus)*(Nuclei.SpinOperators{3}{1}(:,:,RAISE)^2 + Nuclei.SpinOperators{3}{1}(:,:,LOWER)^2) ...
%                        - 1/2 * acommute(Nuclei.SpinOperators{3}{1}(:,:,RAISE), Nuclei.SpinOperators{3}{1}(:,:,LOWER) ) ...
%                      );
%                    
%               Rot = getInverseSpinRotationMatrix( quadrupoleAxis{inucleus},Nuclei.SpinOperators{3}{1}(:,:,Z), ...
%                                                   -1i/2*( Nuclei.SpinOperators{3}{1}(:,:,RAISE) - Nuclei.SpinOperators{3}{1}(:,:,LOWER) )   ...
%                                                  );
%               
%               Nuclei.HNQ{nucleiCounter} = Rot'*HNQ*Rot;
% %               Nuclei.HNQ{nucleiCounter} = Nuclei.HNQ{nucleiCounter}/System.hbar/2/pi; 
%               Nuclei.HNQ{nucleiCounter} = (Nuclei.HNQ{nucleiCounter}+Nuclei.HNQ{nucleiCounter}')/2;
            end
            
            % C ============================================================
          elseif strcmp(type,'C') && System.carbon
            nucleiCounter = nucleiCounter +1;
            Nuclei.Index(nucleiCounter) = nucleiCounter;
            Nuclei.Tile{nucleiCounter} = [sign(Delta_A(1))*norm(Delta_A)/norm(ABC(1,:) ),sign(Delta_B(2))*norm(Delta_B)/norm(ABC(2,:) ),sign(Delta_C(3))*norm(Delta_C)/norm(ABC(3,:) )];
            Nuclei.Type{nucleiCounter} = '13C';
            Nuclei.Element{nucleiCounter} = type;
            Nuclei.Connected{nucleiCounter} = Conect;
            Nuclei.Spin(nucleiCounter) = 0.5; % hbar
            Nuclei.StateMultiplicity(nucleiCounter) = 2*Nuclei.Spin(nucleiCounter) +1;
            Nuclei.Nuclear_g(nucleiCounter) = 1.4048;
            Nuclei.Coordinates((nucleiCounter),:) = NuclearCorrdinates;
            Nuclei.PDBCoordinates((nucleiCounter),:)= Coordinates(inucleus,:);
            Nuclei.pdbID(nucleiCounter) = pdbID(inucleus);
            Nuclei.NumberStates(nucleiCounter) = int8(2);
            Nuclei.valid(nucleiCounter)= true;
            
            Nuclei.Abundance(nucleiCounter) = 0.0107;
            
            Nuclei.Hyperfine(nucleiCounter) = 0;
            
  
          % N =============================================================  
          elseif strcmp(type,'N') && System.nitrogen  && ~System.limitToSpinHalf
            nucleiCounter = nucleiCounter +1;
            Nuclei.Index(nucleiCounter) = nucleiCounter;
            Nuclei.Tile{nucleiCounter} = [sign(Delta_A(1))*norm(Delta_A)/norm(ABC(1,:) ),sign(Delta_B(2))*norm(Delta_B)/norm(ABC(2,:) ),sign(Delta_C(3))*norm(Delta_C)/norm(ABC(3,:) )];
            Nuclei.Type{nucleiCounter} = '14N';
            Nuclei.Element{nucleiCounter} = type;
            Nuclei.Connected{nucleiCounter} = Conect;
            Nuclei.Spin(nucleiCounter) = 1; % hbar
            Nuclei.StateMultiplicity(nucleiCounter) = 2*Nuclei.Spin(nucleiCounter) +1;
            Nuclei.Nuclear_g(nucleiCounter) = 0.403761;
            Nuclei.Coordinates((nucleiCounter),:) = NuclearCorrdinates;
            Nuclei.PDBCoordinates((nucleiCounter),:)= Coordinates(inucleus,:);
            Nuclei.pdbID(nucleiCounter) = pdbID(inucleus);
            Nuclei.NumberStates(nucleiCounter) = int8(3);
            Nuclei.valid(nucleiCounter)= true;
            Nuclei.Abundance(nucleiCounter) = 0.99632;
            if isfield(System.Electron,'P')
              x = Nuclei.Coordinates(nucleiCounter,1);
              y = Nuclei.Coordinates(nucleiCounter,2);
              z = Nuclei.Coordinates(nucleiCounter,3);
              c = (2*System.mu0/3)*System.Electron.g*System.muB*Nuclei.Nuclear_g(nucleiCounter)*System.muN/(2*pi*System.hbar);
              Nuclei.Hyperfine(nucleiCounter) = System.Electron.P(x,y,z)*c;
            elseif norm(Nuclei.Coordinates((nucleiCounter))) < System.angstrom*1e-2
              
              Nuclei.Hyperfine(nucleiCounter) = 1811e6; % Hz.
            else
              Nuclei.Hyperfine(nucleiCounter) = 0;
            end
          % Si ============================================================  
          elseif strcmp(type,'Si') && System.silicon
            nucleiCounter = nucleiCounter +1;
            Nuclei.Index(nucleiCounter) = nucleiCounter;
            Nuclei.Tile{nucleiCounter} = [sign(Delta_A(1))*norm(Delta_A)/norm(ABC(1,:) ),sign(Delta_B(2))*norm(Delta_B)/norm(ABC(2,:) ),sign(Delta_C(3))*norm(Delta_C)/norm(ABC(3,:) )];
            Nuclei.Type{nucleiCounter} = '29Si';
            Nuclei.Element{nucleiCounter} = type;
            Nuclei.Connected{nucleiCounter} = Conect;
            Nuclei.Spin(nucleiCounter) = 0.5; % hbar
            Nuclei.StateMultiplicity(nucleiCounter) = 2*Nuclei.Spin(nucleiCounter) +1;
            Nuclei.Nuclear_g(nucleiCounter) = -1.11058;
            Nuclei.Coordinates((nucleiCounter),:) = NuclearCorrdinates;
            Nuclei.PDBCoordinates((nucleiCounter),:)= Coordinates(inucleus,:);
            Nuclei.pdbID(nucleiCounter) = pdbID(inucleus);
            Nuclei.NumberStates(nucleiCounter) = int8(2);
            Nuclei.valid(nucleiCounter)= true;
            
            Nuclei.Abundance(nucleiCounter) = 0.046832;
            if isfield(System.Electron,'P')
              x = Nuclei.Coordinates(nucleiCounter,1);
              y = Nuclei.Coordinates(nucleiCounter,2);
              z = Nuclei.Coordinates(nucleiCounter,3);
              c = (2*System.mu0/3)*System.Electron.g*System.muB*Nuclei.Nuclear_g(nucleiCounter)*System.muN/(2*pi*System.hbar);
              Nuclei.Hyperfine(nucleiCounter) = System.Electron.P(x,y,z)*c;
            elseif norm(Nuclei.Coordinates((nucleiCounter))) < System.angstrom*1e-2
              
              Nuclei.Hyperfine(nucleiCounter) = 2.96292e6; % Hz.
            else
              Nuclei.Hyperfine(nucleiCounter) = 0;
            end
            
          % electron ======================================================  
          elseif strcmp(type,'e')
            
            % electron, not a nucleus
            nucleiCounter = nucleiCounter +1;
            Nuclei.Index(nucleiCounter) = nucleiCounter;
            Nuclei.Tile{nucleiCounter} = [sign(Delta_A(1))*norm(Delta_A)/norm(ABC(1,:) ),sign(Delta_B(2))*norm(Delta_B)/norm(ABC(2,:) ),sign(Delta_C(3))*norm(Delta_C)/norm(ABC(3,:) )];
            Nuclei.Type{nucleiCounter} = 'e';
            Nuclei.Element{nucleiCounter} = type;
            Nuclei.Connected{nucleiCounter} = Conect;
            Nuclei.Spin(nucleiCounter) = 0.5; % hbar
            Nuclei.StateMultiplicity(nucleiCounter) = 2*Nuclei.Spin(nucleiCounter) +1;
            Nuclei.Nuclear_g(nucleiCounter) = 2.0023;
            Nuclei.Coordinates((nucleiCounter),:) = NuclearCorrdinates;
            Nuclei.PDBCoordinates((nucleiCounter),:)= Coordinates(inucleus,:);
            Nuclei.pdbID(nucleiCounter) = pdbID(inucleus);
            Nuclei.NumberStates(nucleiCounter) = int8(2);
            Nuclei.valid(nucleiCounter)= true;
            Nuclei.Abundance(nucleiCounter) = 1;
            Nuclei.Hyperfine(nucleiCounter) = 0;
             
          elseif System.allAtoms
            nucleiCounter = nucleiCounter +1;
            Nuclei.Index(nucleiCounter) = nucleiCounter;
            Nuclei.Tile{nucleiCounter} = [sign(Delta_A(1))*norm(Delta_A)/norm(ABC(1,:) ),sign(Delta_B(2))*norm(Delta_B)/norm(ABC(2,:) ),sign(Delta_C(3))*norm(Delta_C)/norm(ABC(3,:) )];
            Nuclei.Type{nucleiCounter} = type;
            Nuclei.Element{nucleiCounter} = type;
            Nuclei.Connected{nucleiCounter} = Conect;
            Nuclei.Spin(nucleiCounter) = 0; % hbar
            Nuclei.StateMultiplicity(nucleiCounter) = 2*Nuclei.Spin(nucleiCounter) +1;
            Nuclei.Nuclear_g(nucleiCounter) = 0;
            Nuclei.Coordinates((nucleiCounter),:) = NuclearCorrdinates;
            Nuclei.PDBCoordinates((nucleiCounter),:)= Coordinates(inucleus,:);
            Nuclei.NumberStates(nucleiCounter) = int8(2);
            Nuclei.valid(nucleiCounter)= true;
            Nuclei.Abundance(nucleiCounter) = 1;
            Nuclei.Hyperfine(nucleiCounter) = 0;
          end

        end
        
        % increment  C edge vector of unit cell
        Delta_C = -Delta_C + 0.5*(1-sign(Delta_C(3))+1-abs(sign(Delta_C(3))) )*ABC(3,:);
      end
      
      % increment  B edge vector of unit cell
      Delta_B = -Delta_B + 0.5*(1-sign(Delta_B(2))+1-abs(sign(Delta_B(2))) )*ABC(2,:);
    end
    % increment  A edge vector of unit cell
    Delta_A = -Delta_A + 0.5*(1-sign(Delta_A(1))+1-abs(sign(Delta_A(1))) )*ABC(1,:);
  end
  
end

% translate origin to electron
System.Electron.Coordinates = [0,0,0];

% get number of nuclei
try
  Nuclei.number = uint32(size(Nuclei.Index,2));
catch
  return;
end

Nuclei.DistanceMatrix = zeros(Nuclei.number);

for ispin = 1:Nuclei.number
  for jspin = ispin+1:Nuclei.number
    Nuclei.DistanceMatrix(ispin,jspin) = norm(Nuclei.Coordinates(ispin,:) - Nuclei.Coordinates(jspin,:));
    Nuclei.DistanceMatrix(jspin,ispin) = Nuclei.DistanceMatrix(ispin,jspin);
  end
end
Nuclei = getPairwiseStatistics(System, Nuclei);
Nuclei.maxSpin = max(Nuclei.Spin);

Nuclei.ValidPair = Nuclei.valid'*Nuclei.valid;
Nuclei.ValidPair = Nuclei.ValidPair.*Nuclei.Same_g;
% exclude_spins = find( isnan(Nuclei.Spin));
% Nuclei.ValidPair(exclude_spins,:) = 0;
% Nuclei.ValidPair(:,exclude_spins) = 0;

num_criteria = numel(Method.Criteria);
for ii = 1:num_criteria
  switch Method.Criteria{ii}
    case 'neighbor'
      Max_Neighbor_ = Nuclei.DistanceMatrix <= Method.r0;
      Min_Neighbor_ = Nuclei.DistanceMatrix > Method.r_min;
      Nuclei.ValidPair = Nuclei.ValidPair.*Max_Neighbor_.*Min_Neighbor_;
  
    case 'modulation'
      Min_Mod_ = Nuclei.ModulationDepth >= Method.cutoff.modulation;
      Nuclei.ValidPair = Nuclei.ValidPair.*Min_Mod_;
    
    case 'dipole'
      Min_dipole_ = abs(Nuclei.Nuclear_Dipole) >= Method.cutoff.dipole;
      Nuclei.ValidPair = Nuclei.ValidPair.*Min_dipole_;
  
    case 'minimum-frequency'  
      Min_Freq_ = abs(Nuclei.Frequency_Pair) >  Method.cutoff.minimum_frequency;
      Nuclei.ValidPair = Nuclei.ValidPair.*Min_Freq_;
  
    case 'hyperfine'
      Min_DeltaA_ = abs(Nuclei.DeltaHyperfine) > Method.cutoff.hyperfine_inf;
      Max_DeltaA_ = abs(Nuclei.DeltaHyperfine) < Method.cutoff.hyperfine_sup;
      Nuclei.ValidPair = Nuclei.ValidPair.*Min_DeltaA_.*Max_DeltaA_;
  end
  
end



Nuclei.startSpin = max(1, floor(Method.startSpin));
Nuclei.endSpin = min(Nuclei.number, floor(Method.endSpin));

if Nuclei.startSpin > Nuclei.endSpin
  disp('Starting cluster spin cannot be greater than ending spin.  Swapping assignment.');
  Nuclei.startSpin = max(0, floor(Method.endSpin));
  Nuclei.endSpin = min(Nuclei.number, floor(Method.startSpin));
end

Nuclei.numberStartSpins = min(Nuclei.number, Nuclei.endSpin - Nuclei.startSpin + 1);

% check if unit vectors have been specified by the user
if (isfield(System,'X')&&isfield(System,'Z')) || (isfield(System,'X')&&isfield(System,'Y')) || (isfield(System,'Y')&&isfield(System,'Z'))
  
  % translate unit vector to the elecron centered frame
  if isfield(System,'Z')&&(length(System.Z)==1)
    System.Z = Coordinates(System.Z,:)*1e-10 - Initial_Electron_Coordinates;
  end
  if isfield(System,'Y')&&(length(System.Y)==1)
    System.Y = Coordinates(System.Y,:)*1e-10 - Initial_Electron_Coordinates;
  end
  if isfield(System,'X')&&(length(System.X)==1)
    System.X = Coordinates(System.X,:)*1e-10 - Initial_Electron_Coordinates;
  end
  
  % find the x and z unit vectors if one was not specified
  if ~isfield(System,'X')
    System.X = cross(System.Y,System.Z);
  elseif ~isfield(System,'Z')
    System.Z = cross(System.X,System.Y);
  end
  
  % rotate system
  Rotation = alignCoordinates(System.X',System.Z');
  Nuclei.Coordinates = Nuclei.Coordinates*Rotation';
  
end

% set thermal energy
Nuclei.kT = System.kT;

% set thermal equilibrium state
[Nuclei.State, Nuclei.ZeemanStates]= setState(System,Nuclei);

end

% ========================================================================
% New Function
% ========================================================================

function [Coordinates,Type,UnitCell,Connected, Indices_nonWater pdbID,numberH] = parsePDB(filename,System)
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
          continue;
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

% ========================================================================
% New Function
% ========================================================================

function R = alignCoordinates(MolecularX,MolecularZ)

NormZ = MolecularZ/norm(MolecularZ);
NormX = MolecularX/norm(MolecularX);
if abs(NormX'*NormZ) > 1e-12
  %warning('Coordinate axes are not orthogonal.  Removing the z-component from the x-direction...')
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
Z = R*MolecularZ;
Z = Z/norm(Z);
X = R*MolecularX;
X = X/norm(X);
X = X - (X'*Z)*Z;
X = X/norm(X);
Z = abs(Z - [0;0;1]);
X = abs(X - [1;0;0]);
if (sum(X) + sum(Z)) > 1e-12
  disp(sum(X) + sum(Z));
  error('Coordinates not aligned.');
end

end

% ========================================================================
% New Function
% ========================================================================

function Rx = rotateX(gamma)
Rx  =   [1,  0,            0; ...
  0,  cos(gamma),  -sin(gamma); ...
  0,  sin(gamma),   cos(gamma )];
end

% ========================================================================
% New Function
% ========================================================================

function Ry = rotateY(theta)
Ry =   [ cos(theta),  0,         sin(theta); ...
  0,           1,         0;          ...
  -sin(theta),  0,         cos(theta)];
end

% ========================================================================
% New Function
% ========================================================================

function Rz = rotateZ(phi)
Rz  =   [cos(phi ), -sin(phi ), 0; ...
  sin(phi ),  cos(phi ), 0; ...
  0,          0,         1];
end

% ========================================================================
% Set thermal states.
% ========================================================================

% Sets the initial bath state with a Boltzmann distribution.
% Both output variables contain the same information, but are formated
% differently.
function [State,ZeemanStates] = setState(System,Nuclei)

ZeemanStates = zeros(Nuclei.number, 2*Nuclei.maxSpin+1);

% Loop through all nuclei.
for iinucleus = Nuclei.number:-1:1 
 
  spin_multiplicity = Nuclei.NumberStates(iinucleus);
  
  if isfield(Nuclei,'State') && (length(Nuclei.State) >= iinucleus) && ~isempty(Nuclei.State{iinucleus})
  
    % Assign state. 
    State{iinucleus} = Nuclei.State{iinucleus};
    
    ZeemanStates(iinucleus,1:length(Nuclei.State{iinucleus})) = Nuclei.State{iinucleus};
    continue;
  
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
    deltaE = -mI*System.magneticField*System.muN*Nuclei.Nuclear_g(iinucleus)- en0;
    
    
    % Set population according to the Boltzmann distribution.
    State{iinucleus}(ithresh) = State{iinucleus}(ithresh) + exp(-deltaE/Nuclei.kT/2);
    ZeemanStates(iinucleus,ithresh)  = ZeemanStates(iinucleus,ithresh) + exp(-deltaE/Nuclei.kT/2);
  end
  
  % Normalize states.
  normalization = sqrt(State{iinucleus}'*State{iinucleus});
  State{iinucleus} = State{iinucleus}/normalization; 
  ZeemanStates(iinucleus,:) = ZeemanStates(iinucleus,:)./normalization;
  
end
end

% ========================================================================
% New Function
% ========================================================================

function Connected = formConnection(Connected_,Indices_nonWater)
N = length(Indices_nonWater);
if N<1
  Connected = {};
  return;
end
Connected{N} = {};
if  isempty(Connected_)
  return;
end
for ii=1:Indices_nonWater
  if  isempty(Connected_{ii})
    continue;
  end
  referenceNuclei = Connected_{ii}(1);
  Connected{referenceNuclei} = [Connected{referenceNuclei}, Connected_{ii}(2:end)];
  Connected{referenceNuclei} = unique(Connected{referenceNuclei});
  for jj = 1:length(Connected{referenceNuclei})
    conenctTo = Connected{referenceNuclei}(jj);
    Connected{conenctTo} = [Connected{conenctTo},referenceNuclei];
    if iscell(Connected{conenctTo})
      Connected{conenctTo} = Connected{conenctTo}{:};
    end
    Connected{conenctTo} = unique(Connected{conenctTo});
  end
end

end

% ========================================================================
% New Function
% ========================================================================

function System = setDefault(System)
if ~isfield(System,'protium')
  System.protium = true;
end
if ~isfield(System,'protiumFractionSolvent')
  System.protiumFractionSolvent = 1;
end
if ~isfield(System,'protiumFractionProtein')
  System.protiumFractionProtein = 1;
end
if ~isfield(System,'deuterium')
  System.deuterium = true;
end
if isfield(System,'hydrogen')
  System.protium = System.hydrogen;
  System.dueterium = System.hydrogen;
end
if ~isfield(System,'carbon')
  System.carbon = true;
end
if ~isfield(System,'nitrogen')
  System.nitrogen = true;
end
if ~isfield(System,'oxygen')
  System.oxygen = false;
end
if ~isfield(System,'silicon')
  System.silicon = true;
end
if ~isfield(System,'allAtoms')
  System.allAtoms = false;
end
if ~isfield(System,'scale')
  System.scale = 1;
end

end

% ========================================================================
% New Function
% ========================================================================

function  [Methyl_Data,Coordinates,Type,UnitCell,Connected,Indices_nonWater] = findMethyls(Coordinates,Type,UnitCell,Connected,Indices_nonWater)

Number_Nuclei = length(Type);
Methyl_Data.number_methyls = 0;
Methyl_Data.Hydron_Coordinates = cell(1);
Methyl_Data.ID = [];

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
      indexH = indexH + (indexH >=ii);
      break;
    elseif strcmp(methyl,'HCHH')
      indexH = Connected{ispin}([ 1,3,4 ]);
      indexH = indexH + (indexH >=ii);
      break;
    elseif strcmp(methyl,'HHCH')
      indexH = Connected{ispin}([ 1,2,4]);
      indexH = indexH + (indexH >=ii);
      break;
    elseif strcmp(methyl,'HHHC')
      indexH = Connected{ispin}([ 1,2,3 ]);
      indexH = indexH + (indexH >=ii);
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
  Coordinates(new_index,:) = Center_of_Mass_Coor;
  for iH = 1:3
    % Save hydron coordinates for Hamiltonian construction.
    Methyl_Data.Hydron_Coordinates{new_index}(iH,:) = Coordinates(indexH(iH),:);
  end
  % Change hydrogens into methy pseudo-particles.
  Type{new_index} = 'CH3';
  Type{indexH(1)} = 'CH3_H_source';
  Type{indexH(2)} = 'CH3_H_source';
  Type{indexH(3)} = 'CH3_H_source';
  
  
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
  Methyl_Data.Projection_5,Methyl_Data.Projection_6] = getMethylProjections(P,E);

%{
PEE = zeros(8);
for ii=1:4
  PEE = PEE +Methyl_Data.Transform(ii,:)'*Methyl_Data.Transform(ii,:);
end
PEAB = zeros(8);
PEBA = zeros(8);
PEBA = PEBA +Methyl_Data.Transform(1,:)'*Methyl_Data.Transform(4,:);
PEAB = PEAB +Methyl_Data.Transform(4,:)'*Methyl_Data.Transform(1,:);
PEAB = PEAB +Methyl_Data.Transform(2,:)'*Methyl_Data.Transform(3,:);
PEBA = PEBA +Methyl_Data.Transform(3,:)'*Methyl_Data.Transform(2,:);

PE0 = zeros(8);
for ii=1:4
  for jj=1:4
    PE0=PE0 +Methyl_Data.Transform(ii,:)'*Methyl_Data.Transform(jj,:);
  end
end

PEa = zeros(8);
for ii=1:2
  for jj=1:2
    PEa=PEa +Methyl_Data.Transform(ii,:)'*Methyl_Data.Transform(jj,:);
  end
end
PEb = zeros(8);
for ii=3:4
  for jj=3:4
    PEb=PEb +Methyl_Data.Transform(ii,:)'*Methyl_Data.Transform(jj,:);
  end
end
%}

end

% ========================================================================
% New Function
% ========================================================================

function [State, gA]= getMethylState(System)

Nucleus_.kT=System.kT;
Nucleus_.Nuclear_g(1) = 5.58569;
Nucleus_.Spin(1) = 3/2;
Nucleus_.maxSpin = 3/2;
Nucleus_.NumberStates(1) = int8(4);
Nucleus_.number = 1;
State_A = setState(System,Nucleus_);

Nucleus_.Spin(1) = 1/2;
Nucleus_.NumberStates(1) = int8(2);
State_E = setState(System,Nucleus_);

Spin_Density = [State_E{1}; State_E{1}; State_A{1}].^2;
e_EkT= exp(-System.hbar*System.Methyl.tunnel_splitting/System.Methyl.kT);
Density = Spin_Density;%.*[e_EkT;e_EkT;e_EkT;e_EkT; 1;1;1;1];
Density = Density/sum(Density);
State = sqrt(Density);
gA = 1/(1+e_EkT);
end

% ========================================================================
% New Function
% ========================================================================

function Nuclei = getPairwiseStatistics(System, Nuclei)

Hyperfine = zeros(1,Nuclei.number);

Nuclei.ModulationDepth = zeros(Nuclei.number);
Nuclei.Hyperfine = zeros(Nuclei.number,1);
Nuclei.Nuclear_Dipole = zeros(Nuclei.number);
Nuclei.Frequency_Pair = zeros(Nuclei.number);
Nuclei.DeltaHyperfine = zeros(Nuclei.number);
Nuclei.SameType = eye(Nuclei.number);
Nuclei.Same_g = eye(Nuclei.number);

Nuclei.SpinSet = unique(Nuclei.Spin);

% loop over all nuclei
for inucleus = 1:Nuclei.number
  

  % Calculating hyperfine coupling
  
  gamma_e = -System.Electron.g*System.muB/System.hbar;
  Rn = norm(System.Electron.Coordinates - Nuclei.Coordinates(inucleus,:) );
  gamma_n = Nuclei.Nuclear_g(inucleus)*System.muN/System.hbar;

  R = System.Electron.Coordinates - Nuclei.Coordinates(inucleus,:);
  cosTheta2 = (R*[0;0;1]/norm(R))^2;
  
  Hyperfine(inucleus) = Nuclei.Hyperfine(inucleus)-(System.mu0/4/pi)*gamma_n*gamma_e*System.hbar*(1-3*cosTheta2)*Rn^-3;
  Nuclei.Hyperfine(inucleus) = Hyperfine(inucleus)/(2*pi); % Hz;
  % Calculating bath coupling
  for jnucleus = 1:inucleus-1 
    
    if strcmp(Nuclei.Type{inucleus},Nuclei.Type{jnucleus})
      Nuclei.SameType(inucleus,jnucleus) = true;
      Nuclei.SameType(jnucleus,inucleus) = true;
    end
    
    if abs(Nuclei.Nuclear_g(inucleus)-Nuclei.Nuclear_g(jnucleus)) < 1e-9
      Nuclei.Same_g(inucleus,jnucleus) = true;
      Nuclei.Same_g(jnucleus,inucleus) = true;
    end
    
    % calculate dipolar coupling
    delta_r = Nuclei.Coordinates(inucleus,:) -Nuclei.Coordinates(jnucleus,:);
    cosThetaSquared = ( delta_r(3)/norm(delta_r) )^2;
    b = 0.25*(System.mu0/4/pi)*Nuclei.Nuclear_g(inucleus)*Nuclei.Nuclear_g(jnucleus)*System.muN^2; % J m^3.
    r = norm(Nuclei.Coordinates(inucleus,:) - Nuclei.Coordinates(jnucleus,:));
    r3 = r^3;
    b = -b*(3*cosThetaSquared - 1)/r3; % J.
    b = b/(System.hbar); % rad/s.

    
    c = ( Hyperfine(inucleus)-Hyperfine(jnucleus) )/(4*b);
    
    mod_amp = 4*c^2/(1+c^2)^2;
    
    Nuclei.Frequency_Pair(inucleus,jnucleus) = 2*b*sqrt(1+c^2)/(2*pi);
    Nuclei.Frequency_Pair(jnucleus,inucleus)=    Nuclei.Frequency_Pair(inucleus,jnucleus);
    
    Nuclei.ModulationDepth(inucleus,jnucleus) = mod_amp;
    Nuclei.ModulationDepth(jnucleus,inucleus) = Nuclei.ModulationDepth(inucleus,jnucleus);
    
    Nuclei.Nuclear_Dipole(inucleus,jnucleus) = 4*b/(2*pi); % Hz.
    Nuclei.Nuclear_Dipole(jnucleus,inucleus) = Nuclei.Nuclear_Dipole(inucleus,jnucleus);
    
    Nuclei.DeltaHyperfine(inucleus,jnucleus) = (Nuclei.Hyperfine(inucleus) - Nuclei.Hyperfine(jnucleus)); % Hz
    Nuclei.DeltaHyperfine(jnucleus,inucleus) = -Nuclei.DeltaHyperfine(inucleus,jnucleus);
  end
  
end
end
