function Isotopologue = newHydronIsotopologue(Nuclei,System)
NuclearList = 1:Nuclei.number;
isSolvent = Nuclei.isSolvent(NuclearList);

switch System.HydrogenExchange
  case 'OH'
    for iNuc = NuclearList
      
      if Nuclei.Exchangeable(iNuc)
        deuteriumFraction = System.deuteriumFraction;
      else
        deuteriumFraction = System.deuteriumFraction_nonExchangeable;
      end
      type = Nuclei.Type{iNuc};
      if isSolvent(iNuc) && (strcmp(type,'1H') || strcmp(type,'2H'))
        
        if rand() > deuteriumFraction
          Nuclei.Type{iNuc} = '1H';
          Nuclei.Spin(iNuc) = 0.5; % hbar
          Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
          Nuclei.Nuclear_g(iNuc) = 5.58569;
          Nuclei.NumberStates(iNuc) = int8(2);
          
          
          % D =============================================================
        else
          Nuclei.Type{iNuc} = '2H';
          Nuclei.Spin(iNuc) = 1; % hbar
          Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
          Nuclei.Nuclear_g(iNuc) = 0.857438;
          Nuclei.NumberStates(iNuc) = int8(3);
          
        end
        
        if ~Nuclei.Exchangeable(iNuc)
          sameMolecule = NuclearList(Nuclei.MoleculeID==Nuclei.MoleculeID(iNuc));
          
          for jNuc = sameMolecule 
            if Nuclei.Exchangeable(jNuc) || iNuc==jNuc
              continue;
            end
            Nuclei.Type{jNuc} = Nuclei.Type{iNuc};
            Nuclei.Spin(jNuc) = Nuclei.Spin(iNuc); % hbar
            Nuclei.StateMultiplicity(jNuc) = Nuclei.StateMultiplicity(iNuc);
            Nuclei.Nuclear_g(jNuc) = Nuclei.Nuclear_g(iNuc);
            Nuclei.NumberStates(jNuc) = Nuclei.NumberStates(iNuc);
          end
          
        end
        
      end
    end
  case 'full'
    for iNuc = 1:Nuclei.number
      
      type = Nuclei.Type{iNuc};
      if isSolvent(iNuc) && (strcmp(type,'1H') || strcmp(type,'2H'))
        
        if rand() > System.deuteriumFraction
          Nuclei.Type{iNuc} = '1H';
          Nuclei.Spin(iNuc) = 0.5; % hbar
          Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
          Nuclei.Nuclear_g(iNuc) = 5.58569;
          Nuclei.NumberStates(iNuc) = int8(2);
          
          
          % D =============================================================
        else
          Nuclei.Type{iNuc} = '2H';
          Nuclei.Spin(iNuc) = 1; % hbar
          Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
          Nuclei.Nuclear_g(iNuc) = 0.857438;
          Nuclei.NumberStates(iNuc) = int8(3);
          
        end
        
      end
    end
    
  case 'none'
    for iMol = Nuclei.MoleculeIDunique
      MolecularNuclei = find(Nuclei.MoleculeID== iMol);
      if isempty(MolecularNuclei)
        continue;
      end
      
      testNuc = MolecularNuclei(1);
      if isSolvent(testNuc)
        
        for iNuc = MolecularNuclei
          
          type = Nuclei.Type{iNuc};
          if isSolvent(iNuc) && (strcmp(type,'1H') || strcmp(type,'2H'))
            
            if rand() > System.deuteriumFraction
              Nuclei.Type{iNuc} = '1H';
              Nuclei.Spin(iNuc) = 0.5; % hbar
              Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
              Nuclei.Nuclear_g(iNuc) = 5.58569;
              Nuclei.NumberStates(iNuc) = int8(2);
              
              
              % D =============================================================
            else
              Nuclei.Type{iNuc} = '2H';
              Nuclei.Spin(iNuc) = 1; % hbar
              Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
              Nuclei.Nuclear_g(iNuc) = 0.857438;
              Nuclei.NumberStates(iNuc) = int8(3);
              
            end
            
          end
        end
      end
    end
  otherwise
    error('Could not generate isotopologue.');
end

is1H = cellfun(@(x)isequal(x,'1H'),Nuclei.Type) & isSolvent;
is2H = cellfun(@(x)isequal(x,'2H'),Nuclei.Type) & isSolvent;
isExch = Nuclei.Exchangeable;
n1H = sum(is1H); % number of solvent protons.
n2H = sum(is2H); % number of solvent deuterons.
n1Hx = sum(is1H & isExch); % number of exchangeable solvent protons.
n2Hx = sum(is2H & isExch); % number of exchangeable solvent deuterons. 
n1Hn = sum(is1H & ~isExch);  % number of non-exchangeable solvent protons.
n2Hn = sum(is2H & ~isExch);  % number of non-exchangeable solvent deuterons. 
nH = (n1H+n2H); % number of solvent hydrons.
nHx = n1Hx + n2Hx; % number of exchangeable solvent hydrons.
nHn = n1Hn + n2Hn; % number of non-exchangeable solvent hydrons.

Ptarget = nHx/nH*System.deuteriumFraction + nHn/nH*System.deuteriumFraction_nonExchangeable;
P = n2H/(n1H+n2H);
Px =  n2Hx/nHx;
Pn =  n2Hn/nHn;


fprintf('Target: \n  P(D) = %d/%d = %d. \n  P(D|OH) = %d/%d = %d.\n  P(D|CH) = %d/%d = %d.\n',...
  round(nHx*System.deuteriumFraction + nHn*System.deuteriumFraction_nonExchangeable), nH, Ptarget, ...
  round(nHx*System.deuteriumFraction), nHx, System.deuteriumFraction,...
  round(nHn*System.deuteriumFraction_nonExchangeable), nHn, System.deuteriumFraction_nonExchangeable);

fprintf('Initialized %d hydrons: \n  P(D) = %d/%d = %d. \n  P(D|OH) = %d/%d = %d.\n  P(D|CH) = %d/%d = %d.\n',nH,...
  n2H, (n1H+n2H), P, ...
  n2Hx, nHx, Px, ...
  n2Hn, nHn, Pn);


Isotopologue = Nuclei;
Isotopologue.Isotopologue.Name = {'All','OH','CH'};
Isotopologue.Isotopologue.TypeNumber = [nH,nHx,nHn];
Isotopologue.Isotopologue.TargetFraction = [Ptarget,System.deuteriumFraction,System.deuteriumFraction_nonExchangeable];
Isotopologue.Isotopologue.Instance_2H_Number = [n2H,n2Hx,n2Hn];
Isotopologue.Isotopologue.InstanceFraction = [P,Px,Pn];
end