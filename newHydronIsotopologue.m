function Isotopologue = newHydronIsotopologue(Nuclei,System)
isSolvent = Nuclei.isSolvent;

switch System.HydrogenExchange
  case 'OH'
    NuclearList = 1:Nuclei.number;
    for iNuc = NuclearList
      
      if Nuclei.Exchangable(iNuc)
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
        
        if ~Nuclei.Exchangable(iNuc)
          sameMolecule = NuclearList(Nuclei.MoleculeID==Nuclei.MoleculeID(iNuc));
          
          for jNuc = sameMolecule 
            if Nuclei.Exchangable(jNuc) || iNuc==jNuc
              continue;
            end
            Nuclei.Type{jNuc} = Nuclei.Type{iNuc};
            Nuclei.Spin(jNuc) = Nuclei.Spin(iNuc); % hbar
            Nuclei.StateMultiplicity(jNuc) = Nuclei.StateMultiplicity(iNuc);
            Nuclei.Nuclear_g(jNuc) = Nuclei.Nuclear_g(iNuc);
            Nuclei.Nuclear_g(jNuc) = Nuclei.Nuclear_g(iNuc);
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

Isotopologue = Nuclei;
end