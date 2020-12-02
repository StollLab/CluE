function Isotopologue = newHydronIsotopologue(Nuclei,System)
isSolvent = Nuclei.isSolvent;

switch System.HydrogenExchange
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