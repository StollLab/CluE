function Isotopologue = newHydronIsotopologue(Nuclei,System)

isSolvent = Nuclei.isSolvent;

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
    elseif strcmp(type,'D') && System.deuterium && ~System.limitToSpinHalf
      Nuclei.Type{iNuc} = '2H';
      Nuclei.Spin(iNuc) = 1; % hbar
      Nuclei.StateMultiplicity(iNuc) = 2*Nuclei.Spin(iNuc) +1;
      Nuclei.Nuclear_g(iNuc) = 0.857438;
      Nuclei.NumberStates(iNuc) = int8(3);
      
      % Set up quadrupole tensors for water deuterons
      if System.nuclear_quadrupole
        
        
        if isempty(Conect)
          error('Nucleus %d is not connected to anything - cannot build NQ tensor.',iNuc);
        end
        for iconnect = Conect
          switch Type{iconnect}
            case {'O','C'}
              zQ = ElectronCenteredCoordinates(iconnect,:) - NuclearCoordinates;
            case {'M','D'}
              xQ = ElectronCenteredCoordinates(iconnect,:) - NuclearCoordinates;
          end
        end
        isOnWater = ~any(iNuc==Indices_nonWater);
        
        if isOnWater
          
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
        Nuclei = setQuadrupoleTensor(e2qQh_,eta_,zQ,xQ,iNuc,Nuclei);
        
      end
    end
    
  end
end


Isotopologue = Nuclei;
end