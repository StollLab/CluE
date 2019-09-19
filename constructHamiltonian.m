% This function creates a a set of matrices (each member correspond to a particular pairwise coupling) for each nuclei for each electronic state.
% The full Hamiltonian (not used) would be the direct product over all nuclei and
% electronic states of the sum over all couplings.  
% A cluster Hamiltonian, which is used, is the direct product over all
% nuclei in the cluster of the sum over all couplings in the cluster.

function [H_diagonal, H_offDiagonal,zeroIndex] = constructHamiltonian(System,Nuclei,Cluster)
  zeroIndex = min(Cluster) - 1;
  maxNuclearIndex = max(Cluster) - zeroIndex;
  Indices = fliplr(Cluster); 
  H_offDiagonal = [];

  for eState = (2*System.Electron.spin+1):-1:1
    % Electron Zeeman
    ms = eState - 1 - System.Electron.spin; 
    electronZeeman = constructElectronZeeman(System,eState);
    H_diagonal{eState, maxNuclearIndex+1, maxNuclearIndex+1} = electronZeeman;
    for i_index_nucleus = Indices
      inucleus = i_index_nucleus - zeroIndex;
      % Nuclear Zeeman
      H_diagonal{eState, inucleus, inucleus} =  constructNuclearZeeman(System, Nuclei,i_index_nucleus);
      
      % Hyperfine
      H_diagonal{eState, inucleus, inucleus} = H_diagonal{eState, inucleus, inucleus} + constructHyperfine(System, ms, Nuclei, i_index_nucleus);
      
      % Nucleus-Nucleus Coupling
    
      for j_index_nucleus = Indices
        jnucleus = j_index_nucleus - zeroIndex;
        if jnucleus == inucleus
          continue;
        end
 
        % Dipole Coupling
        H_diagonal{eState, inucleus, jnucleus} = constructNuclearBath(System, Nuclei, i_index_nucleus, j_index_nucleus);
        H_offDiagonal{eState, inucleus, jnucleus} = constructNuclearBathFlipFlop(System, Nuclei, i_index_nucleus, j_index_nucleus);
        
%         if max(isnan(H_diagonal{eState, inucleus, jnucleus})) || max(isinf(H_diagonal{eState, inucleus, jnucleus})) || max(isnan(H_offDiagonal{eState, inucleus, jnucleus})) || max(isinf(H_offDiagonal{eState, inucleus, jnucleus}))
%           error('Hamiltonian could not be calculated.');
%         end
      end
    end
  end   
end

function electronZeeman = constructElectronZeeman(System,eState)
  % eState counts the electron state from lowest to highest.
  % Spin 1/2: ms = -1/2 --> eState = 1, ms = +1/2 --> eStat e= 2.
  % Spin 1: ms = -1 --> eState=1, ms = 0 --> eState = 2, ms = +1 --> eState = 3 .
  ms = eState - 1 - System.Electron.spin; 
  ge = System.Electron.g;
  electronZeeman = ms*ge*System.muB*System.magneticField; % J.
  electronZeeman =electronZeeman/(2*pi*System.hbar); % Hz. 
end

function NuclearZeeman = constructNuclearZeeman(System, Nuclei, i_index_nucleus)
  NuclearZeeman = ones(2*Nuclei.Spin(i_index_nucleus) + 1,1);
  for ispinState = 1:(2*Nuclei.Spin(i_index_nucleus) + 1)
    %mI = ispinState - -1 -Nuclei.Spin(i_index_nucleus); 
    mI = ispinState - 1 -Nuclei.Spin(i_index_nucleus); 
    NuclearZeeman(ispinState) = -Nuclei.Nuclear_g(i_index_nucleus)*System.magneticField*System.muN*mI; % J.
    NuclearZeeman(ispinState) = NuclearZeeman(ispinState)/(2*pi*System.hbar);
  end
end

function Hyperfine = constructHyperfine(System, ms, Nuclei, i_index_nucleus)

  cosTheta2 = cosFieldAngle(System.Electron.Coordinates,Nuclei.Coordinates(i_index_nucleus,:));
  cosTheta2 = cosTheta2*cosTheta2;
  r = norm(System.Electron.Coordinates - Nuclei.Coordinates(i_index_nucleus,:));
  r3 = r^3;
  km = System.mu0/(4*pi);
  ge = System.Electron.g;
  muB = System.muB;
  gn = Nuclei.Nuclear_g(i_index_nucleus);
  muN = System.muN;
  a = km*ge*muB*gn*muN/r3*(3*cosTheta2-1); % J.
  a = a/(2*pi*System.hbar); % Hz.;
  a = a + Nuclei.Hyperfine(i_index_nucleus);
  
  Hyperfine = ones(2*Nuclei.Spin(i_index_nucleus) + 1,1);
  for ispinState = 1:(2*Nuclei.Spin(i_index_nucleus) + 1)
    mI = ispinState - 1 -Nuclei.Spin(i_index_nucleus);
    Hyperfine(ispinState) = a*ms*mI; % Hz.
  end
end

function Bath = constructNuclearBath(System, Nuclei, i_index_nucleus, j_index_nucleus)
  % The z-axis is defined as the direction of the applied magnetic field.
 if Nuclei.Index(i_index_nucleus) == Nuclei.Index(j_index_nucleus)
   Bath = []; 
   return;
 elseif Nuclei.Index(i_index_nucleus) < Nuclei.Index(j_index_nucleus)
   
   b = 0;
   if abs(Nuclei.Nuclear_g(i_index_nucleus)-Nuclei.Nuclear_g(j_index_nucleus))<1e-6
     
     cosTheta2 = cosFieldAngle(Nuclei.Coordinates(i_index_nucleus,:),Nuclei.Coordinates(j_index_nucleus,:));
     cosTheta2 = cosTheta2*cosTheta2;
     
     km = System.mu0/(4*pi);
     muN = System.muN;
     gni = Nuclei.Nuclear_g(i_index_nucleus);
     gnj = Nuclei.Nuclear_g(j_index_nucleus);
     r = norm(Nuclei.Coordinates(i_index_nucleus,:)-Nuclei.Coordinates(j_index_nucleus,:));
     r3 = r^3;
     b = km*gni*muN*gnj*muN/r3*(1-3*cosTheta2); % J.
     b = b/(2*pi*System.hbar); % Hz.
   end
   Bath = b*spinZ(Nuclei.Spin(i_index_nucleus));
 elseif Nuclei.Index(i_index_nucleus) > Nuclei.Index(j_index_nucleus)
   Bath = spinZ(Nuclei.Spin(i_index_nucleus));
 end
 Bath = diag(Bath);
end

function Bath = constructNuclearBathFlipFlop(System,Nuclei, i_index_nucleus, j_index_nucleus)
 % The z-axis is defined as the direction of the applied magnetic field.
 if Nuclei.Index(i_index_nucleus) == Nuclei.Index(j_index_nucleus)
   Bath = []; 
   return;
 elseif Nuclei.Index(i_index_nucleus) < Nuclei.Index(j_index_nucleus)
   
   b = 0;
   if abs(Nuclei.Nuclear_g(i_index_nucleus)-Nuclei.Nuclear_g(j_index_nucleus))<1e-6
     
     cosTheta2 = cosFieldAngle(Nuclei.Coordinates(i_index_nucleus,:),Nuclei.Coordinates(j_index_nucleus,:));
     cosTheta2 = cosTheta2*cosTheta2;
     
     km = System.mu0/(4*pi);
     muN = System.muN;
     gni = Nuclei.Nuclear_g(i_index_nucleus);
     gnj = Nuclei.Nuclear_g(j_index_nucleus);
     r = norm(Nuclei.Coordinates(i_index_nucleus,:)-Nuclei.Coordinates(j_index_nucleus,:));
     r3 = r^3;
     b = -0.25*km*gni*muN*gnj*muN*(1-3*cosTheta2)/r3; % J.
     b = b/(2*pi*System.hbar); % Hz.
   end
   Bath = b*spinRaise(Nuclei.Spin(i_index_nucleus));
   Bath = [Bath(1:end-1,2:end)];
   Bath = diag(Bath);
   Bath(end+1) = +1;
 elseif Nuclei.Index(i_index_nucleus) > Nuclei.Index(j_index_nucleus)
   Bath = spinLower(Nuclei.Spin(i_index_nucleus));
   Bath = [Bath(2:end,1:end-1)];
   Bath = diag(Bath);
   Bath(end+1) = -1;
 end
end

function Iz = spinZ(spin)
Iz = eye(2*spin+1);
for ii = [1:(2*spin+1)]
  Iz(ii,ii) = spin+ 1 - ii;
end
end

function Iminus = spinLower(spin)
Iminus = zeros(2*spin+1);
for ii = [1:(2*spin+1)-1]
  Iminus(ii+1,ii) = sqrt(spin*(spin+1)-(spin - ii)*(spin+ 1 - ii));
end
end

function Iplus = spinRaise(spin)
Iplus = zeros(2*spin+1);
for ii = [1:(2*spin+1)-1]
  Iplus(ii,ii+1) = sqrt(spin*(spin+1)-(spin - ii)*(spin+ 1 - ii));
end
end

function cos_theta = cosFieldAngle(ParticleCoordinates_1,ParticleCoordinates_2)
   R = ParticleCoordinates_2 - ParticleCoordinates_1;
   cos_theta = R*[0;0;1]/norm(R); 
end
