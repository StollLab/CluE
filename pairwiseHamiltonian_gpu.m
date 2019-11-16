% This function creates a a set of matrices (each member correspond to a particular pairwise coupling) for each nuclei for each electronic state.
% The full Hamiltonian (not used) would be the direct product over all nuclei and
% electronic states of the sum over all couplings.  
% A cluster Hamiltonian, which is used, is the direct product over all
% nuclei in the cluster of the sum over all couplings in the cluster.

% CHANGES
% Hamiltonian -> zeros(3,3, N+1,N+1)
% Nuclei -> arrays Nuclei_g, Nuclei_Coordinates
% System -> magneticField, ge, muB, muN, mu0, hbar
%{
  Nuclei_Coordinates = Nuclei.Coordinates;
  Nuclei_g = Nuclei.Nuclear_g;
  ge=System.gMatrix;
  magneticField = System.magneticField;
  muB = System.muB;
  muN = System.muN;
  mu0 = System.mu0;
  hbar = System.hbar;
%}
function [Hamiltonian,zeroIndex] = pairwiseHamiltonian_gpu(Nuclei_g, Nuclei_Coordinates,Cluster,magneticField, ge, muB, muN, mu0, hbar,useHamiltonian,MethylID)

  zeroIndex = min(Cluster) - 1;
  Indices = fliplr(Cluster);
  N = size(Cluster,2);
  Hamiltonian = zeros(3,3, N+1,N+1); % nspins by nspins
  inucleus = 0;
  MethylSelectionRules = [0,0,0;0,0,0;0,0,1];
  
  % ENUM
  EZEEMAN = 1; NZEEMAN = 2; HYPERFINE = 3; NDIPOLE = 4;
  
  % Electron Zeeman
  if useHamiltonian(EZEEMAN)
    Hamiltonian(:,:,1,1) =  constructElectronZeeman(magneticField,ge, muB, hbar);
  end
  for i_index_nucleus = Indices
    %     inucleus = i_index_nucleus - zeroIndex;
    inucleus  = inucleus + 1;
    % Nuclear Zeeman
    if useHamiltonian(NZEEMAN)
      Hamiltonian(:,:,1+inucleus,1+inucleus) =  constructNuclearZeeman(Nuclei_g,i_index_nucleus,magneticField, muN, hbar);
    end
    
    % Hyperfine
    if useHamiltonian(HYPERFINE)
      Hamiltonian(:,:,1,1+inucleus )= constructHyperfine(Nuclei_g,Nuclei_Coordinates, i_index_nucleus,ge, muB, muN, mu0, hbar);
    end
    
    % Nucleus-Nucleus Coupling
    if ~useHamiltonian(NDIPOLE)
      continue;
    end
    
    jnucleus = 0;
    for j_index_nucleus = Indices
      % jnucleus = j_index_nucleus - zeroIndex;
      jnucleus = jnucleus + 1;
      if jnucleus <= inucleus
        continue;
      end
 
      % Dipole Coupling
      Hamiltonian(:,:,1+inucleus,1+jnucleus) = constructNuclearDipoleCoupling(Nuclei_g,Nuclei_Coordinates,i_index_nucleus,j_index_nucleus, muN, mu0, hbar);      
%       if MethylID(inucleus) > 0 && MethylID(inucleus) == MethylID(jnucleus)
%         Hamiltonian(:,:,1+inucleus,1+jnucleus) = MethylSelectionRules.*Hamiltonian(:,:,1+inucleus,1+jnucleus);
%       end
    end
  end
   
end

function electronZeeman = constructElectronZeeman(magneticField,ge, muB, hbar) 
  electronZeeman = ge*muB*magneticField; % J.
  electronZeeman = electronZeeman/(2*pi*hbar); % Hz. 
end

function NuclearZeeman = constructNuclearZeeman(Nuclei_g, i_index_nucleus, magneticField, muN, hbar)
   
    NuclearZeeman = Nuclei_g(i_index_nucleus)*magneticField*muN*eye(3); % J.
    NuclearZeeman = NuclearZeeman/(2*pi*hbar); % Hz
  
end

function Hyperfine = constructHyperfine(Nuclei_g,Nuclei_Coordinates, i_index_nucleus,ge, muB, muN, mu0, hbar)

  gni = Nuclei_g(i_index_nucleus);
  r = Nuclei_Coordinates(i_index_nucleus,:);
  if size(r,2)==3
    r=r';
  end
  n = r/norm(r);
  nnt = n*n'; 
  r3 = norm(r)^3;
  Hdd = mu0/(4*pi)*ge*muB*gni*muN/r3*(eye(3)-3*nnt);
  Hyperfine = -Hdd/(2*pi*hbar); % Hz.

  return;

end

function Hdd = constructNuclearDipoleCoupling(Nuclei_g,Nuclei_Coordinates,i_index_nucleus,j_index_nucleus, muN, mu0, hbar)
   
  gni = Nuclei_g(i_index_nucleus);
  gnj = Nuclei_g(j_index_nucleus);
  r = Nuclei_Coordinates(i_index_nucleus,:)-Nuclei_Coordinates(j_index_nucleus,:);
  if size(r,2)==3
    r=r';
  end
  n = r/norm(r);
  nnt = n*n'; 
  r3 = norm(r)^3;
  Hdd = -mu0/(4*pi)*gni*gnj*muN^2/r3*(eye(3)-3*nnt);
  Hdd = Hdd/(2*pi*hbar); % Hz.
  
end
