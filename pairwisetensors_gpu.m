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
function [tensors,zeroIndex] = pairwisetensors_gpu(Nuclei_g, Nuclei_Coordinates,...
  Cluster,HF_tensor, magneticField, ge, geff, muB, muN, mu0, hbar,theory)

zeroIndex = min(Cluster) - 1;
Indices = fliplr(Cluster);
N = size(Cluster,2);
tensors = zeros(3,3, N+1,N+1); % nspins by nspins
MethylSelectionRules = [0,0,0;0,0,0;0,0,1];

% ENUM
useEZ       = theory(1);
useNZ       = theory(2);
useHF       = any(theory(3:4));
useNucDD    = any(theory(5:8));

% Electron Zeeman
if useEZ
  tensors(3,3,1,1) = constructElectronZeeman(magneticField,geff, muB, hbar);
end

inucleus = 0;
for i_index_nucleus = Indices
  % inucleus = i_index_nucleus - zeroIndex;
  inucleus  = inucleus + 1;
  % Nuclear Zeeman
  if useNZ
    tensors(3,3,1+inucleus,1+inucleus) = constructNuclearZeeman(Nuclei_g,i_index_nucleus,magneticField, muN, hbar);
  end
  
  % Hyperfine
  if useHF
    hf_tensor = HF_tensor(:,:,i_index_nucleus);
    if any(hf_tensor~=0)
      tensors(:,:,1,1+inucleus ) =  hf_tensor;
    else
    tensors(:,:,1,1+inucleus )= constructHyperfine(Nuclei_g,Nuclei_Coordinates, i_index_nucleus,ge, muB, muN, mu0, hbar);
    end
  end
  
  % Nucleus-Nucleus Coupling
  if ~useNucDD, continue; end
  
  jnucleus = 0;
  for j_index_nucleus = Indices
    % jnucleus = j_index_nucleus - zeroIndex;
    jnucleus = jnucleus + 1;
    if jnucleus <= inucleus, continue; end
    
    % Dipole Coupling
    tensors(:,:,1+inucleus,1+jnucleus) = constructNuclearDipoleCoupling(Nuclei_g,Nuclei_Coordinates,i_index_nucleus,j_index_nucleus, muN, mu0, hbar);
    %       if MethylID(inucleus) > 0 && MethylID(inucleus) == MethylID(jnucleus)
    %         Hamiltonian(:,:,1+inucleus,1+jnucleus) = MethylSelectionRules.*Hamiltonian(:,:,1+inucleus,1+jnucleus);
    %       end
  end
end

end

function electronZeeman = constructElectronZeeman(magneticField,ge, muB, hbar)
electronZeeman = ge*muB*magneticField; % J.
electronZeeman = electronZeeman/(2*pi*hbar); % J -> Hz.
end

function NuclearZeeman = constructNuclearZeeman(Nuclei_g, i_index_nucleus, magneticField, muN, hbar)
gn = Nuclei_g(i_index_nucleus);
NuclearZeeman = -gn*muN*magneticField; % J.
NuclearZeeman = NuclearZeeman/(2*pi*hbar); % J -> Hz
end

function Hyperfine = constructHyperfine(Nuclei_g,Nuclei_Coordinates, i_index_nucleus,ge, muB, muN, mu0, hbar)

gni = Nuclei_g(i_index_nucleus);
r = Nuclei_Coordinates(i_index_nucleus,:)';
n = r/norm(r);
nnt = n*n';
r3 = norm(r)^3;
Hhf = mu0/(4*pi)*ge*muB*gni*muN/r3*(3*nnt - eye(3));
Hyperfine = Hhf/(2*pi*hbar); % Hz. 

return

end

function Hdd = constructNuclearDipoleCoupling(Nuclei_g,Nuclei_Coordinates,i_index_nucleus,j_index_nucleus, muN, mu0, hbar)

gni = Nuclei_g(i_index_nucleus);
gnj = Nuclei_g(j_index_nucleus);
r = Nuclei_Coordinates(i_index_nucleus,:)'-Nuclei_Coordinates(j_index_nucleus,:)';
n = r/norm(r);
nnt = n*n';
r3 = norm(r)^3;
Hdd = -mu0/(4*pi)*gni*gnj*muN^2/r3*(3*nnt - eye(3));
Hdd = Hdd/(2*pi*hbar); % Hz.

end
