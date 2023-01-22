% This function creates a a set of matrices (each member correspond to a 
% particular pairwise coupling) for each nuclei for each electronic state.
% The full Hamiltonian (not used) would be the direct product over all 
% nuclei and electronic states of the sum over all couplings.
% A cluster Hamiltonian, which is used, is the direct product over all
% nuclei in the cluster of the sum over all couplings in the cluster.
function tensors = get_tensors(Nuclei, System,theory,...
  mean_Dipole_z_Z, mean_Dipole_x_iy_Z)



N = Nuclei.number;
Indices = 1:N;
tensors = zeros(3,3, N+1,N+1); % nspins by nspins

% ENUM
useEZ       = theory(1);
useNZ       = theory(2);
useHF       = any(theory(3:4));
useNucDD    = any(theory(5:8));
% useRF       = abs(B1x)>0 || abs(B1y) >0; 

% Electron Zeeman
if useEZ
  tensors(3,3,1,1) = constructElectronZeeman(...
    System.magneticField,System.gMatrix(3,3), System.muB, System.hbar);
end


inucleus = 0;
for i_index_nucleus = Indices
  % inucleus = i_index_nucleus - zeroIndex;
  inucleus  = inucleus + 1;
  % Nuclear Zeeman
  if useNZ
    tensors(3,3,1+inucleus,1+inucleus) = constructNuclearZeeman(...
      Nuclei.Nuclear_g,i_index_nucleus,System.magneticField, System.muN, System.hbar);
  end

  % RF (rotatinge fram approximation)
  if useNZ
    % z
    tensors(3,3,1+inucleus,1+inucleus) = constructNuclearZeemanRotatingFrame(...
      Nuclei.Nuclear_g,i_index_nucleus,System.magneticField, System.muN, ...
      System.hbar,System.RF.nuRF);
    
    % x
    tensors(1,1,1+inucleus,1+inucleus) = constructNuclearZeeman(...
      Nuclei.Nuclear_g,i_index_nucleus,System.RF.B1x, System.muN, System.hbar);
    % y
    tensors(2,2,1+inucleus,1+inucleus) = constructNuclearZeeman(...
      Nuclei.Nuclear_g,i_index_nucleus,System.RF.B1y, System.muN, System.hbar);
  end

  
  % Add mean fields.
  if ~isempty( mean_Dipole_z_Z)
    tensors(3,3,1+inucleus,1+inucleus) = tensors(3,3,1+inucleus,1+inucleus)...
      + mean_Dipole_z_Z(i_index_nucleus);
  end
  if ~isempty(mean_Dipole_x_iy_Z)
    % x
    tensors(1,1,1+inucleus,1+inucleus) = tensors(1,1,1+inucleus,1+inucleus)...
      + real(mean_Dipole_x_iy_Z(i_index_nucleus));
    % y
    tensors(2,2,1+inucleus,1+inucleus) = tensors(2,2,1+inucleus,1+inucleus)...
      + imag(mean_Dipole_x_iy_Z(i_index_nucleus));
  end

  
  % Hyperfine
  if useHF
    hf_tensor = reshape( full( Nuclei.Atensor(i_index_nucleus,:))',3,3) ;
    if any(hf_tensor~=0)
      tensors(:,:,1,1+inucleus ) =  hf_tensor;
    else
    tensors(:,:,1,1+inucleus )= constructHyperfine(...
      Nuclei.Nuclear_g,Nuclei.Coordinates, i_index_nucleus,System.ge, ...
      System.muB, System.muN, System.mu0, System.hbar);
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
    tensors(:,:,1+inucleus,1+jnucleus) = ...
      constructNuclearDipoleCoupling(Nuclei.Nuclear_g,Nuclei.Coordinates,...
      i_index_nucleus,j_index_nucleus, System.muN, System.mu0, System.hbar);
  end
end

end

function electronZeeman = constructElectronZeeman(magneticField,ge, muB, hbar)
electronZeeman = ge*muB*magneticField; % J.
electronZeeman = electronZeeman/(2*pi*hbar); % J -> Hz.
end

function NuclearZeeman = constructNuclearZeeman(Nuclei_g, i_index_nucleus, ...
  magneticField, muN, hbar)
gn = Nuclei_g(i_index_nucleus);
NuclearZeeman = -gn*muN*magneticField; % J.
NuclearZeeman = NuclearZeeman/(2*pi*hbar); % J -> Hz
end

function NuclearZeeman = constructNuclearZeemanRotatingFrame(Nuclei_g, ...
  i_index_nucleus, magneticField, muN, hbar,nuRF)
omega_n = Nuclei_g(i_index_nucleus)*muN*magneticField/hbar;
omegaRF = 2*pi*nuRF;
Omega = omega_n - omegaRF;
NuclearZeeman = -Omega; % J.
NuclearZeeman = NuclearZeeman/(2*pi); % rad/s -> Hz
end

function Hyperfine = constructHyperfine(Nuclei_g,Nuclei_Coordinates, ...
  i_index_nucleus,ge, muB, muN, mu0, hbar)

gni = Nuclei_g(i_index_nucleus);
r = Nuclei_Coordinates(i_index_nucleus,:)';
n = r/norm(r);
nnt = n*n';
r3 = norm(r)^3;
Hhf = mu0/(4*pi)*ge*muB*gni*muN/r3*(3*nnt - eye(3));
Hyperfine = Hhf/(2*pi*hbar); % Hz. 

return

end

function Hdd = constructNuclearDipoleCoupling(Nuclei_g,Nuclei_Coordinates, ...
  i_index_nucleus,j_index_nucleus, muN, mu0, hbar)

gni = Nuclei_g(i_index_nucleus);
gnj = Nuclei_g(j_index_nucleus);
r = Nuclei_Coordinates(i_index_nucleus,:)' ...
  - Nuclei_Coordinates(j_index_nucleus,:)';
n = r/norm(r);
nnt = n*n';
r3 = norm(r)^3;
Hdd = -mu0/(4*pi)*gni*gnj*muN^2/r3*(3*nnt - eye(3));
Hdd = Hdd/(2*pi*hbar); % Hz.

end
