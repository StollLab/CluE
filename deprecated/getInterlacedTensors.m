% This funtions adds < Iz = mI| bIzJz | Iz = mI > to the Zeeman term to Jz.   
function intTensors = getInterlacedTensors(tensors,zeroIndex, thisCluster,ClusterComplement,iState,...
  Nuclei_g,state_multiplicity,Nuclei_Coordinates, muN, mu0, hbar)

% Initialize the interlaced externally aware tensors.
intTensors = tensors;

complementSize = length(ClusterComplement);
complementMultiplicity = state_multiplicity(ClusterComplement);

% Get the magnetic quantum number for the external nucleus.
mI = get_mI(complementMultiplicity,iState,complementSize);

% Loop over all nuclei in the primary cluster.
for jj = 1:length(thisCluster)
  jnucleus = thisCluster(jj);
  % Loop over all nuclei in the complement to the primary cluster.
  for ii = 1:length(ClusterComplement)
    inucleus = ClusterComplement(ii);
    %  Get the magnetic coupling between the external and internal nuclei.
    Hdd = constructNuclearDipoleCoupling(Nuclei_g,Nuclei_Coordinates,jnucleus,inucleus, muN, mu0, hbar);
    
    % Adjust the Zeeman term.
    intTensors(:,:,1+jj,1+jj) = intTensors(:,:,1+jj,1+jj) + mI(ii)*Hdd;

  end
end
end

function mI = get_mI(complementMultiplicity,iState,complementSize)

% Initialize output.
mI = zeros(1,complementSize);

% Determine spins.
complementSpin = (complementMultiplicity -1)/2;

% Switch to 0-based counting.
ipsi = (iState-1);

% Loop over external nuclei.
for inuc = 1:complementSize
  
  % Extract multiplicity.
  Omega_ = complementMultiplicity(inuc);
  
  % Get mI index.
  m = mod(ipsi,Omega_) + 1;
  
  % Update ipsi.
  ipsi = floor(ipsi/Omega_);

  % Find mI.  
  mI(inuc) = (m -1 ) - complementSpin(inuc) ;
end

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


