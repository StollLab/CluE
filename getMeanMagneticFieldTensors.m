function newTensors = getMeanMagneticFieldTensors(tensors,spinVector, thisCluster,ClusterComplement,iState,...
  Nuclei_g,Nuclei_Coordinates, muN, mu0, hbar)

% Initialize the interlaced externally aware tensors.
newTensors = tensors;

complementSize = length(ClusterComplement);
clusterSize = length(thisCluster);

% Loop over all nuclei in the primary cluster.
for jj = 1:clusterSize
  jnucleus = thisCluster(jj);
  % Loop over all nuclei in the complement to the primary cluster.
  for ii = 1:complementSize
    inucleus = ClusterComplement(ii);
    %  Get the magnetic coupling between the external and internal nuclei.
    Hdd = constructNuclearDipoleCoupling(Nuclei_g,Nuclei_Coordinates,jnucleus,inucleus, muN, mu0, hbar);
    
    % Adjust the Zeeman term.
    newTensors(:,:,1+jj,1+jj) = newTensors(:,:,1+jj,1+jj) + diag(spinVector(:,iState, inucleus)'*Hdd);

  end
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