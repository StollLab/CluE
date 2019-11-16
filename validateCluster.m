function isvalid = validateCluster(Cluster,Nuclei_ValidPair,graphCriterion,clusterSize)

% ENUM
% CONNECTED = 0;  
COMPLETE = 1;

% bool valid for each spin

Adjacency = Nuclei_ValidPair(Cluster,Cluster);

if all(Adjacency==0)
  isvalid = false;
  return;
end

Degree = diag( Adjacency*ones(clusterSize,1) );

Laplacian = Degree - Adjacency;

if graphCriterion == COMPLETE
  % Check if the cluster forms a complete graph.
  
  % Determine the size of the cluster.
  clusterSize = length(Cluster);
  if min(diag(Laplacian)) == clusterSize-1
    isvalid = true;
  else
    isvalid = false;
  end
  return;
end

% Check if the cluster forms a connected graph.

[~,eigenvalues] = eig(Laplacian);
eigenvalues = diag(eigenvalues);

number_of_zero_eigenvalues = sum(abs(eigenvalues) < 1e-12);

if abs(number_of_zero_eigenvalues-1) < 1e-12
  isvalid = true;
elseif number_of_zero_eigenvalues > 1
  isvalid = false;
else
  error('Laplacian matrix has no zero eigenvalues.');  
end

end