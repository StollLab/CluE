function isvalid = validateCluster(Cluster,Nuclei_ValidPair,completeGraph)

clusterSize  = numel(Cluster);
if clusterSize == 1
  isvalid = true;
  return;
end

Adjacency = Nuclei_ValidPair(Cluster,Cluster,clusterSize);

if all(Adjacency(:)==0)
  isvalid = false;
  return
end

Degree = diag(sum(Adjacency,2));
Laplacian = Degree - Adjacency;

if completeGraph
  % Check if the cluster forms a complete graph.
  
  clusterSize = length(Cluster);
  isvalid = min(diag(Laplacian)) == clusterSize-1;

else

  % Check if the cluster forms a connected graph.
  eigenvalues = eig(Laplacian);

  nZeroEigenvalues = sum(abs(eigenvalues) < 1e-12);

  if abs(nZeroEigenvalues-1) < 1e-12
    isvalid = true;
  elseif nZeroEigenvalues > 1
    isvalid = false;
  else
    error('Laplacian matrix has no zero eigenvalues.');  
  end

end

end
