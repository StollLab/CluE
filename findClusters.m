function Clusters = findClusters(Nuclei,order)

% initialize cluster array
blockSize = 1e5; % number of rows to pre-allocate at a time
Clusters = zeros(blockSize,order);

% loop over all clusters
currCluster = 1:order; % start with first cluster
nClusters = 0;
while ~isempty(currCluster)
  
  % add cluster if it is valid
  if isConnectedCluster(currCluster,Nuclei)
    nClusters = nClusters + 1;
    if nClusters > size(Clusters,1)
      % pre-allocate next block of rows
      Clusters(end+blockSize,:) = 0;
    end
    Clusters(nClusters,:) = currCluster;
  end
  
  % get next cluster
  currCluster = getNextCluster(currCluster,Nuclei.number);
  
end

% trim excess allocation
Clusters = Clusters(1:nClusters,:);

end


function newCluster = getNextCluster(Cluster,numberSpins)

% Initialize the output.
newCluster = Cluster;

% Determine the size of the cluster.
clusterSize = length(Cluster);

% Start at final spin.
for idx = clusterSize:-1:1
  
  % Check if the element can be incremented
  if Cluster(idx) < (numberSpins + idx - clusterSize)
    
    % Increment element.
    newCluster(idx) = newCluster(idx) + 1;
    for jcluster = idx+1:clusterSize
      % Reset the elements after the incremented element.
      newCluster(jcluster) = newCluster(jcluster-1) + 1;
    end
    
    % Return output.
    return
  end
  
end

% There is no next cluster; return the empty set.
newCluster = [];

end


function isvalid = isConnectedCluster(Cluster,Nuclei)

Adjacency = Nuclei.ValidPair(Cluster,Cluster);
Degree = diag(sum(Adjacency,2));
Laplacian = Degree - Adjacency;

if strcmp(Nuclei.graphCriterion,'complete')
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
eigenvalues = eig(Laplacian);
number_of_zero_eigenvalues = sum(abs(eigenvalues) < 1e-12);

if abs(number_of_zero_eigenvalues-1) < 1e-12
  isvalid = true;
elseif number_of_zero_eigenvalues > 1
  isvalid = false;
else
  error('Laplacian matrix has no zero eigenvalues.');
end

end
