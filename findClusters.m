% findClusters   Generate list of cluster of a given order
%
%  Clusters = findClusters(Nuclei,order)
%
% Inputs:
%   Nuclei           ... structure with fields
%    .number         ... number of nuclei
%    .ValidPair      ... adjacency matrix
%    .graphCriterion ... method for establishing cluster connectivity
%                        'complete', ''
%   order            ... ckuster order: 1, 2, 3, etc
%
% Outputs:
%   Clusters         ... Mxorder array of clusters of size order, one
%                        custer per row

function Clusters = findClusters(Nuclei,order)

if Nuclei.number<order
  error('Cluster order cannot be larger than number of nuclear spins.');
end

completeGraph = strcmp(Nuclei.graphCriterion,'complete');
nSpins = Nuclei.number;

% initialize cluster array
blockSize = 1e5; % number of rows to pre-allocate at a time
Clusters = zeros(blockSize,order);

% loop over all clusters
currCluster = 1:order; % start with first cluster
nClusters = 0;
while ~isempty(currCluster)
  
  % add cluster if it is valid
  if isConnectedCluster(currCluster,Nuclei.ValidPair,completeGraph)
    nClusters = nClusters + 1;
    if nClusters > size(Clusters,1)
      % pre-allocate next block of rows
      Clusters(end+blockSize,:) = 0;
    end
    Clusters(nClusters,:) = currCluster;
  end
  
  % get next cluster
  currCluster = getNextCluster(currCluster,nSpins);
  
end

% trim excess allocation
Clusters = Clusters(1:nClusters,:);

end


function nextCluster = getNextCluster(Cluster,nSpins)

% Initialize the output.
nextCluster = Cluster;

% Determine the size of the cluster.
clusterSize = length(Cluster);

% Start at final spin.
for idx = clusterSize:-1:1
  
  % Check if the element can be incremented
  if Cluster(idx) < nSpins + idx - clusterSize
    
    % Increment element.
    nextCluster(idx) = nextCluster(idx) + 1;
    for jcluster = idx+1:clusterSize
      % Reset the elements after the incremented element.
      nextCluster(jcluster) = nextCluster(jcluster-1) + 1;
    end
    
    % Return output.
    return
  end
  
end

% There is no next cluster; return the empty set.
nextCluster = [];

end


function isvalid = isConnectedCluster(Cluster,Adjacency,completeGraph)

Adjacency = Adjacency(Cluster,Cluster);
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
