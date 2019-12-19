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

if ~isfield(Nuclei,'graphCriterion')
  Nuclei.graphCriterion = 'connected';
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
  if validateCluster(currCluster,Nuclei.ValidPair,completeGraph)
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
