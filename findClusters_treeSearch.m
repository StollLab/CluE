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

function Clusters = findClusters_treeSearch2(Nuclei,order)

% Initialize clusters.
Clusters = cell(order,1);

% Get number of nuclei.
N = Nuclei.number;

% Get list of nuclear spins.
C1 = [1:N];

% Set 1-clusters.
Clusters{1} = C1';

%Get adjacency matrix.
Adjacency = Nuclei.ValidPair(:,:,1);

% Get average number of neighbors.
q = mean(sum(Adjacency));

% Loop over cluster sizes.
for clusterSize = 2:order
  
  % Estimate the number of clusters.
  estNum_ = ceil(1/factorial(clusterSize)*double(N)*q^(clusterSize-1));
  
  % Initialize cluster list.
  Clusters{clusterSize} = uint32(zeros(estNum_,clusterSize));
  
  % Get number of clusters for clusterSize -1.
  numC_ = size(Clusters{clusterSize-1},1);
  
  % Initialize template for adding new vertices.
  newVertex = uint32(zeros(1,clusterSize));
  newVertex(clusterSize)  =1;
  
  % Initialize cluster counters.
  cluster_counter0 = 1;
  cluster_counter1 = 1;
  
  % Loop over all clusters of size clusterSize -1.
  for icluster = 1:numC_
    cluster = Clusters{clusterSize-1}(icluster,:);
    
    % Find all vertices connected to spins in cluster.
    neighbors = any(Adjacency(cluster,:),1);
    
    
    neighbors(cluster) = false;
    neighbors(1:cluster(1)) = false;

    neighborspins = C1(neighbors);
    if isempty(neighborspins)
      continue;
    end
    
    % Form clusters that contain cluster. 
    newclusters = [repmat(cluster,numel(neighborspins),1), neighborspins(:)];
    
    % update cluster counters.
    cluster_counter0 = cluster_counter1;
    cluster_counter1 = cluster_counter0 + size(newclusters,1);
    
    % Reallot memory if necessary. 
    if size(Clusters{clusterSize},1) < cluster_counter1
      Clusters{clusterSize}(end + estNum_ ,:) = zeros(1,clusterSize);
    end  
    
    % Add clusters to cluster list.
    Clusters{clusterSize}(cluster_counter0:cluster_counter1-1,1:clusterSize) = newclusters;
    
  end
  
  % Remove multicounted clusters.

  C = Clusters{clusterSize};
  % Remove zero clusters.
  C(C(:,1)==0,:) = [];

  
  % Sort clusters by spin index.
  MS = sort(C,2);
  
  % Remove duplicates.
  MC = unique( MS ,'rows');
  
  Clusters{clusterSize} = sortrows(MC);
  
end



end


