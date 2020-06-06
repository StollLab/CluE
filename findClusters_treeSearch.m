% findClusters   Generate list of cluster of a given order
%
%   Clusters = findClusters(Nuclei,order)
%
% Inputs:
%   Nuclei       ... structure with fields
%    .ValidPair  ... adjacency matrix
%   order        ... ckuster order: 1, 2, 3, etc
%
% Outputs:
%   Clusters     ... Mxorder array of clusters of size order, one
%                    cluster per row

function Clusters = findClusters_treeSearch(Nuclei,order)

% Get adjacency matrix and convert to logical
Adjacency = logical(Nuclei.ValidPair(:,:,1));
% Get number of nuclei
N = size(Adjacency,1);

% Initialize output array
Clusters = cell(order,1);

intDataType = 'uint16';

% Get list of 1-clusters
C1 = cast(1:N,intDataType).';

% Set 1-clusters
Clusters{1} = C1;

% Get average number of neighbors
q = mean(sum(Adjacency));

% Loop over cluster sizes
for clusterSize = 2:order
  
  % Estimate the number of clusters
  estNum_ = ceil(1/factorial(clusterSize)*N*q^(clusterSize-1));
  
  % Initialize cluster list
  Clusters{clusterSize} = zeros(estNum_,clusterSize,intDataType);
  
  % Get number of clusters for clusterSize-1
  numC_ = size(Clusters{clusterSize-1},1);
  
  % Initialize template for adding new vertices
  newVertex = zeros(1,clusterSize,intDataType);
  newVertex(clusterSize) = 1;
  
  % Initialize cluster counters
  cluster_counter0 = 1;
  cluster_counter1 = 1;
  
  % Loop over all clusters of size clusterSize-1
  for icluster = 1:numC_
    cluster = Clusters{clusterSize-1}(icluster,:);
    
    % Find all vertices connected to nuclei in cluster
    neighbors = any(Adjacency(:,cluster),2);
    
    neighbors(cluster) = false;
    neighbors(1:cluster(1)) = false;
    
    neighbornuclei = C1(neighbors);
    if isempty(neighbornuclei)
      continue
    end
    
    % Form expanded clusters that contain parent cluster
    newclusters = [repmat(cluster,numel(neighbornuclei),1), neighbornuclei];
    
    % Update cluster counters
    cluster_counter0 = cluster_counter1;
    cluster_counter1 = cluster_counter0 + size(newclusters,1);
    
    % Reallot memory if necessary
    if size(Clusters{clusterSize},1) < cluster_counter1
      Clusters{clusterSize}(end+estNum_,clusterSize) = 0;
    end  
    
    % Add clusters to cluster list
    Clusters{clusterSize}(cluster_counter0:cluster_counter1-1,:) = newclusters;
    
  end
  
  % Trim array
  C = Clusters{clusterSize}(1:cluster_counter1-1,:);
  
  % Sort nuclei in each cluster by nuclei index
  C = sort(C,2);
  
  % Sort clusters
  C = sortrows(C);
  
  % Remove duplicates
  keep = [true; any(C(1:end-1,:)~=C(2:end,:),2)];
  Clusters{clusterSize} = C(keep,:);
  
end

end
