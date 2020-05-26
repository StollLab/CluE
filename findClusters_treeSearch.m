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

function Clusters = findClusters_treeSearch(Nuclei,order)

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
    next_vertices = C1(  any( Adjacency(cluster,:)>0,1) );
    if isempty(next_vertices)
      continue;
    end
    
    % Form clusters that contain cluster. 
    C_ = [cluster,uint32(0)] + newVertex.*next_vertices';
    
    % update cluster counters.
    cluster_counter0 = cluster_counter1;
    cluster_counter1 = cluster_counter0 + size(C_,1);
    
    % Reallot memory if necessary. 
    if size(Clusters{clusterSize},1) < cluster_counter1
      Clusters{clusterSize}(end + estNum_ ,:) = zeros(1,clusterSize);
    end  
    
    % Add clusters to cluster list.
    Clusters{clusterSize}(cluster_counter0:cluster_counter1-1,1:clusterSize) = C_;
    
  end
  
  % Remove multicounted clusters.
  Clusters{clusterSize} = reduceClusters(Clusters{clusterSize},true);
  
end



end




function MC = reduceClusters(C,removeZeros)

% Remove zero clusters.
if removeZeros
  C(C(:,1)==0,:) = [];
end

% Sort clusters by spi index.
MS = sort(C,2);

% Remove duplicates.
MC = unique( MS ,'rows');

% Get number of clusters.
numMC = size(MC,1);

% Initialize list of clusters to remove.
removeRows = false(1,numMC);

% Loop through all rows.
for irow = 1:numMC
  
  % Check for repeated indices within the cluster.
  if numel(unique( MC(irow,:)  ))~=numel( MC(irow,:)  )
    
    % Mark cluster for removal.
    removeRows(irow)=true;
  end
end

% Remove clusters that list the same spin multiple times.
MC(removeRows,:) = [];

% Sort cluster lexicographically.
MC = sortrows(MC);

end