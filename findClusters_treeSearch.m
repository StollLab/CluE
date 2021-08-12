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

function Clusters = findClusters_treeSearch(Nuclei,order,adjacencyOrder,...
  inClusters, Method)


% Get number of nuclei
N = Nuclei.number;

inOrder = numel(inClusters);

% Get adjacency matrix and convert to logical
if (inOrder < order) && (inOrder>=2)
  Adjacency = clusters2Adjacency(inClusters{2},N) | logical(Nuclei.Adjacency(:,:,1));
else
  Adjacency = logical(Nuclei.Adjacency(:,:,adjacencyOrder));
  AntiAdjacency = Nuclei.AntiAdjacency(:,:,adjacencyOrder);
end

% Initialize output array
Clusters = cell(order,1);

intDataType = 'uint16';

% Get list of 1-clusters
all1Clusters = cast(1:N,intDataType).';
valid1Clusters = cast(all1Clusters(Nuclei.valid),intDataType);
remove = setdiff(all1Clusters,valid1Clusters);
% Set 1-clusters
if inOrder >= 1
  Clusters{1} = inClusters{1};
else
  Clusters{1} = valid1Clusters;
end


% Get average number of neighbors
q = mean(sum(Adjacency));

% Loop over cluster sizes
for clusterSize = 2:order

  if inOrder >= clusterSize
    Clusters{clusterSize} = inClusters{clusterSize};
    continue;
  end
    
  % Estimate the number of clusters
  estNum_ = ceil(...
    1000 + 1/factorial(clusterSize)*N*q^(clusterSize-1) ...
  );
  
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
    
    neighbors(remove) = false;
    neighbors(cluster) = false;
    neighbors(1:cluster(1)) = false;
    
    neighbornuclei = all1Clusters(neighbors);
    if isempty(neighbornuclei)
      continue
    end
    
    % Remove forbidden clusters  
    for ispin = cluster
      neighbornuclei(AntiAdjacency(ispin, neighbornuclei))=[];
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
  if isempty(C)
    error(['Error in findClusters_treeSearch(): ', ...
      'no cluster found of size ', num2str(clusterSize), ...
      '.  Try relaxing a cutoff.'])
  end
  Clusters{clusterSize} = C(keep,:);
  
end


methylCoupledOnly = any(strcmp(Method.Criteria,'methyl coupled only'));


if methylCoupledOnly
  
  % Methyl Groups
  isMethylCarbon = strcmp(Nuclei.Type,'CH3');
  isMethylHydron = strcmp(Nuclei.Type,'CH3_1H');
  isMethyl = isMethylCarbon | isMethylHydron;
  
  for clusterSize = 1:order
    if ~Method.cutoff.methylCoupledOnly(clusterSize)
      continue;
    end
    
    C = Clusters{clusterSize};
    
    keep  = any(isMethyl(C),2);
    
    Clusters{clusterSize} = C(keep,:);
  end
end

if Method.includeAllSubclusters
  Clusters{order};
  for clusterSize = order-1:-1:1
    
    C = Clusters{clusterSize};
    for ii=1:clusterSize+1
      C = [C; [Clusters{clusterSize+1}(:,1:ii-1), ...
        Clusters{clusterSize+1}(:,ii+1:clusterSize+1) ]  ];
    end
    
    % Sort nuclei in each cluster by nuclei index
    C = sort(C,2);
    
    % Sort clusters
    C = sortrows(C);
    
    % Remove duplicates
    keep = [true; any(C(1:end-1,:)~=C(2:end,:),2)];
    Clusters{clusterSize} = C(keep,:);
    
  end
end

end



function Adjacency = clusters2Adjacency(C2,N1)

Adjacency = zeros(N1);
N2 = size(C2,1);

for ii =1:N2
  m_ = C2(ii,1);
  n_ = C2(ii,2);
  Adjacency(m_,n_) = 1;
  Adjacency(n_,m_) = 1;
end

end


