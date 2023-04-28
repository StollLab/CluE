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

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function Clusters = findClusters_treeSearch(Nuclei,order,adjacencyOrder,...
  inClusters, Method)


% Get number of nuclei
N = Nuclei.number;

inOrder = numel(inClusters);

% Get adjacency matrix and convert to logical
if (inOrder < order) && (inOrder>=2)
  Adjacency ...
  = clusters2Adjacency(inClusters{2},N) | logical(Nuclei.Adjacency(:,:,1));
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

  % 
  useCompleteGraph = strcmp(Method.graphCriterion,'complete');
  if useCompleteGraph && clusterSize > 2
    for iCluster = 1:size(C,1)
      if ~keep(iCluster), continue; end
      
      Cluster = C(iCluster,:);
      keep(iCluster) = validateCluster(Cluster,Adjacency,useCompleteGraph);

    end
  end

  if ~isempty(C)
    
    Clusters{clusterSize} = C(keep,:);

  elseif Method.emptyClusterSetsOkay
  
    Clusters{clusterSize} = C;

  else
    error(['Error in findClusters_treeSearch(): ', ...
      'no cluster found of size ', num2str(clusterSize), ...
      '.  Try relaxing a cutoff or setting',...
      '\n  Method.emptyClusterSetsOkay = true;'])
  end

end


%doMethylCoupledOnly = any(strcmp(Method.Criteria,'methyl coupled only'));
doMethylCoupledOnly = any(Method.neighborCutoff.methylCoupledOnly);


if doMethylCoupledOnly
  
  % Methyl Groups
  isMethylCarbon = strcmp(Nuclei.Type,'CH3');
  isMethylHydron = strcmp(Nuclei.Type,'CH3_1H');
  isMethyl = isMethylCarbon | isMethylHydron;

  for clusterSize = 3:order
    if ~Method.neighborCutoff.methylCoupledOnly(clusterSize)
      continue;
    end
    
    C = Clusters{clusterSize};
    
    keep  = sum(isMethyl(C),2) >= Method.neighborCutoff.methylCoupledOnlyNumber;
    for ii=1:length(keep)
      if ~keep(ii), continue; end
      cluster = C(ii,:);
      IDs = Nuclei.MethylID(cluster);
      IDs = IDs(IDs>0);
      uniqueIDs = unique(IDs);
      keep(ii) = all(sum(uniqueIDs == IDs',1)==3);
    end


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


if Method.useMethylPseudoParticles

  for cluster_size = 1:order
    % Initialize keep seletor.
    keep = true( size(Clusters{cluster_size},1),1 );

    % Select all clusters with at least one methyl hydron.
    if cluster_size == 1
      sele = Nuclei.MethylID(Clusters{cluster_size})' > 0;
      keep(sele) = false;
    else
      sele = any(Nuclei.MethylID(Clusters{cluster_size}) > 0,2);
      % Remove clusters that do not contain all 3 hydron of a methyl
      keep(sele) = remove_incomplete_methyls(...
        Nuclei.MethylID(Clusters{cluster_size}(sele,:)));
    end

    % Finalize removal.
    Clusters{cluster_size} = Clusters{cluster_size}(keep,:);
  end
end


if Method.randomize_cluster_ordering
  for cluster_size = 1:order

    n_cluster = size(Clusters{cluster_size},1);
    new_ordering = randperm(n_cluster);
    Clusters{cluster_size} = Clusters{cluster_size}(new_ordering,:);

  end
end

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% cluster_methylIDs is an array of methyl ids, where each row is mapped
% from a cluster to the methyl ids of the corresponding atom.
function keep = remove_incomplete_methyls(methylIDs)

% Get the number of cluster.
n_clusters = size(methylIDs,1);

% 1 and 2-clusters cannot hold a full methyl group.
if size(methylIDs,2)<3
  keep = false(n_clusters,1);
  return;
end

% Initialize selector.
keep = true(n_clusters,1);

% Loop through clusters.
for icluster = 1:n_clusters
  % Find unique ids.
  unique_ids = unique( methylIDs(icluster,:) );

  % Remove clusters that do not contain all or none of the protons from every 
  % methyl group.
  keep(icluster) = all( ...
    sum(methylIDs(icluster,:) == unique_ids',2)==3 ...
    | unique_ids' ==0);
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function isvalid = validateCluster(Cluster,Nuclei_ValidPair,completeGraph)

clusterSize  = numel(Cluster);
if clusterSize == 1
  isvalid = true;
  return;
end

if size(Nuclei_ValidPair,3) >= clusterSize
  Adjacency = Nuclei_ValidPair(Cluster,Cluster,clusterSize);
else
  Adjacency = Nuclei_ValidPair(Cluster,Cluster,end);
end

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
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>