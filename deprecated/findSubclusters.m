% ========================================================================
% Subcluster Function gpu
% ========================================================================

function Indices = findSubclusters(Clusters,clusterSize,iCluster,reference_clusterSize)
% Indices = Indices{subCluster_size} = SubclusterIndices_clusterSize
% SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
% the jth cluster of size subCluster_size that is a subcluster of
% the ith cluster of size clusterSize.
% In a valid cluster, all indices must be natural numbers.

% SubclusterIndices{clustersize}(:,subclusterSize,iCluster)

% Initialize Indices assuming that all n choose k subsets are subclusters. 
Indices = zeros(NchooseK(reference_clusterSize, ceil(reference_clusterSize/2)),reference_clusterSize);

% Each cluster is a subset of itself.
Indices(1,clusterSize)=iCluster;

% Vector of number of subsets for each subset size.
numberSubClusters = NchooseK(clusterSize,1:clusterSize);

% Check all non-empty subsets have been found.
if clusterSize == 1
    return;
end

% Initialize lists of all possible subclusters.
for iSize = 1:clusterSize
  possibleSubClusters{iSize} = zeros(numberSubClusters(iSize),clusterSize);
end

% Initialize subcluster counters.
subcluster_count = zeros(1,clusterSize);

% Loop over all proper subsets. 
for isc = 1:2^(clusterSize)-1
  
  % Get the binary representation of the subset.
  subcluster_str = dec2bin(isc);
  
  % Determine the number of filler zeros required.
  nZeros = clusterSize - length(subcluster_str);
  
  % Generate str representing the subcluster.
  subcluster_str = [repmat('0',1,nZeros) subcluster_str];
  
  % Determine the size of the subcluster.
  sc_size = sum(subcluster_str=='1');
  
  % Increment subcluster counter.
  subcluster_count(sc_size) = subcluster_count(sc_size) + 1;
  
  % Convert subcluster_str to logical array. 
  possibleSubClusters{sc_size}(subcluster_count(sc_size), ...
    1:clusterSize) = subcluster_str=='1';
   
  
end

% Get the ith cluster of size clusterSize from the list of clusters.
Cluster = Clusters(iCluster, 1:clusterSize ,clusterSize);

% Loop through the elements of the cluster.
for ii =1:length(Cluster)
  % Get spin ii.
  inucleus = Cluster(ii);
  
  % Find the index in Clusters that cooresponds to inucleus 
  % and store it in Indices. 
  Indices(ii,1) = find(Clusters(:,1,1)'==inucleus);
end

% Check if all subcluster sizes were accounted for.
if clusterSize == 2
  return;
end
%--------------------------------------------------------------------------

for isize = 2:clusterSize-1  
  
  % Loop through all subclusters.
  for jCluster = 1:numberSubClusters(isize)
    
    % Get the jth subcluster.
    SubCluster = Cluster(possibleSubClusters{isize}(jCluster,:)==1);
    
    % Search for a valid cluster that equals the subcluster.
    Search = Clusters(:, 1:isize, isize)==SubCluster;
    
    % Locate where the subcluster is.
    subclusterIndex = find( all(Search,2));
    
    % Set invalid subcluster indices to zero.
    if isempty(subclusterIndex)
      subclusterIndex = 0;
    end
    
    % Record the subcluster's index.
    Indices(jCluster,isize) = subclusterIndex;
  end
end

end
