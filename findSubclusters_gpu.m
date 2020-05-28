% ========================================================================
% Subcluster Function gpu
% ========================================================================

function Indices = findSubclusters_gpu(Clusters,clusterSize,iCluster,reference_clusterSize)
% Indices = Indices{subCluster_size} = SubclusterIndices_clusterSize
% SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
% the jth cluster of size subCluster_size that is a subcluster of
% the ith cluster of size clusterSize.
% In a valid cluster, all indices must be natural numbers.

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
if clusterSize > 2
  possibleSubClusters_2 = zeros(numberSubClusters(2),clusterSize);
  
  if clusterSize > 3
    possibleSubClusters_3 = zeros(numberSubClusters(3),clusterSize);
    
    if clusterSize > 4
      possibleSubClusters_4 = zeros(numberSubClusters(4),clusterSize);
      
      if clusterSize > 5
        possibleSubClusters_5 = zeros(numberSubClusters(5),clusterSize);
        
        if clusterSize > 6
          possibleSubClusters_6 = zeros(numberSubClusters(5),clusterSize);
          
          if clusterSize > 7
            possibleSubClusters_7 = zeros(numberSubClusters(5),clusterSize);
            
            if clusterSize > 8
              possibleSubClusters_8 = zeros(numberSubClusters(5),clusterSize);
              
              if clusterSize > 9
                possibleSubClusters_9 = zeros(numberSubClusters(5),clusterSize);
                
                if clusterSize > 10
                  possibleSubClusters_10 = zeros(numberSubClusters(5),clusterSize);
                  
                  if clusterSize > 11
                    possibleSubClusters_11 = zeros(numberSubClusters(5),clusterSize);
                    
                    if clusterSize > 12
                      error('Cluster size not supported.');
                    end; end; end; end; end; end; end; end; end; end; end

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
  switch sc_size
    
    case 2
      possibleSubClusters_2(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 3
      possibleSubClusters_3(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 4
      possibleSubClusters_4(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 5
      possibleSubClusters_5(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 6
      possibleSubClusters_6(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 7
      possibleSubClusters_7(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 8
      possibleSubClusters_8(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 9
      possibleSubClusters_9(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 10
      possibleSubClusters_10(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 11
      possibleSubClusters_11(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
  end
  
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
  
  % Get list of possible subclusters.
  switch isize
    case 2
      possibleSubClusters_ = possibleSubClusters_2;
    case 3
      possibleSubClusters_ = possibleSubClusters_3;
    case 4  
      possibleSubClusters_ = possibleSubClusters_4;
    case 5
      possibleSubClusters_ = possibleSubClusters_5;
    case 6  
      possibleSubClusters_ = possibleSubClusters_6;
    case 7  
      possibleSubClusters_ = possibleSubClusters_7;
    case 8  
      possibleSubClusters_ = possibleSubClusters_8;
    case 9  
      possibleSubClusters_ = possibleSubClusters_9;
    case 10  
      possibleSubClusters_ = possibleSubClusters_10;
    case 11  
      possibleSubClusters_ = possibleSubClusters_11;
  end

  % Loop through all subclusters.
  for jCluster = 1:numberSubClusters(isize)
    
    % Get the jth subcluster.
    SubCluster = Cluster(possibleSubClusters_(jCluster,:)==1);
    
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
