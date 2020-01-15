% ========================================================================
% Subcluster Function gpu
% ========================================================================

function Indices = findSubclusters_gpu(Clusters,clusterSize,iCluster,reference_clusterSize)
% from Indices{subCluster_size} = list of all jCluster such that Clusters{subCluster_size}(jCluster,:) is a subcluster of Clusters{clusterSize}(iCluster,:)
% to
% Given the ith cluster of size clusterSize, 
% Indices(jCluster,subCluster_size) = jth cluster of size subCluster_size
% that is a subcluster of the ith cluster of size clusterSize
%
% SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
% the jth cluster of size subCluster_size that is a subcluster of
% the ith ccluster of size clusterSize.
% In a valid cluster, all indices must be natural numbers.


Indices = zeros(nchoosek(reference_clusterSize, ceil(reference_clusterSize/2)),reference_clusterSize);
jCluster = 0;

% Each cluster is a subset of itself.
Indices(1,clusterSize)=iCluster;
    
numberSubClusters = NchooseK(clusterSize,1:clusterSize);



switch clusterSize
  case 1
    % All non-empty subsets have been found.
    return;
    
  case 3
    possibleSubClusters_2 = zeros(numberSubClusters(2),clusterSize);
%     possibleSubClusters_2 = [0,1,1; 1, 0,1; 1,1,0];
    
  case 4
    possibleSubClusters_2 = zeros(numberSubClusters(2),clusterSize);
    possibleSubClusters_3 = zeros(numberSubClusters(3),clusterSize);
%     
%     possibleSubClusters_2 = [0,0,1,1; 0,1,0,1; 1,0,0,1; 0,1,1,0; 1,0,1,0; 1,1,0,0];
%     possibleSubClusters_3 = [0,1,1,1; 1,0,1,1; 1,1,0,1; 1,1,1,0];
   
  case 5
    possibleSubClusters_2 = zeros(numberSubClusters(2),clusterSize);
    possibleSubClusters_3 = zeros(numberSubClusters(3),clusterSize); 
    possibleSubClusters_4 = zeros(numberSubClusters(4),clusterSize); 
    
  case 6
    possibleSubClusters_2 = zeros(numberSubClusters(2),clusterSize);
    possibleSubClusters_3 = zeros(numberSubClusters(3),clusterSize); 
    possibleSubClusters_4 = zeros(numberSubClusters(4),clusterSize);  
    possibleSubClusters_5 = zeros(numberSubClusters(5),clusterSize); 
end  

subcluster_count = zeros(1,clusterSize);
for isc = 1:2^(clusterSize)-1
  
  subcluster_str = dec2bin(isc);
  subcluster_str = pad(subcluster_str,clusterSize,'left','0');
  sc_size = sum(subcluster_str=='1');
  
  switch sc_size
    
    case 2
      subcluster_count(sc_size) = subcluster_count(sc_size) + 1; 
      possibleSubClusters_2(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 3
      subcluster_count(sc_size) = subcluster_count(sc_size) + 1; 
      possibleSubClusters_3(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 4
      subcluster_count(sc_size) = subcluster_count(sc_size) + 1; 
      possibleSubClusters_4(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 5
      subcluster_count(sc_size) = subcluster_count(sc_size) + 1; 
      possibleSubClusters_5(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
  end
  
end

Cluster = Clusters(iCluster, 1:clusterSize ,clusterSize);
Indices(1:clusterSize,1) = Cluster;

% for jCluster = 1:numberSubClusters(1)
%   
%   SubCluster = Cluster(possibleSubClusters_1==1);
%   
%   % Search for a valid cluster that equals the subcluster.
%   Search = Clusters(:, 1:clusterSize -1, clusterSize -1)==SubCluster;
%  
%   % Locate where the subcluster is.
%   subclusterIndex = find( all(Search,2));
%   
%   Indices(jCluster,1) = subclusterIndex;
% end

if clusterSize == 2
  return;
end
  
for jCluster = 1:numberSubClusters(2)
  
  SubCluster = Cluster(possibleSubClusters_2(jCluster,:)==1);
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:2, 2)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2));
 
  if isempty(subclusterIndex)
    subclusterIndex = 0;
  end
  
  Indices(jCluster,2) = subclusterIndex;
end

if clusterSize == 3
  return;
end
  
for jCluster = 1:numberSubClusters(3)
  
  SubCluster = Cluster(possibleSubClusters_3(jCluster,:)==1);
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:3, 3)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2));
 
  if isempty(subclusterIndex)
    subclusterIndex = 0;
  end
  
  Indices(jCluster,3) = subclusterIndex;
end

%--------------------------------------------------------------------------

isize = 4;  
if clusterSize == isize
  return;
end

for jCluster = 1:numberSubClusters(isize)
  
  SubCluster = Cluster(possibleSubClusters_4(jCluster,:)==1);
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:isize, isize)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2));
 
  if isempty(subclusterIndex)
    subclusterIndex = 0;
  end
  
  Indices(jCluster,isize) = subclusterIndex;
end
%--------------------------------------------------------------------------

isize = 5;  
if clusterSize == isize
  return;
end

for jCluster = 1:numberSubClusters(isize)
  
  SubCluster = Cluster(possibleSubClusters_5(jCluster,:)==1);
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:isize, isize)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2));
 
  if isempty(subclusterIndex)
    subclusterIndex = 0;
  end
  
  Indices(jCluster,isize) = subclusterIndex;
end

end

function Indices = findSubclusters_gpu_TBD(Clusters,clusterSize,iCluster,reference_clusterSize)
% from Indices{subCluster_size} = list of all jCluster such that Clusters{subCluster_size}(jCluster,:) is a subcluster of Clusters{clusterSize}(iCluster,:)
% to
% Given the ith cluster of size clusterSize, 
% Indices(jCluster,subCluster_size) = jth cluster of sizesubCluster_size
% that is a subcluster of the ith cluster of size clusterSize
%
% SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
% the jth cluster of size subCluster_size that is a subcluster of
% the ith ccluster of size clusterSize.
% In a valid cluster, all indices must be natural numbers.


Indices = zeros(nchoosek(reference_clusterSize, ceil(reference_clusterSize/2)),reference_clusterSize);
jCluster = 0;

% Each cluster is a subset of itself.
Indices(1,clusterSize)=iCluster;

switch clusterSize
  case 1
    % All non-empty subsets have been found.
    return;
    
  case 3
    
    numberSubClusters = [3,3,1];
   
    possibleSubClusters_2 = [0,1,1; 1, 0,1; 1,1,0];
    
  case 4
    
    numberSubClusters = [4,6,4,1];
    
    possibleSubClusters_2 = [0,0,1,1; 0,1,0,1; 1,0,0,1; 0,1,1,0; 1,0,1,0; 1,1,0,0];
    possibleSubClusters_3 = [0,1,1,1; 1,0,1,1; 1,1,0,1; 1,1,1,0];
end  

Cluster = Clusters(iCluster, 1:clusterSize ,clusterSize);
Indices(1:clusterSize,1) = Cluster;

% for jCluster = 1:numberSubClusters(1)
%   
%   SubCluster = Cluster(possibleSubClusters_1==1);
%   
%   % Search for a valid cluster that equals the subcluster.
%   Search = Clusters(:, 1:clusterSize -1, clusterSize -1)==SubCluster;
%  
%   % Locate where the subcluster is.
%   subclusterIndex = find( all(Search,2));
%   
%   Indices(jCluster,1) = subclusterIndex;
% end

if clusterSize == 2
  return;
end
  
for jCluster = 1:numberSubClusters(2)
  
  SubCluster = Cluster(possibleSubClusters_2(jCluster,:)==1);
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:2, 2)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2));
 
  if isempty(subclusterIndex)
    subclusterIndex = 0;
  end
  
  Indices(jCluster,2) = subclusterIndex;
end

if clusterSize == 3
  return;
end
  
for jCluster = 1:numberSubClusters(3)
  
  SubCluster = Cluster(possibleSubClusters_3(jCluster,:)==1);
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:3, 3)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2));
 
  if isempty(subclusterIndex)
    subclusterIndex = 0;
  end
  
  Indices(jCluster,3) = subclusterIndex;
end


end