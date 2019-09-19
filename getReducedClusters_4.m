function [ReducedClusterArray, SubclusterIndices_2,SubclusterIndices_3,SubclusterIndices_4] = getReducedClusters_4(order)

ReducedClusterArray = zeros(nchoosek(order,ceil(order/2)),order,order);

numClusters = zeros(1,4);

numClusters(1) = order;

switch order
  case 1
    ReducedClusterArray = 1;
  case 2
    ReducedClusterArray(1:2,1,1) = 1:2;
    ReducedClusterArray(1,:,2) = 1:2;
    
  case 3
    ReducedClusterArray(1:3,1,1) = 1:3;
    ReducedClusterArray(:,:,2) = [1,2,0; 1,3,0; 2,3,0];
    ReducedClusterArray(1,:,3) = 1:3;
    
  case 4
    ReducedClusterArray(1:4,1,1) = 1:4;
    ReducedClusterArray(1:3,1,2) = 1;
    ReducedClusterArray(1:3,2,2) = 2:4;
    ReducedClusterArray(4:5,1,2) = 2;
    ReducedClusterArray(4:5,2,2) = 3:4;
    ReducedClusterArray(6,1:2,2) = [3,4];
    ReducedClusterArray(1,1:3,3) = 1:3;
    ReducedClusterArray(2,1:3,3) = [1,2,4];
    ReducedClusterArray(3,1:3,3) = [1,3,4];
    ReducedClusterArray(4,1:3,3) = 2:4;
    ReducedClusterArray(1,1:4,4) = 1:4;
    
end

%--------------------------------------------------------------------------
clusterSize = 2;

if order < clusterSize
  SubclusterIndices_2 = [];
  SubclusterIndices_3 = [];
  SubclusterIndices_4 = [];
  return;
end
numClusters(clusterSize) = nchoosek(order, clusterSize);
SubclusterIndices_2 = zeros(nchoosek(clusterSize,1), clusterSize , numClusters(clusterSize));

for iCluster = 1:numClusters(clusterSize)
  SubclusterIndices_2(:,:,iCluster) = findSubclusters_gpu(ReducedClusterArray,clusterSize,iCluster,clusterSize);
end

%--------------------------------------------------------------------------
clusterSize = 3;

if order < clusterSize
  SubclusterIndices_3 = [];
  SubclusterIndices_4 = [];
  return
end

numClusters(clusterSize) = nchoosek(order, clusterSize);
SubclusterIndices_3 = zeros(nchoosek(clusterSize,1), clusterSize , numClusters(clusterSize));

for iCluster = 1:numClusters(clusterSize)
  SubclusterIndices_3(:,:,iCluster) = findSubclusters_gpu(ReducedClusterArray,clusterSize,iCluster,clusterSize);
end

%--------------------------------------------------------------------------
clusterSize = 4;

if order < clusterSize
  SubclusterIndices_4 = [];
  return
end

numClusters(clusterSize) = nchoosek(order, clusterSize);
SubclusterIndices_4 = zeros(nchoosek(clusterSize,ceil(clusterSize/2)), clusterSize , numClusters(clusterSize));

for iCluster = 1:numClusters(clusterSize)
  SubclusterIndices_4(:,:,iCluster) = findSubclusters_gpu(ReducedClusterArray,clusterSize,iCluster,clusterSize);
end

end
