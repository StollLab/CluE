function [ClusterArray, SubclusterIndices_2,SubclusterIndices_3,SubclusterIndices_4,SubclusterIndices_5,SubclusterIndices_6] = cluster4(order)

ClusterArray = zeros(nchoosek(order,ceil(order/2)),order,order);

numClusters = zeros(1,4);
numClusters(1) = order;

switch order
  case 2
    ClusterArray(1:2,1,1) = 1:2;
    ClusterArray(1,:,2) = 1:2;
    
  case 3
    ClusterArray(1:3,1,1) = 1:3;
    ClusterArray(:,:,2) = [1,2,0; 1,3,0; 2,3,0];
    
    ClusterArray(1,:,3) = 1:3;
    
  case 4
    
    
    ClusterArray(1:4,1,1)=1:4;
    
    ClusterArray(1:3,1,2)=1; ClusterArray(1:3,2,2)=2:4;
    ClusterArray(4:5,1,2)=2; ClusterArray(4:5,2,2)=3:4;
    ClusterArray(6,1:2,2) = [3,4];
    
    ClusterArray(1,1:3,3) = 1:3;
    ClusterArray(2,1:3,3) = [1,2,4];
    ClusterArray(3,1:3,3) = [1,3,4];
    ClusterArray(4,1:3,3) = 2:4;
    
    ClusterArray(1,1:4,4) = 1:4;
    
    
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
  SubclusterIndices_2(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
end

%--------------------------------------------------------------------------
clusterSize = 3;

if order < clusterSize
  SubclusterIndices_3 = [];
  SubclusterIndices_4 = [];
  return;
end

numClusters(clusterSize) = nchoosek(order, clusterSize);
SubclusterIndices_3 = zeros(nchoosek(clusterSize,1), clusterSize , numClusters(clusterSize));

for iCluster = 1:numClusters(clusterSize)
  SubclusterIndices_3(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
end

%--------------------------------------------------------------------------
clusterSize = 4;

if order < clusterSize
  SubclusterIndices_4 = [];
  return;
end

numClusters(clusterSize) = nchoosek(order, clusterSize);
SubclusterIndices_4 = zeros(nchoosek(clusterSize,ceil(clusterSize/2)), clusterSize , numClusters(clusterSize));

for iCluster = 1:numClusters(clusterSize)
  SubclusterIndices_4(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
end
%--------------------------------------------------------------------------
clusterSize = 5;

if order < clusterSize
  SubclusterIndices_5 = [];
  return;
end

numClusters(clusterSize) = nchoosek(order, clusterSize);
SubclusterIndices_5 = zeros(nchoosek(clusterSize,ceil(clusterSize/2)), clusterSize , numClusters(clusterSize));

for iCluster = 1:numClusters(clusterSize)
  SubclusterIndices_5(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
end

%--------------------------------------------------------------------------
clusterSize = 6;

if order < clusterSize
  SubclusterIndices_6 = [];
  return;
end

numClusters(clusterSize) = nchoosek(order, clusterSize);
SubclusterIndices_6 = zeros(nchoosek(clusterSize,ceil(clusterSize/2)), clusterSize , numClusters(clusterSize));

for iCluster = 1:numClusters(clusterSize)
  SubclusterIndices_6(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
end

end
