function [ReducedClusterArray, SubclusterIndices_2, ...
          SubclusterIndices_3,SubclusterIndices_4,...
          SubclusterIndices_5,SubclusterIndices_6] = getReducedClusters_4(order)

ReducedClusterArray = zeros(nchoosek(order,ceil(order/2)),order,order);
maxSize = 6;

numClusters = NchooseK(order,1:maxSize);

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
    
  case 5
 
    ReducedClusterArray(1:5,1,1) = 1:5;
    
    ReducedClusterArray(1:10,1:2,2) = [1,2; 1,3; 1,4; 1,5; ...
                                  2,3; 2,4; 2,5; ...
                                  3,4; 3,5; ...
                                  4,5];
                                
    ReducedClusterArray(1:10,1:3,3) = [ 1,2,3; 1,2,4; 1,2,5; ...
                                   1,3,4; 1,3,5; 1,4,5; ...
                                   2,3,4; 2,3,5; 2,4,5; ...
                                   3,4,5];
                                 
    ReducedClusterArray(1:5,1:4,4) = [ 1,2,3,4; 1,2,3,5; 1,2,4,5; ...
                                   1,3,4,5; 2,3,4,5];
                                 
    ReducedClusterArray(1,1:5,5) = 1:5;
    
    case 6
 
    ReducedClusterArray(1:6,1,1) = 1:6;
    
    ReducedClusterArray(1:15,1:2,2) = [1,2; 1,3; 1,4; 1,5; 1,6;...
                                  2,3; 2,4; 2,5; 2,6;...
                                  3,4; 3,5; 3,6;...
                                  4,5; 4,6; ...
                                  5,6];
                                
    ReducedClusterArray(1:20,1:3,3) = [ 1,2,3; 1,2,4; 1,2,5; 1,2,6;...
                                        1,3,4; 1,3,5; 1,3,6; ...
                                        1,4,5; 1,4,6;...
                                        1,5,6; ...
                                        2,3,4; 2,3,5; 2,3,6; ...
                                        2,4,5; 2,4,6; ...
                                        2,5,6; ...
                                        3,4,5; 3,4,6; ...
                                        3,5,6; ...
                                        4,5,6];
                                 
    ReducedClusterArray(1:15,1:4,4) = [ 1,2,3,4; 1,2,3,5; 1,2,3,6;...
                                      1,2,4,5; 1,2,4,6; ...
                                      1,2,5,6; ...
                                      1,3,4,5; 1,3,4,6; ...
                                      1,3,5,6; ...
                                      1,4,5,6; ...
                                      2,3,4,5; 2,3,4,6; ...
                                      2,3,5,6;
                                      2,4,5,6;...
                                      3,4,5,6];
                                    
    ReducedClusterArray(1:6,1:5,5) = [1,2,3,4,5; 1,2,3,4,6; 1,2,3,5,6;...
                                      1,2,4,5,6; 1,3,4,5,6; ...
                                      2,3,4,5,6];   
                                    
    ReducedClusterArray(1,1:6,6) = 1:6;
    
    
end

%--------------------------------------------------------------------------
clusterSize = 2;

if order < clusterSize
  SubclusterIndices_2 = [];
  SubclusterIndices_3 = [];
  SubclusterIndices_4 = [];
  SubclusterIndices_5 = [];
  SubclusterIndices_6 = [];
  return;
end

SubclusterIndices_2 = zeros(nchoosek(clusterSize,1), clusterSize , numClusters(clusterSize));

for iCluster = 1:numClusters(clusterSize)
  SubclusterIndices_2(:,:,iCluster) = findSubclusters_gpu(ReducedClusterArray,clusterSize,iCluster,clusterSize);
end

%--------------------------------------------------------------------------
clusterSize = 3;

if order < clusterSize
  SubclusterIndices_3 = [];
  SubclusterIndices_4 = [];
  SubclusterIndices_5 = [];
  SubclusterIndices_6 = [];
  return
end
 
SubclusterIndices_3 = zeros(nchoosek(clusterSize,1), clusterSize , numClusters(clusterSize));

for iCluster = 1:numClusters(clusterSize)
  SubclusterIndices_3(:,:,iCluster) = findSubclusters_gpu(ReducedClusterArray,clusterSize,iCluster,clusterSize);
end

%--------------------------------------------------------------------------
clusterSize = 4;

if order < clusterSize
  SubclusterIndices_4 = [];
  SubclusterIndices_5 = [];
  SubclusterIndices_6 = [];
  return
end
 
SubclusterIndices_4 = zeros(nchoosek(clusterSize,ceil(clusterSize/2)), clusterSize , numClusters(clusterSize));

for iCluster = 1:numClusters(clusterSize)
  SubclusterIndices_4(:,:,iCluster) = findSubclusters_gpu(ReducedClusterArray,clusterSize,iCluster,clusterSize);
end

%--------------------------------------------------------------------------
clusterSize = 5;

if order < clusterSize
  SubclusterIndices_5 = [];
  SubclusterIndices_6 = [];
  return
end
 
SubclusterIndices_5 = zeros(nchoosek(clusterSize,ceil(clusterSize/2)), clusterSize , numClusters(clusterSize));

for iCluster = 1:numClusters(clusterSize)
  SubclusterIndices_5(:,:,iCluster) = findSubclusters_gpu(ReducedClusterArray,clusterSize,iCluster,clusterSize);
end

%--------------------------------------------------------------------------
clusterSize = 6;

if order < clusterSize
  SubclusterIndices_6 = [];
  return
end
 
SubclusterIndices_6 = zeros(nchoosek(clusterSize,ceil(clusterSize/2)), clusterSize , numClusters(clusterSize));

for iCluster = 1:numClusters(clusterSize)
  SubclusterIndices_6(:,:,iCluster) = findSubclusters_gpu(ReducedClusterArray,clusterSize,iCluster,clusterSize);
end

end
