function ClusterGeo = getClusterGeoStats(Clusters,Coordinates, DistanceMatrix,numberClusters,Method_order) 

ClusterGeo = cell(1,Method_order);

for isize = Method_order:-1:1
  
  nC_ = numberClusters(isize);
  C_ = 1:isize;
  ClusterGeo{isize}.number = nC_;
  ClusterGeo{isize}.cluster_indices = 1:nC_;
  ClusterGeo{isize}.cluster_distance = max(Coordinates.r(Clusters{isize}),[],2); 
  

  ClusterGeo{isize}.cluster_proximity = min(Coordinates.r(Clusters{isize}),[],2);
  ClusterGeo{isize}.cluster_distance_indices = diag( ...
    Clusters{isize}(ClusterGeo{isize}.cluster_indices', ... 
      sum(  (Coordinates.r(Clusters{isize}) == ClusterGeo{isize}.cluster_distance) ... 
        .*C_ , 2 ) ...
      ) ...
    );
  ClusterGeo{isize}.cluster_proximity_indices = diag(Clusters{isize}(ClusterGeo{isize}.cluster_indices', sum(  (Coordinates.r(Clusters{isize}) == ClusterGeo{isize}.cluster_proximity).*C_ ,2 )));

  r1_ = ClusterGeo{isize}.cluster_distance;
  ind1_ = ClusterGeo{isize}.cluster_distance_indices;
  
  r2_ = ClusterGeo{isize}.cluster_proximity;
  ind2_ = ClusterGeo{isize}.cluster_proximity_indices;
  if isize==1
    continue;
  end
  
  ClusterGeo{isize}.cluster_separation = diag(DistanceMatrix(ind1_, ind2_));

  R_ = ClusterGeo{isize}.cluster_separation;
  % R^2 = r1^2 + r2^2 -2*r1*r2*cos(theta).
  % 2*r1*r2*cos(theta) = r1^2 + r2^2 -R^2.
  % cos(theta) = (r1^2 + r2^2 -R^2)/(2*r1*r2).
  ClusterGeo{isize}.cluster_cos = (r1_.^2 + r2_.^2 - R_.^2)./(2*r1_.*r2_);

end

end