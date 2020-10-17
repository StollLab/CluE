function SuperclusterUnion = findSuperClusterUnion(ClusterArray,clusterSize,iCluster)

% Extract cluster of interest.
thisCluster = ClusterArray(iCluster,1:clusterSize,clusterSize);

% Get maximum cluster size in ClusterArray. 
maxClusterSize = size(ClusterArray,3);

% Assume all cluster could be super-clusters.
potentialCluster = true(size(ClusterArray));
for isize = 1:clusterSize
  potentialCluster(:,:,isize) = false;
end
potentialCluster(iCluster,1:clusterSize,clusterSize) = true;

% Loop over all cluster sizes greater than thisCluster.
for isize = clusterSize+1:maxClusterSize
  
  % Loop over all spins in thisCluster.
  for ii = 1:clusterSize
    % Get the spin inddex.
    ispin = thisCluster(ii);
    
    % Remove all clusters that do not contain ispin from examination.
    potentialCluster(:,1:isize,isize) = potentialCluster(:,1:isize,isize) & any(ClusterArray(:,1:isize,isize)==ispin,2);
  end

end

SuperclusterUnion = sort(unique(ClusterArray(potentialCluster)))';
if SuperclusterUnion(1)==0
  SuperclusterUnion = SuperclusterUnion(2:end);
end


end