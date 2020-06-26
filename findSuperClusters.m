function SuperClusterIndices = findSuperClusters(ClusterArray,clusterSize,iCluster)

% Extract cluster of interest.
thisCluster = ClusterArray(iCluster,1:clusterSize,clusterSize);

% Get maximum cluster size in ClusterArray. 
maxClusterSize = size(ClusterArray,3);

% Determine the highest possible number of super-cluster.
maxNumberOfClusters = size(ClusterArray,1);

% Number line over cluster index range.
Indices = (1:maxNumberOfClusters)';

% Initialize SuperClusterIndices.
SuperClusterIndices_ = zeros(maxNumberOfClusters,maxClusterSize);

% Each cluster is a superset of itself.
SuperClusterIndices_(1,clusterSize)=iCluster;

Length = 1;

% Loop over all cluster sizes greater than thisCluster.
for isize = clusterSize+1:maxClusterSize
  
  % Assume all cluster could be super-clusters.
  potenstialCluster = true(maxNumberOfClusters,1);
  
  % Loop over all spins in thisCluster.
  for ii = 1:clusterSize
    % Get the spin inddex.
    ispin = thisCluster(ii);
    
    % Remove all clusters that do not contain ispin from examination.
    potenstialCluster = potenstialCluster & any(ClusterArray(:,1:isize,isize)==ispin,2);
  end
  
  % Generate list of indices for the super-clusters of size isize.
  idx_ = Indices(potenstialCluster);
  size_idx_ = size(idx_,1);
  Length = max(size_idx_,Length);
  % Record indices.
  SuperClusterIndices_(1:size_idx_, isize) = idx_;

end

% Prune SuperClusterIndices.
% SuperClusterIndices_( all(SuperClusterIndices_  == zeros(1,maxClusterSize),2),:  ) = [];
SuperClusterIndices = SuperClusterIndices_(1:Length,:);
end