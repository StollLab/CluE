% ========================================================================
% New Function
% ========================================================================

% Subset Enumeration via the Banker Method
function newCluster = getNextCluster(Cluster,numberSpins)

% Initialize the output.
newCluster = Cluster;

% Determine the size of the cluster.
clusterSize = length(Cluster);

% Start at final spin.
for icluster = clusterSize:-1:1
  
  % Check if the element can be incremented
  if Cluster(icluster) < (numberSpins + icluster- clusterSize )
    
    % Increment element.
    newCluster(icluster) = newCluster(icluster) + 1;
    for jcluster = icluster+1:clusterSize
      % Reset the elements after the incremented element.
      newCluster(jcluster) = newCluster(jcluster-1) + 1;
    end
    
    % Return output.
    return;
  end
end

% There is no next cluster; return the empty set.
newCluster = [];

end