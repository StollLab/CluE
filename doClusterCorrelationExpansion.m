function [Signals, AuxiliarySignals] = ...
  doClusterCorrelationExpansion(...
  Coherences, Clusters, SubclusterIndices, ...
  timepoints,dimensionality, order,numberClusters)
%   Coherences{n(iCluster,timepoints)

% SubclusterIndices{clusterSize} (jCluster,subCluster_size, iCluster) =
% the jth cluster of size subCluster_size that is a subcluster of
% the ith cluster of size clusterSize.
% In a valid cluster, all indices must be natural numbers.


% Initialize data.
AuxiliarySignals = Coherences;
Signals = ones(order,timepoints^dimensionality);

%--------------------------------------------------------------------------
% n-Clusters 
%--------------------------------------------------------------------------

for isize = 1:order
  
  for iCluster=1:numberClusters(isize)
    
    % Remove m-cluster correlations
    for jsize = 1:(isize-1)
      
      for subcluster_index = SubclusterIndices{isize}(:,jsize,iCluster)'
        
        if subcluster_index <= 0
          continue;
        end
        
        AuxiliarySignals{isize}(iCluster,:) = ...
          AuxiliarySignals{isize}(iCluster,:) ...
          ./AuxiliarySignals{jsize}(subcluster_index,:);
        
      end
      
      
    end
    
    % Update cluster signals.
    for iorder=isize:order
      this_cluster = Clusters(iCluster,1:isize,isize);
      if this_cluster(1)==0
        continue;
      end
      Signals(iorder,:) = Signals(iorder,:).*AuxiliarySignals{isize}(iCluster,:);
    end
    
  end
end

end