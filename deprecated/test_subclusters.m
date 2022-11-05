function pass = test_subclusters(ClusterArray, SubclusterIndices_2,SubclusterIndices_3,SubclusterIndices_4,verbose)
% SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
% the jth cluster of size subCluster_size that is a subcluster of
% the ith cluster of size clusterSize.
% In a valid cluster, all indices must be natural numbers.

% [number,maxSize,maxClusterSize]=size(ClusterArray);
maxClusterSize=size(ClusterArray,3);
hline='------------------------------------------------------------------';
pass = false;
pass_count = 0;
zero_count = 0;
fail_count = 0;

clusterSize = 2;
if maxClusterSize < clusterSize
  pass = true;
end

% Loop through all subclusters of supercluster size 2.
for jCluster = 1:size(SubclusterIndices_2,1)
  
  % Loop through all subclusters sizes.
  for subCluster_size = 1:size(SubclusterIndices_2,2)
    
    % Loop through all supercluster with the particular subcluster.
    for iCluster = 1:size(SubclusterIndices_2,3)
      
      % Get the index of the the jth subCluster_size cluster 
      % of the ith supercluster.
      index = SubclusterIndices_2(jCluster,subCluster_size,iCluster);
      
      % Skip zero indices. 
      if index ==0
        zero_count = zero_count + 1; 
        continue;
      end
      
      % Get subcluster.
      Subcluster = ClusterArray(index,:,subCluster_size);
      
      % Remove zeros.
      Subcluster(Subcluster==0) = [];
      
      % Get supercluster.
      Cluster = ClusterArray(iCluster, :, clusterSize);
      
      % Check if Subcluster is contained in Cluster.  
      if ~isSubcluster(Subcluster,Cluster)
        if verbose
          disp(hline);
          fprintf('|C| = %d.\n',clusterSize);
          fprintf('index = %d.\n',index);
          disp('Subcluster');
          disp(Subcluster);
          disp('Cluster');
          disp(Cluster);
          disp(hline);
        end
        fail_count = fail_count +1;
      else
        pass_count = pass_count + 1;
      end
      
    end
  end
end

% Loop through all subclusters of supercluster size 3.
clusterSize = 3;
if maxClusterSize < clusterSize
  pass = true;
end

% Loop through all supercluster with the particular subcluster.
for jCluster = 1:size(SubclusterIndices_3,1)
  
  % Loop through all subclusters sizes.
  for subCluster_size = 1:size(SubclusterIndices_3,2)
    
    % Loop through all supercluster with the particular subcluster.
    for iCluster = 1:size(SubclusterIndices_3,3)
      
      % Get the index of the the jth subCluster_size cluster
      % of the ith supercluster.
      index = SubclusterIndices_3(jCluster,subCluster_size,iCluster);
      
      % Skip zero indices. 
      if index ==0
        zero_count = zero_count + 1; 
        continue;
      end
      
      % Get subcluster.
      Subcluster = ClusterArray(index,:,subCluster_size);
      
      % Remove zeros.
      Subcluster(Subcluster==0) = [];
      
      % Get supercluster.
      Cluster = ClusterArray(iCluster, :, clusterSize);
      
      % Check if Subcluster is contained in Cluster.
      
      if ~isSubcluster(Subcluster,Cluster)
        if verbose
          disp(hline);
          fprintf('|C| = %d.\n',clusterSize);
          fprintf('index = %d.\n',index);
          disp('Subcluster');
          disp(Subcluster);
          disp('Cluster');
          disp(Cluster);
          disp(hline);
        end
        fail_count = fail_count +1;
      else
        pass_count = pass_count + 1;
      end
      
    end
  end
end

clusterSize = 4;
if maxClusterSize < clusterSize
  pass = true;
end

% Loop through all supercluster with the particular subcluster.
for jCluster = 1:size(SubclusterIndices_4,1)
  
  % Loop through all subclusters sizes.
  for subCluster_size = 1:size(SubclusterIndices_4,2)
    
    % Loop through all supercluster with the particular subcluster.
    for iCluster = 1:size(SubclusterIndices_4,3)
      
      % Get the index of the the jth subCluster_size cluster
      % of the ith supercluster.
      index = SubclusterIndices_4(jCluster,subCluster_size,iCluster);
     
      % Skip zero indices. 
      if index ==0
        zero_count = zero_count + 1; 
        continue;
      end
      
      % Get subcluster.
      Subcluster = ClusterArray(index,:,subCluster_size);
      
      % Remove zeros.
      Subcluster(Subcluster==0) = [];
      
      % Get supercluster.
      Cluster = ClusterArray(iCluster, :, clusterSize);
      
      % Check if Subcluster is contained in Cluster.
      if ~isSubcluster(Subcluster,Cluster)
        if verbose
          disp(hline);
          fprintf('|C| = %d.\n',clusterSize);
          fprintf('index = %d.\n',index);
          disp('Subcluster');
          disp(Subcluster);
          disp('Cluster');
          disp(Cluster);
          disp(hline);
        end
        fail_count = fail_count +1;
      else
        pass_count = pass_count + 1;
      end
      
    end
  end
end

% Print summary.
if verbose
  disp(hline);
  fprintf('Passes = %d.\n',pass_count);
  fprintf('Fails = %d.\n',fail_count);
  fprintf('Zeros = %d.\n',zero_count);
  disp(hline);
end

% Determine if subclusters indices are valid.
if fail_count == 0
  pass = true;
end

end