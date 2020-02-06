function pass = test_findSubClusters_gpu()
tic
pass = true;
oldpath = path;
path('../',oldpath);

% ENUM
FULL = 1;  RANDOM = 2;
for order = 1:6
  N = 2^(9-order);
  Nuclei.number = N;
  
  for adjType = FULL:RANDOM
    switch adjType
      case FULL
        Nuclei.ValidPair = ones(Nuclei.number);
      case RANDOM
        % Get random adjacency matrix.
        Nuclei.ValidPair = (rand(Nuclei.number)-rand() + eye(Nuclei.number))>0;
        Nuclei.ValidPair = (Nuclei.ValidPair + Nuclei.ValidPair') - rand() > 0;
        if mma(Nuclei.ValidPair - Nuclei.ValidPair')>0
          disp(mma(Nuclei.ValidPair - Nuclei.ValidPair'));
          error('Nuclei.ValidPair is not symmetric.');
        end
    end
    
    % Get clusters.
    numClusters = NchooseK(N,1:order);
    Nc = numClusters;
    Clusters = zeros(numClusters(order),order, order);
    for clusterSize = 1:order
      c_ = findClusters(Nuclei,clusterSize);
      Nc(clusterSize) = size(c_,1);
      Clusters(1:Nc(clusterSize), 1:clusterSize, clusterSize ) = c_;
    end
    Clusters(max(Nc)+1:end,:,:) = [];
    
    if adjType == FULL
      if ~all(Nc == numClusters)
        pass = false;
        disp('Not all subclusters were found.');
        return;
      end
    end
    
    
    % Loop over all super cluster sizes.
    for super_clusterSize = 1:order
      
      % Loop over possible cluster sizes.
      for clusterSize = 1:super_clusterSize
        
        % Loop over clusters.
        for iCluster = 1:Nc(clusterSize)
          
          % Given the ith cluster of size clusterSize,
          % Indices(jCluster,subCluster_size) = jth cluster of size subCluster_size
          % that is a subcluster of the ith cluster of size clusterSize
          SubclusterList = findSubclusters_gpu(Clusters,clusterSize,iCluster,super_clusterSize);
          
          % Get cluster for testing.
          test_cluster = Clusters(iCluster,1:clusterSize,clusterSize);
          
          % Loop over subcluster sizes.
          for jsize = 1:clusterSize
            
            % Loop over subclusters of size jsize.
            for jCluster = 1:size(SubclusterList,1)
              
              % Get index pointing to the jth subcluster.
              subcluster_index = SubclusterList(jCluster,jsize);
              
              if subcluster_index == 0
                continue;
              end
              
              subcluster = Clusters(subcluster_index,1:jsize,jsize);
              
              if ~isSubcluster(subcluster,test_cluster)
                pass = false;
                fprintf('test_findSubClusters_gpu() failed.\n')
                return;
              end
              
            end
          end
          
        end
      end
    end
  end
end
toc
end
function a = mma(A)
a = max(max(abs(A)));
end