function [Clusters,Nuclei] = getClusters(Nuclei,Method,Data)

Clusters = cell(Method.extraOrder,1);
for clusterSize = 1:Method.extraOrder
  Clusters{clusterSize} = [];
end

inClusters = {};
doFindClusters = ~Method.Ori_cutoffs;
if strcmp(Method.method,'count clusters')
  doFindClusters = true;
end

if ~isempty(Data.ClusterData)  || isfield(Method,'Clusters')
  Clusters = {};
  if ~isempty(Data.ClusterData)
    try
      load(Data.ClusterData,'Clusters');      
    catch
      disp('Could not load clusters.')      
      if Data.exitOnFailedLoad
        error('Could not load clusters.')
      end
    end
  else
    Clusters = Method.Clusters;
  end
   
  inOrder = numel(Clusters); 
   
  for clusterSize = 1:inOrder
    
    Nuclei.numberClusters(clusterSize) = size(Clusters{clusterSize},1);
    if Method.verbose
      fprintf('Loaded %d clusters of size %d.\n', Nuclei.numberClusters(clusterSize),clusterSize);
    end
   
  end
  if inOrder < Method.order
    inClusters = Clusters;
  else
    doFindClusters = isempty(Clusters);
  end
end


% Loop over cluster sizes, start at the largest (most time consuming) size
if doFindClusters

  if Method.verbose
    fprintf('Finding clusters of up to size %d.\n',Method.extraOrder);
  end

  if Method.combineClusters
    Clusters = findClusters_treeSearch(Nuclei,Method.extraOrder,1,{});
    for clusterSize = 1:min(Method.order, numel(inClusters),Method)
      % Combine arrays.
      C = [Clusters{clusterSize}; inClusters{clusterSize}];

      % Sort clusters
      C = sortrows(C);

      % Remove duplicates
      keep = [true; any(C(1:end-1,:)~=C(2:end,:),2)];
      Clusters{clusterSize} = C(keep,:);
    end
  else
    if ~Method.neighborCutoff.sizeDependent
      Clusters = findClusters_treeSearch(Nuclei,Method.extraOrder,1,...
        inClusters, Method);
    else
      Clusters = cell(1,Method.extraOrder);
      for adjacencyOrder = Method.extraOrder:-1:1
        Clusters_ = findClusters_treeSearch(Nuclei,Method.extraOrder,...
          adjacencyOrder,inClusters, Method);
        Clusters{adjacencyOrder} = Clusters_{adjacencyOrder};
      end
    end
  end
  Nuclei.numberClusters = zeros(1,Method.extraOrder);

  for clusterSize = 1:Method.extraOrder
    Nuclei.numberClusters(clusterSize) = size(Clusters{clusterSize},1);

    fprintf('  Found %d clusters of size %d.\n', ...
      Nuclei.numberClusters(clusterSize),clusterSize);
  end


end
end