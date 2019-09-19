function [Signal, AuxiliarySignal, order_n_signal] = doCCE(Coherences,Clusters, order)

% find auxiliary signals
AuxiliarySignal = getCorrelations(Coherences,Clusters, order);

% initialize full coherence signal
Signal = ones(size(AuxiliarySignal{1,1}));

% initialize nth order coherence signal
order_n_signal{order} = Signal;

% calculate observable signal
% loop over cluster sizes
for iClusterSize = 1:order
  
  % mulyiply over all auxiliary signals from size iClusterSize clusters
  for jCluster = 1:size(Clusters{iClusterSize},1)
    Signal = Signal.*AuxiliarySignal{iClusterSize,jCluster};
  end
  
  % save the n-CCE signal
  order_n_signal{iClusterSize} = Signal;
end

end

function Correlations = getCorrelations(Coherences,Clusters, order)
Correlations = Coherences;
for iClusterSize = 2:order
  for jCluster = 1:size(Clusters{iClusterSize},1)
    Correlations{iClusterSize,jCluster} = removeSubClusterCorrelations(Correlations,Clusters,iClusterSize,jCluster);
  end
end
end

function Correlation = removeSubClusterCorrelations(CorrelationSet,Clusters,clusterSize,jCluster)

Correlation = CorrelationSet{clusterSize,jCluster};
Indices = findSubclusters(Clusters,clusterSize,jCluster);

for iSize = 1:(clusterSize-1)
  for index = Indices{iSize}
    Correlation = Correlation./CorrelationSet{iSize,index};
  end
end
end


function Indices = findSubclusters(Clusters,clusterSize,iCluster)
% Indices{size} = list of all jCluster such that Clusters{size}(jCluster,:) is a subcluster of Clusters{clusterSize}(iCluster,:)
Indices{clusterSize}=iCluster;
for index = 1:clusterSize
  SubCluster = [Clusters{clusterSize}(iCluster,1:index-1),Clusters{clusterSize}(iCluster,index+1:end)];
  % Search for the subcluster.
  Search = Clusters{clusterSize -1}==SubCluster;
  for ii = 2:size(Search,2)
    Search(:,1) = Search(:,1).*Search(:,ii);
  end
  subclusterIndex = find(Search(:,1)==1);
  
  if isempty(subclusterIndex)
    continue;
  end
  Indices{clusterSize -1} = [Indices{clusterSize -1} , subclusterIndex];
  
  if clusterSize > 2
    Subindices = findSubclusters(Clusters,clusterSize -1,subclusterIndex);
    for isize = (clusterSize-2):-1:1
      if ~isempty(Subindices{isize})
        Indices{isize} = [Indices{isize},Subindices{isize}];
      end
    end
  end
  
end

for isize = 1:clusterSize
  Indices{isize} = unique(Indices{isize});
end
end

