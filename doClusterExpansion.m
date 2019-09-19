function [Signal,AuxiliarySignal,Order_n_Signal] = doClusterExpansion(Coherences,Clusters, SubclusterIndices,order)

AuxiliarySignal = Coherences;
ContributionSum = zeros(size(Coherences{1,1}));
Order_n_Signal{1} = ones(size(Coherences{1,1}));

% Calculate signal from 2-clusters
if order>=2
  clusterSize = 2;
  for iCluster = 1:size(Clusters{clusterSize},1)
    if ~checkCluster(Clusters, SubclusterIndices,iCluster,clusterSize)
      AuxiliarySignal{clusterSize,iCluster} = 0*AuxiliarySignal{clusterSize,iCluster};
      continue
    end
    sci = SubclusterIndices{clusterSize,iCluster};
    v{2} = Coherences{clusterSize-1, sci{clusterSize - 1}(2) };
    v{1} = Coherences{clusterSize-1, sci{clusterSize - 1}(1)};
    AuxiliarySignal{clusterSize,iCluster} = Coherences{clusterSize,iCluster} - v{1}.*v{2};
    ContributionSum = ContributionSum + AuxiliarySignal{clusterSize,iCluster};
    AuxiliarySignal{clusterSize,iCluster} = exp(AuxiliarySignal{clusterSize,iCluster});
    clear v
  end
end

% calculate 2-CE
Order_n_Signal{2} = exp(ContributionSum);

% Calculate signal from 3-clusters
if order>=3
  clusterSize = 3;
  for iCluster = 1:size(Clusters{clusterSize},1)
    if ~checkCluster(Clusters, SubclusterIndices,iCluster,clusterSize)
      AuxiliarySignal{clusterSize,iCluster} = 0*AuxiliarySignal{clusterSize,iCluster};
      continue
    end
    sci = SubclusterIndices{clusterSize,iCluster};
    v12 = AuxiliarySignal{clusterSize-1, sci{clusterSize - 1}(1)};
    v13 = AuxiliarySignal{clusterSize-1, sci{clusterSize - 1}(2)};
    v23 = AuxiliarySignal{clusterSize-1, sci{clusterSize - 1}(3)};
    v3  = AuxiliarySignal{clusterSize-2, sci{clusterSize - 2}(3)};
    v2  = AuxiliarySignal{clusterSize-2, sci{clusterSize - 2}(2)};
    v1  = AuxiliarySignal{clusterSize-2, sci{clusterSize - 2}(1)};
    AuxiliarySignal{clusterSize,iCluster} = Coherences{clusterSize,iCluster} - v12.*v3 - v13.*v2 - v23.*v1 - v1.*v2.*v3;
    ContributionSum = ContributionSum + AuxiliarySignal{clusterSize,iCluster};
    AuxiliarySignal{clusterSize,iCluster} = exp(AuxiliarySignal{clusterSize,iCluster});
    clear v1 v2 v3 v12 v13 v23
  end
end

% calculate 3-CE
Order_n_Signal{3} = exp(ContributionSum);

% Calculate signal from 4-clusters
if order>=4
  clusterSize = 4;
  for iCluster = 1:size(Clusters{clusterSize},1)
    if ~checkCluster(Clusters, SubclusterIndices,iCluster,clusterSize)
      AuxiliarySignal{clusterSize,iCluster} = 0*AuxiliarySignal{clusterSize,iCluster};
      continue
    end
    sci = SubclusterIndices{clusterSize,iCluster};
    v123 = AuxiliarySignal{clusterSize-1, sci{clusterSize - 1}(1)};
    v124 = AuxiliarySignal{clusterSize-1, sci{clusterSize - 1}(2)};
    v134 = AuxiliarySignal{clusterSize-1, sci{clusterSize - 1}(3)};
    v234 = AuxiliarySignal{clusterSize-1, sci{clusterSize - 1}(4)};
    v12  = AuxiliarySignal{clusterSize-2, sci{clusterSize - 2}(1)};
    v13  = AuxiliarySignal{clusterSize-2, sci{clusterSize - 2}(2)};
    v14  = AuxiliarySignal{clusterSize-2, sci{clusterSize - 2}(3)};
    v23  = AuxiliarySignal{clusterSize-2, sci{clusterSize - 2}(4)};
    v24  = AuxiliarySignal{clusterSize-2, sci{clusterSize - 2}(5)};
    v34  = AuxiliarySignal{clusterSize-2, sci{clusterSize - 2}(6)};
    v1   = AuxiliarySignal{clusterSize-3, sci{clusterSize - 3}(1)};
    v2   = AuxiliarySignal{clusterSize-3, sci{clusterSize - 3}(2)};
    v3   = AuxiliarySignal{clusterSize-3, sci{clusterSize - 3}(3)};
    v4   = AuxiliarySignal{clusterSize-3, sci{clusterSize - 3}(4)};
    AuxiliarySignal{clusterSize,iCluster} = Coherences{clusterSize,iCluster} ...
      - v123.*v4-v124.*v3 - v134.*v2 - v234.*v1 ...
      - v12.*v34 - v13.*v24 - v14.*v23 ...
      - v12.*v3.*v4 - v13.*v2.*v4 - v14.*v2.*v3 - v23.*v1.*v4 - v24.*v1.*v3 - v34.*v1.*v2 ...
      - v1.*v2.*v3.*v4;
    ContributionSum = ContributionSum + AuxiliarySignal{clusterSize,iCluster};
    AuxiliarySignal{clusterSize,iCluster} = exp(AuxiliarySignal{clusterSize,iCluster});
    clear v1 v2 v3 v4 v12 v13 v14 v23 v24 v34 v123 v124 v134 v234;
  end
end

% calculate 4-CE
Order_n_Signal{4} = exp(ContributionSum);

Signal = exp(ContributionSum);

end

function valid = checkCluster(Clusters, SubclusterIndices,iCluster,clusterSize)
valid = true;
n_subclusters = 2^(clusterSize) - 1;
for iSize = 1:clusterSize
  n_subclusters = n_subclusters - size(SubclusterIndices{clusterSize,iCluster}{iSize},2);
end
if n_subclusters~=0
  valid = false;
  return
end
Test = Clusters{clusterSize}(iCluster,:)'==Clusters{1}(SubclusterIndices{clusterSize,iCluster}{1}(:));
if ~all(all(Test))
  valid = false;
  return
end

end