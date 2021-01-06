function [Signals, ...
          AuxiliarySignal_1,AuxiliarySignal_2] = ...
  doHydrogenIsotopologueCCE(...
  Coherences_1H,Coherences_2H,Coherences_1D,Coherences_2D, fractions, Clusters, ...
  SubclusterIndices_2H,SubclusterIndices_2D,...
  timepoints,dimensionality, order,numberClusters,Exchangable,MoleculeID)
%   Coherences_n(iCluster,timepoints)

% SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
% the jth cluster of size subCluster_size that is a subcluster of
% the ith cluster of size clusterSize.
% In a valid cluster, all indices must be natural numbers.

% 
% pass = test_subclusters(Clusters, SubclusterIndices_2,SubclusterIndices_3,SubclusterIndices_4,false);
% if ~pass
%   error('Incorrect Subcluster');
% end

% Initialize data.
code_maxSize = 2;
maxSize = min(code_maxSize,order);

SubclusterIndices_2 = max(SubclusterIndices_2H , SubclusterIndices_2D);

PD = fractions;
nP = length(PD);
PH = 1 - PD;
Signals = ones(code_maxSize,timepoints^dimensionality,nP);

AuxiliarySignal_1 = ones([size(Coherences_1H),nP]).*Coherences_1H;
AuxiliarySignal_2 = ones([size(Coherences_2H),nP]).*Coherences_2H;


%--------------------------------------------------------------------------
% 1-Clusters 
%--------------------------------------------------------------------------
cluster_order = 1;

if order<cluster_order, return; end

for iCluster=1:numberClusters(cluster_order)
  
  % Update cluster signals.
  for iorder=cluster_order:maxSize
    
    % The isotope probability is used to determine the likelihood a
    % particular cluster will contain the assumed spins.
    % If the cluster does not contain all spinful nuclei, 
    % the auxiliarary signal will be unity.
    
    
    this_cluster = Clusters(iCluster,1:cluster_order,cluster_order);
    if this_cluster(1)==0
      continue;
    end
    auxH_ = Coherences_1H(iCluster,:);
    auxD_ = Coherences_1D(iCluster,:);
    for ifrac = 1:nP
      AuxiliarySignal_1(iCluster,:,ifrac) = PH(ifrac).*auxH_ + PD(ifrac).*auxD_;
    end
    % v(t) = v'(t)*(  (1-p)  + p*v''(t) ) = v'(t)*(  1  + p*(  v''(t) -1 )  ).
    Signals(iorder,:,:) = Signals(iorder,:,:).*AuxiliarySignal_1(iCluster,:,:);
  end
  
end

%--------------------------------------------------------------------------
% 2-Clusters 
%--------------------------------------------------------------------------
cluster_order = 2;

if order<cluster_order, return; end

for iCluster=1:numberClusters(cluster_order)
  
  auxH_ = Coherences_2H(iCluster,:);
  auxD_ = Coherences_2D(iCluster,:);
  
  a1_ = SubclusterIndices_2(1,1,iCluster);
  b1_ = SubclusterIndices_2(2,1,iCluster);
  
  auxH1a_ = Coherences_1H(a1_,:);
  auxH1b_ = Coherences_1H(b1_,:);
  
  auxD1a_ = Coherences_1D(a1_,:);
  auxD1b_ = Coherences_1D(b1_,:);
  
  this_cluster = Clusters(iCluster,1:cluster_order,cluster_order);
  isSameMolecule = MoleculeID(this_cluster(1)) == MoleculeID(this_cluster(2));
  if isSameMolecule && all(~Exchangable( this_cluster  ))
    for ifrac = 1:nP
      AuxiliarySignal_2(iCluster,:,ifrac) = PH(ifrac).*auxH_ + PD(ifrac).*auxD_ ;
    end
  else
    for ifrac = 1:nP 
      AuxiliarySignal_2(iCluster,:,ifrac) = PH(ifrac)^2.*auxH_ + PD(ifrac)^2.*auxD_ + ...
        PH(ifrac)*PD(ifrac).*(auxH1a_.*auxD1b_ + auxD1a_.*auxH1b_);
    end
  end
  
  % Remove 1-cluster correlations
  for subcluster_index = SubclusterIndices_2(:,1,iCluster)'
    
   if subcluster_index <= 0
    continue;
   end

    AuxiliarySignal_2(iCluster,:,:) = AuxiliarySignal_2(iCluster,:,:)./AuxiliarySignal_1(subcluster_index,:,:);

  end
  
  
  % Update cluster signals.
  for iorder=cluster_order:maxSize
    this_cluster = Clusters(iCluster,1:cluster_order,cluster_order);
    if this_cluster(1)==0
      continue;
    end
    Signals(iorder,:,:) = Signals(iorder,:,:).*AuxiliarySignal_2(iCluster,:,:);
  end
  
end

end