function [Signals, ...
          AuxiliarySignal_1,AuxiliarySignal_2, ...
          AuxiliarySignal_3,AuxiliarySignal_4,...
          AuxiliarySignal_5,AuxiliarySignal_6] = ...
  doClusterCorrelationExpansion(...
  Coherences_1,Coherences_2,Coherences_3,Coherences_4,Coherences_5,Coherences_6,Clusters, ...
  SubclusterIndices_2,SubclusterIndices_3,SubclusterIndices_4,SubclusterIndices_5,SubclusterIndices_6,...
  timepoints,dimensionality, order,numberClusters, Nuclei_Abundance)
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
maxSize = 6;
AuxiliarySignal_1 = Coherences_1;
AuxiliarySignal_2 = Coherences_2;
AuxiliarySignal_3 = Coherences_3;
AuxiliarySignal_4 = Coherences_4;
AuxiliarySignal_5 = Coherences_5;
AuxiliarySignal_6 = Coherences_6;

Signals = ones(maxSize,timepoints^dimensionality);

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
      break;
    end
    isotopeProbability =prod(Nuclei_Abundance( this_cluster ));
    % v(t) = v'(t)*(  (1-p)  + p*v''(t) ) = v'(t)*(  1  + p*(  v''(t) -1 )  ).
    Signals(iorder,:) = Signals(iorder,:).*(1 + isotopeProbability*( AuxiliarySignal_1(iCluster,:) -1) );
  end
  
end

%--------------------------------------------------------------------------
% 2-Clusters 
%--------------------------------------------------------------------------
cluster_order = 2;

if order<cluster_order, return; end

for iCluster=1:numberClusters(cluster_order)
  
  % Remove 1-cluster correlations
  for subcluster_index = SubclusterIndices_2(:,1,iCluster)'
    
   if subcluster_index <= 0
    break;
   end
   
    AuxiliarySignal_2(iCluster,:) = AuxiliarySignal_2(iCluster,:)./AuxiliarySignal_1(subcluster_index,:);

  end
  
  
  % Update cluster signals.
  for iorder=cluster_order:maxSize
    this_cluster = Clusters(iCluster,1:cluster_order,cluster_order);
    if this_cluster(1)==0
      break;
    end
    isotopeProbability =prod(Nuclei_Abundance( this_cluster ));
    Signals(iorder,:) = Signals(iorder,:).*(1 + isotopeProbability*( AuxiliarySignal_2(iCluster,:) -1) );
  end
  
end

%--------------------------------------------------------------------------
% 3-Clusters 
%--------------------------------------------------------------------------
cluster_order = 3;

if order<cluster_order, return; end

for iCluster=1:numberClusters(cluster_order)
  % Remove 1-cluster correlations
  for subcluster_index = SubclusterIndices_3(:,1,iCluster)'
    
   if subcluster_index <= 0
    break;
   end
   
    AuxiliarySignal_3(iCluster,:) = AuxiliarySignal_3(iCluster,:)./AuxiliarySignal_1(subcluster_index,:);
 
  end
  
  % Remove 2-cluster correlations
  for subcluster_index = SubclusterIndices_3(:,2,iCluster)'
    
   if subcluster_index <= 0
     break;
   end
   
    AuxiliarySignal_3(iCluster,:) = AuxiliarySignal_3(iCluster,:)./AuxiliarySignal_2(subcluster_index,:);
 
    
  end
  
  % Update cluster signals.
  for iorder=cluster_order:maxSize
    this_cluster = Clusters(iCluster,1:cluster_order,cluster_order);
    if this_cluster(1)==0
      break;
    end
    isotopeProbability =prod(Nuclei_Abundance( this_cluster )); 
    Signals(iorder,:) = Signals(iorder,:).*(1 + isotopeProbability*( AuxiliarySignal_3(iCluster,:) -1) );
  end
end

%--------------------------------------------------------------------------
% 4-Clusters 
%--------------------------------------------------------------------------
cluster_order = 4;

if order<cluster_order, return; end

for iCluster=1:numberClusters(cluster_order)
    
  % Remove 1-cluster correlations
  for subcluster_index = SubclusterIndices_4(:,1,iCluster)'
        
   if subcluster_index <= 0
    break;
   end
   
   AuxiliarySignal_4(iCluster,:) = AuxiliarySignal_4(iCluster,:)./AuxiliarySignal_1(subcluster_index,:);

  end
  
  % Remove 2-cluster correlations  
  for subcluster_index = SubclusterIndices_4(:,2,iCluster)'
    
   if subcluster_index <= 0
     break;
   end
   
    AuxiliarySignal_4(iCluster,:) = AuxiliarySignal_4(iCluster,:)./AuxiliarySignal_2(subcluster_index,:);

  end
  
  % Remove 3-cluster correlations
  for subcluster_index = SubclusterIndices_4(:,3,iCluster)'
    
   if subcluster_index <= 0
     break;
   end
   
    AuxiliarySignal_4(iCluster,:) = AuxiliarySignal_4(iCluster,:)./AuxiliarySignal_3(subcluster_index,:);

  end

  % Update cluster signals.
  for iorder=cluster_order:maxSize
    this_cluster = Clusters(iCluster,1:cluster_order,cluster_order);
    if this_cluster(1)==0
      break;
    end
    isotopeProbability =prod(Nuclei_Abundance( this_cluster ));
    Signals(iorder,:) = Signals(iorder,:).*(1 + isotopeProbability*( AuxiliarySignal_4(iCluster,:) -1) );
  end
end

%--------------------------------------------------------------------------
% 5-Clusters 
%--------------------------------------------------------------------------
cluster_order = 5;

if order<cluster_order, return; end

for iCluster=1:numberClusters(cluster_order)
    
  % Remove 1-cluster correlations
  for subcluster_index = SubclusterIndices_5(:,1,iCluster)'
        
   if subcluster_index <= 0
    break;
   end
   
   AuxiliarySignal_5(iCluster,:) = AuxiliarySignal_5(iCluster,:)./AuxiliarySignal_1(subcluster_index,:);

  end
  
  % Remove 2-cluster correlations  
  for subcluster_index = SubclusterIndices_5(:,2,iCluster)'
    
   if subcluster_index <= 0
     break;
   end
   
    AuxiliarySignal_5(iCluster,:) = AuxiliarySignal_5(iCluster,:)./AuxiliarySignal_2(subcluster_index,:);

  end
  
  % Remove 3-cluster correlations
  for subcluster_index = SubclusterIndices_5(:,3,iCluster)'
    
   if subcluster_index <= 0
     break;
   end
   
    AuxiliarySignal_5(iCluster,:) = AuxiliarySignal_5(iCluster,:)./AuxiliarySignal_3(subcluster_index,:);

  end
  
  % Remove 4-cluster correlations
  for subcluster_index = SubclusterIndices_5(:,4,iCluster)'
    
   if subcluster_index <= 0
     break;
   end
   
    AuxiliarySignal_5(iCluster,:) = AuxiliarySignal_5(iCluster,:)./AuxiliarySignal_4(subcluster_index,:);

  end
  % Update cluster signals.
  for iorder=cluster_order:maxSize
    this_cluster = Clusters(iCluster,1:cluster_order,cluster_order);
    if this_cluster(1)==0
      break;
    end
    isotopeProbability =prod(Nuclei_Abundance( this_cluster ));
    Signals(iorder,:) = Signals(iorder,:).*(1 + isotopeProbability*( AuxiliarySignal_5(iCluster,:) -1) );
  end
end


%--------------------------------------------------------------------------
% 6-Clusters 
%--------------------------------------------------------------------------
cluster_order = 6;

if order<cluster_order, return; end

for iCluster=1:numberClusters(cluster_order)
    
  % Remove 1-cluster correlations
  for subcluster_index = SubclusterIndices_6(:,1,iCluster)'
        
   if subcluster_index <= 0
    break;
   end
   
   AuxiliarySignal_6(iCluster,:) = AuxiliarySignal_6(iCluster,:)./AuxiliarySignal_1(subcluster_index,:);

  end
  
  % Remove 2-cluster correlations  
  for subcluster_index = SubclusterIndices_6(:,2,iCluster)'
    
   if subcluster_index <= 0
     break;
   end
   
    AuxiliarySignal_6(iCluster,:) = AuxiliarySignal_6(iCluster,:)./AuxiliarySignal_2(subcluster_index,:);

  end
  
  % Remove 3-cluster correlations
  for subcluster_index = SubclusterIndices_6(:,3,iCluster)'
    
   if subcluster_index <= 0
     break;
   end
   
    AuxiliarySignal_6(iCluster,:) = AuxiliarySignal_6(iCluster,:)./AuxiliarySignal_3(subcluster_index,:);

  end
  
  % Remove 4-cluster correlations
  for subcluster_index = SubclusterIndices_6(:,4,iCluster)'
    
   if subcluster_index <= 0
     break;
   end
   
   AuxiliarySignal_6(iCluster,:) = AuxiliarySignal_6(iCluster,:)./AuxiliarySignal_4(subcluster_index,:);

  end
  
  % Remove 5-cluster correlations
  for subcluster_index = SubclusterIndices_6(:,5,iCluster)'
    
    if subcluster_index <= 0
      break;
    end
    
    AuxiliarySignal_6(iCluster,:) = AuxiliarySignal_6(iCluster,:)./AuxiliarySignal_5(subcluster_index,:);
    
  end
  % Update cluster signals.
  for iorder=cluster_order:maxSize
    this_cluster = Clusters(iCluster,1:cluster_order,cluster_order);
    if this_cluster(1)==0
      break;
    end
    isotopeProbability =prod(Nuclei_Abundance( this_cluster ));
    Signals(iorder,:) = Signals(iorder,:).*(1 + isotopeProbability*( AuxiliarySignal_6(iCluster,:) -1) );
  end
end
end