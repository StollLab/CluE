function [Signals, ...
          AuxiliarySignal_1,AuxiliarySignal_2, ...
          AuxiliarySignal_3,AuxiliarySignal_4,...
          AuxiliarySignal_5,AuxiliarySignal_6] = ...
  doGeneralClusterExpansion_gpu(...
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

% ENUM 
CCE = 1; CE2 = 2; SUM = 3;
opID = SUM;
switch opID
  case CCE
    opTerm = 0;
    a0 = 1;
    v0 = 1;
  case CE2
    opTerm = -1;
    a0 = 0;
    v0 = 1;
  case SUM 
    opTerm = 0;
    a0 = 0;
    v0 = 1;
end
% Initialize data.
maxSize = 6;
AuxiliarySignal_1 = Coherences_1;
AuxiliarySignal_2 = Coherences_2;
AuxiliarySignal_3 = Coherences_3;
AuxiliarySignal_4 = Coherences_4;
AuxiliarySignal_5 = Coherences_5;
AuxiliarySignal_6 = Coherences_6;

Signals = v0*ones(maxSize,timepoints^dimensionality);

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
      continue
    end
    isotopeProbability =prod(Nuclei_Abundance( this_cluster ));
    
    %For CCE, v(t) = v'(t)*(  (1-p)  + p*v''(t) ) = v'(t)*(  1  + p*(  v''(t) -1 )  ).
    AuxiliarySignal_1(iCluster,:) = inverseOperator(opID,AuxiliarySignal_1(iCluster,:), a0);
    
     
  end
  
  AuxiliarySignal_1(iCluster,:) = AuxiliarySignal_1(iCluster,:) + opTerm;
  
  for iorder=cluster_order:maxSize
    Signals(iorder,:) = Operator(opID,Signals(iorder,:), AuxiliarySignal_1(iCluster,:),isotopeProbability);
  end
end

%--------------------------------------------------------------------------
% n-Clusters 
%--------------------------------------------------------------------------
if order < 2, return; end

for cluster_order = 2:order
  
switch cluster_order
  case 2
    SubclusterIndices = SubclusterIndices_2;
    AuxSig_ = AuxiliarySignal_2;
  case 3
    SubclusterIndices = SubclusterIndices_3;
    AuxSig_ = AuxiliarySignal_3;
  case 4    
    SubclusterIndices = SubclusterIndices_4;
    AuxSig_ = AuxiliarySignal_4;
  case 5    
    SubclusterIndices = SubclusterIndices_5;
    AuxSig_ = AuxiliarySignal_5;
  case 6    
    SubclusterIndices = SubclusterIndices_6;
    AuxSig_ = AuxiliarySignal_6;
end



for iCluster=1:numberClusters(cluster_order)
  
  for subClusterSize= 1:cluster_order-1
    
    switch subClusterSize
      case 1
        aux_ = AuxiliarySignal_1;
      case 2
        aux_ = AuxiliarySignal_2;
      case 3
        aux_ = AuxiliarySignal_3;
      case 4
        aux_ = AuxiliarySignal_4;
      case 5
        aux_ = AuxiliarySignal_5;
      case 6
        aux_ = AuxiliarySignal_6;
      otherwise
        warning('Requested cluster size exceeds the limits of this program.');
        continue;
    end
    
    for subcluster_index = SubclusterIndices(:,subClusterSize,iCluster)'
      
      if subcluster_index <= 0
        continue;
      end
      
      AuxSig_(iCluster,:) = inverseOperator(opID,AuxSig_(iCluster,:), aux_(subcluster_index,:));
      
    end
    
    switch cluster_order
      case 2
        AuxiliarySignal_2 = AuxSig_;% + opTerm;
      case 3
        AuxiliarySignal_3 = AuxSig_;% + opTerm;
      case 4
        AuxiliarySignal_4 = AuxSig_;% + opTerm;
      case 5
        AuxiliarySignal_5 = AuxSig_;% + opTerm;
      case 6
        AuxiliarySignal_6 = AuxSig_;% + opTerm;
    end
  end
  
  % Update cluster signals.
  for iorder=cluster_order:maxSize
    this_cluster = Clusters(iCluster,1:cluster_order,cluster_order);
    if this_cluster(1)==0
      continue
    end
    isotopeProbability =prod(Nuclei_Abundance( this_cluster )); 
    switch cluster_order
      case 2
        Signals(iorder,:) = Operator(opID,Signals(iorder,:), AuxiliarySignal_2(iCluster,:),isotopeProbability);
      case 3
        Signals(iorder,:) = Operator(opID,Signals(iorder,:), AuxiliarySignal_3(iCluster,:),isotopeProbability);
      case 4
        Signals(iorder,:) = Operator(opID,Signals(iorder,:), AuxiliarySignal_3(iCluster,:),isotopeProbability);
      case 5
        Signals(iorder,:) = Operator(opID,Signals(iorder,:), AuxiliarySignal_5(iCluster,:),isotopeProbability);
      case 6
        Signals(iorder,:) = Operator(opID,Signals(iorder,:), AuxiliarySignal_6(iCluster,:),isotopeProbability);
    end
  end
end

end
end



function V = Operator(opID,v, aux,x)

% ENUM 
CCE = 1; CE2 = 2; SUM = 3;
switch opID
  case CCE
    V = v.*(1 + x*(aux - 1));
  case CE2
    V = v.*( 1 + x*aux );
  case SUM
    V = v + x*aux;
end
end


function V = inverseOperator(opID,vC, aux)

% ENUM 
CCE = 1;  CE2 = 2; SUM = 3;
switch opID
  case CCE
    V = vC./aux;
  case CE2 
    V = vC./(aux + 1); % still need to subtract 1 at end.
  case SUM 
    V = vC - aux; 
end
end

