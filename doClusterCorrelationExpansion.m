function [Signal, AuxiliarySignal,order_n_signal] = doClusterCorrelationExpansion(Coherences,Clusters, SubclusterIndices,order,Nuclei)
  % Coherences{clusterSize,iCluster}(iTime)
  % Clusters{clusterSize}(index,element)
  % SubclusterIndices{hostClusterSize,hostClusterIndex}{subclusterSize}(index)
  
  % Find the first non-empty coherence to establish array size.
  index = 1; 
  while isempty(Coherences{index})
    index = index + 1;
  end
  
  % Initialize data.
  AuxiliarySignal = Coherences;
  Signal = ones(size(Coherences{index}));
  order_n_signal{order} = Signal;
  
%   clusterSize =1;
%   for inucleus = Clusters{clusterSize}'
%     if isempty(AuxiliarySignal{clusterSize,inucleus})
%       AuxiliarySignal{clusterSize,inucleus} = 1;
%       continue;
%     end
%     Signal = Signal.*AuxiliarySignal{clusterSize,inucleus};
%   end
%  
%   order_n_signal{clusterSize} = Signal;
  
 % Loop over non-trivial cluster sizes. 
  for clusterSize = 1:order
    
    % Loop over all clusters of size clusterSize.
    iCluster = 0;
    for ic = 1:size(Clusters{clusterSize},1)
      
      if Clusters{clusterSize}(ic,1)==0
        continue;
      end
      iCluster = iCluster + 1;

      if clusterSize > 1 && isempty(SubclusterIndices{clusterSize,iCluster})
        AuxiliarySignal{clusterSize,iCluster} = 1;
        continue;
      end
      
      % Check to make sure the cluster is made from the union of it size 1
      % subclusters.
      if false %clusterSize > 1
        Test = Clusters{clusterSize}(iCluster,:)'==[Clusters{1}(SubclusterIndices{clusterSize,iCluster}{1}(:))];
        if ~all(all(Test))
          error('Cluster is not made from the union of it size 1 subclusters.');
        end
      end
      
      if isempty(AuxiliarySignal{clusterSize,iCluster})
        AuxiliarySignal{clusterSize,iCluster} = 1;
        continue;
      end
      
      % Loop over non-trivial sub-cluster sizes. 
      for iSize = 1:clusterSize-1

        n_subClusters = size(SubclusterIndices{clusterSize,iCluster}{iSize},2);
        
        % Calculate auxiliary signals.
        for ii = 1:n_subClusters
          AuxiliarySignal{clusterSize,iCluster} =...
            AuxiliarySignal{clusterSize,iCluster}./AuxiliarySignal{  iSize, SubclusterIndices{clusterSize,iCluster}{iSize}(ii)  };
        end
        
      end
      
      % Build decoherence signal.
      
      isotopeProbability = prod(Nuclei.Abundance(Clusters{clusterSize}(iCluster,:)));
      v_ = 1 + isotopeProbability*(AuxiliarySignal{clusterSize,iCluster}-1);
      
      Signal = Signal.*v_;
      
    end
    
    % Save the signal for each cluster size.
    order_n_signal{clusterSize} = Signal;
  end
  

end