function [partial_signal,Cluster_Statistics,NonCluster_Statistics] = ...
  conserve_memory_gpu_loop( ... %System,Method,
  OutputData,... %Nuclei, ...
  par_cluster,Bundle,iorder,iCore,Division_Cluster_Limit,Division_Cluster_Fraction,Division_Cluster_Increment,shuffle,...
  timepoints,dt, ... % Method_order, ...
  EXPERIMENT,dimensionality,System_full_Sz_Hyperfine,total_time, ...  
  Nuclei_Coordinates, Nuclei_ValidPair, Nuclei_Abundance, Nuclei_Spin, Nuclei_g, NumberStates, ZeemanStates, ...
  Spin2Op1, Spin2Op2, Spin2Op3, Spin2Op4, Spin2Op5, Spin2Op6, ...
  Spin3Op1, Spin3Op2, Spin3Op3, Spin3Op4, Spin3Op5, Spin3Op6, ...
  numberClusters,... %ClusterArray,SubclusterIndices_2,SubclusterIndices_3,SubclusterIndices_4, ...
  ge, magneticField, muB, muN, mu0, hbar, ...
  Method_seed, Nuclei_number,graphCriterion,Method_partialSave,Method_MonteCarlo_Threshold,maxPossibleNumberSubClusters, ...
  Reduced_ClusterArray, SubclusterIndices_2,SubclusterIndices_3,SubclusterIndices_4,...
  SubclusterIndices_5,SubclusterIndices_6, ...
  Method_record_clusters,useHamiltonian)

%--------------------------------------------------------------------------

% set signal to unperturbed value
partial_signal = ones(1,timepoints^dimensionality);
partial_signal0 = partial_signal;

Cluster_Statistics = 0;

NonCluster_Statistics = 0;
% Save to file.
if Method_partialSave
  partial_file = ['partial_', OutputData(1:end-4), '_',num2str(iorder),'cce_', num2str(iCore), '.mat'] ;
  if isfile(partial_file)
    try
      load(partial_file,'parsignal','seed');
      if seed==Method_seed
        partial_signal = parsignal;
        return;
      end
    catch
    end
  end
end

% loop over all the clusters in the bundle
iBundle = Bundle(iCore,1);

while iBundle <= Bundle(iCore,2)
  
  % Do not calculate more than Division_Cluster_Limit valid clusters.
  if Cluster_Statistics>Division_Cluster_Limit
    
    clusterFraction = (NonCluster_Statistics + Cluster_Statistics)/Division_Cluster_Possible_Clusters;
    partial_signal1 = partial_signal.^(1/clusterFraction);
    
    delta_sig = (partial_signal1-partial_signal0);
    rmsd = sqrt(  delta_sig*delta_sig'/System.timepoints  );
        
    if rmsd < Method_MonteCarlo_Threshold(iorder) || clusterFraction > Division_Cluster_Fraction  
      break;
    end
    
    partial_signal0 = partial_signal1;

    Division_Cluster_Limit = Division_Cluster_Limit + Division_Cluster_Increment;
  
  end
  
  % shuffle nuclei
  shuffled_cluster = shuffle(par_cluster);
  
  % check cluster validity  
  isvalid_1 = validateCluster(shuffled_cluster,Nuclei_ValidPair,graphCriterion);
  
  if ~isvalid_1
    
    % update internal records.
    NonCluster_Statistics = NonCluster_Statistics + 1;
    
    % get next cluster
    par_cluster = getNextCluster(par_cluster,Nuclei_number);
    
    % check if there are no more clusters
    if isempty(par_cluster)
      break;
    end
    
    % update location index
    iBundle = par_cluster(1);
    
    % start over
    continue;
  end
  

  if Method_record_clusters
    % Generate text fileof clusters.
    record_cluster(shuffled_cluster, OutputData);
  end
  
  % update internal records.
  Cluster_Statistics = Cluster_Statistics + 1;
    
  % re-index
  new_labels = zeros(1,max( shuffled_cluster  ));
  old_labels = zeros(1,max( shuffled_cluster  ));
  for ii=1:iorder
    new_labels( shuffled_cluster (ii) ) = ii;
    old_labels(new_labels( shuffled_cluster (ii) )) = shuffled_cluster (ii);
  end
  
  % Initialize cluster array.
  Reduced_Clusters = Reduced_ClusterArray;
  
  % Assign labels.
  Reduced_Clusters(Reduced_ClusterArray>0) = old_labels( Reduced_ClusterArray(Reduced_ClusterArray>0) );
  
  % Initialize copy of the subcluster indices.
  Reduced_SubclusterIndices_2 = SubclusterIndices_2;
  Reduced_SubclusterIndices_3 = SubclusterIndices_3;
  Reduced_SubclusterIndices_4 = SubclusterIndices_4;
  Reduced_SubclusterIndices_5 = SubclusterIndices_5;
  Reduced_SubclusterIndices_6 = SubclusterIndices_6;
  
  % Loop through subcluster sizes.
  for jsize = 2:iorder-1
    
    % Loop through subclusters.
    for jCluster = 1:maxPossibleNumberSubClusters(iorder,jsize)
      
      cluster_ = Reduced_Clusters(jCluster,1:jsize,jsize);
      zero_cluster = cluster_(1)==0 || ~validateCluster(cluster_ ,Nuclei_ValidPair,graphCriterion);
      
      if zero_cluster
        
        % Zero disconnected clusters.
        Reduced_Clusters(jCluster,1:jsize,jsize) = 0*cluster_;
        
        % SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
        % the jth cluster of size subCluster_size that is a subcluster of
        % the ith cluster of size clusterSize.
        % In a valid cluster, all indices must be natural numbers
        
        % Zero disconnected subclusters.
        switch jsize
          
          case 2
            Reduced_SubclusterIndices_2(:,:,jCluster) = 0*Reduced_SubclusterIndices_2(:,:,jCluster);
           
            if isempty(Reduced_SubclusterIndices_3)
              continue;
            end
            toZero3 = Reduced_SubclusterIndices_3(:,jsize,:); 
            toZero3(toZero3==jCluster) = 0;
            toZero3(toZero3 >= 1) = 1;
            Reduced_SubclusterIndices_3(:,jsize,:) = toZero3.*Reduced_SubclusterIndices_3(:,jsize,:);
        
            if isempty(Reduced_SubclusterIndices_4)
              continue;
            end
            toZero4 = Reduced_SubclusterIndices_4(:,jsize,:); 
            toZero4(toZero4==jCluster) = 0;
            toZero4(toZero4 >= 1) = 1;
            Reduced_SubclusterIndices_4(:,jsize,:) = toZero4.*Reduced_SubclusterIndices_4(:,jsize,:);
            
          case 3
            Reduced_SubclusterIndices_3(:,:,jCluster) = 0*Reduced_SubclusterIndices_3(:,:,jCluster);
          
            if isempty(Reduced_SubclusterIndices_4)
              continue;
            end
            toZero4 = Reduced_SubclusterIndices_4(:,jsize,:); 
            toZero4(toZero4==jCluster) = 0;
            toZero4(toZero4 >= 1) = 1;
            Reduced_SubclusterIndices_4(:,jsize,:) = toZero4.*Reduced_SubclusterIndices_4(:,jsize,:);
          case 4
            Reduced_SubclusterIndices_4(:,:,jCluster) = 0*Reduced_SubclusterIndices_4(:,:,jCluster);
            
            if isempty(Reduced_SubclusterIndices_5)
              continue;
            end
            
            toZero5 = Reduced_SubclusterIndices_5(:,jsize,:);
            toZero5(toZero5==jCluster) = 0;
            toZero5(toZero5 >= 1) = 1;
            Reduced_SubclusterIndices_5(:,jsize,:) = toZero5.*Reduced_SubclusterIndices_5(:,jsize,:);
            
            if isempty(Reduced_SubclusterIndices_6)
              continue;
            end
            
            toZero6 = Reduced_SubclusterIndices_6(:,jsize,:);
            toZero6(toZero6==jCluster) = 0;
            toZero6(toZero6 >= 1) = 1;
            Reduced_SubclusterIndices_6(:,jsize,:) = toZero6.*Reduced_SubclusterIndices_6(:,jsize,:);
         
          case 5
            Reduced_SubclusterIndices_5(:,:,jCluster) = 0*Reduced_SubclusterIndices_5(:,:,jCluster);
            
            if isempty(Reduced_SubclusterIndices_6)
              continue;
            end
            
            toZero6 = Reduced_SubclusterIndices_6(:,jsize,:);
            toZero6(toZero6==jCluster) = 0;
            toZero6(toZero6 >= 1) = 1;
            Reduced_SubclusterIndices_6(:,jsize,:) = toZero6.*Reduced_SubclusterIndices_6(:,jsize,:);
            
          case 6
            Reduced_SubclusterIndices_6(:,:,jCluster) = 0*Reduced_SubclusterIndices_6(:,:,jCluster);
            
            
        end
      end
      
    end   
  end
 
 
  [~, AuxiliarySignal_1,AuxiliarySignal_2,AuxiliarySignal_3,AuxiliarySignal_4,...
          AuxiliarySignal_5,AuxiliarySignal_6,~] ... 
       = calculateSignal_gpu(...
       timepoints,dt, iorder, EXPERIMENT,dimensionality,System_full_Sz_Hyperfine,total_time, ...
       Nuclei_Coordinates,Nuclei_ValidPair, graphCriterion, Nuclei_Abundance, Nuclei_Spin, Nuclei_g, NumberStates, ZeemanStates, ...
       Spin2Op1, Spin2Op2, Spin2Op3, Spin2Op4, Spin2Op5, Spin2Op6, ...
       Spin3Op1, Spin3Op2, Spin3Op3, Spin3Op4, Spin3Op5, Spin3Op6, ...
       numberClusters,Reduced_Clusters, ...
       Reduced_SubclusterIndices_2,Reduced_SubclusterIndices_3, ...
       Reduced_SubclusterIndices_4, ...
       Reduced_SubclusterIndices_5,Reduced_SubclusterIndices_6, ...
       ge, magneticField, muB, muN, mu0, hbar,useHamiltonian);
  
  % collect cluster contribution to a the node output

  isotopeProbability = prod(Nuclei_Abundance( shuffled_cluster ));
  switch iorder
    case 1
      v_ = 1 + isotopeProbability*(AuxiliarySignal_1(1,:) - 1);
    case 2
      v_ = 1 + isotopeProbability*(AuxiliarySignal_2(1,:) - 1);
    case 3
      v_ = 1 + isotopeProbability*(AuxiliarySignal_3(1,:)- 1);
    case 4
      v_ = 1 + isotopeProbability*(AuxiliarySignal_4(1,:) - 1);
    case 5
      v_ = 1 + isotopeProbability*(AuxiliarySignal_5(1,:)- 1);
    case 6
      v_ = 1 + isotopeProbability*(AuxiliarySignal_6(1,:) - 1);
  end
  
  % partial_signal  = partial_signal.*Cluster_AuxiliarySignal{iorder,1};
  partial_signal  = partial_signal.*v_;
  
  % get next cluster
  par_cluster = getNextCluster(par_cluster,Nuclei_number);
  
  % checkk if next cluster exists
  if isempty(par_cluster)
    break;
  end
  
  % update location index
  iBundle = par_cluster(1);
end


% Save to file.
if Method_partialSave
  partial_file = ['partial_', OutputData(1:end-4), '_',num2str(iorder),'cce_', num2str(iCore), '.mat'] ;
  parsavefile = matfile(partial_file,'writable',true);
  parsavefile.parsignal = partial_signal;
  parsavefile.seed = Method_seed;
end


end

