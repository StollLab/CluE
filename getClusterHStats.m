function ClusterH = getClusterHStats(Hyperfine,DeltaHyperfine,Nuclear_Dipole, ModulationDepth,Frequency_Pair,Adjacency,Coordinates,Clusters,ClusterGeo,numberClusters,isize,TM_powder)
 
  ClusterH.ENUM = {'MIN', 'MAX', 'EDGE' ,'CRIT'};
  
  E_ = 1-eye(isize)>0;
  E_ = E_(:);
  
  nC_ = numberClusters(isize);
  ClusterH.ModulationDepth = zeros(nC_,4);
  ClusterH.Hyperfine = zeros(nC_,2);
  ClusterH.DeltaHyperfine = zeros(nC_,4);
  ClusterH.Nuclear_Dipole = zeros(nC_,4);
  ClusterH.Frequency_Pair = zeros(nC_,4);
  
  for icluster = 1:nC_
    cluster_ = Clusters{isize}(icluster,:);
    adj_ = Adjacency(cluster_,cluster_) > 0;
    adj_ = adj_(:);
    adj_ = adj_(E_);
    
    prox_ = cluster_(find(Coordinates.r(cluster_) == ClusterGeo{isize}.cluster_proximity(icluster)));
    dist_ = cluster_(find(Coordinates.r(cluster_) == ClusterGeo{isize}.cluster_distance(icluster)));
    
    moddepth_ = ModulationDepth(cluster_,cluster_);
    moddepth_ = moddepth_(:);
    moddepth_ = moddepth_(E_);
    
    hyperfine_ = Hyperfine(cluster_);
    
    deltahyperfine_ = DeltaHyperfine(cluster_,cluster_);
    deltahyperfine_ = deltahyperfine_(:);
    deltahyperfine_ = deltahyperfine_(E_);
    
    dd_ = Nuclear_Dipole(cluster_,cluster_);
    dd_ = dd_(:);
    dd_ = dd_(E_);
    
    freq_ = Frequency_Pair(cluster_,cluster_);
    freq_ = freq_(:);
    freq_ = freq_(E_);
    
    
    if any(adj_)
      crit_moddepth_ = min(moddepth_(adj_));
      crit_deltahyperfine_ = min(deltahyperfine_(adj_));
      crit_dd_ = min(dd_(adj_));
      crit_freq_ = min(freq_(adj_));
    else
      crit_moddepth_ = -1;
      crit_deltahyperfine_ = -1;
      crit_dd_ = -1;
      crit_freq_ = -1;
    end
    
    
    ClusterH.ModulationDepth(icluster,:) = [min(moddepth_),max(moddepth_), ModulationDepth(prox_,dist_), crit_moddepth_];
    ClusterH.Hyperfine(icluster,:) = [min(abs(hyperfine_)),max(hyperfine_)];
    ClusterH.DeltaHyperfine(icluster,:) = [min(deltahyperfine_),max(deltahyperfine_), DeltaHyperfine(prox_,dist_), crit_deltahyperfine_];
    ClusterH.Nuclear_Dipole(icluster,:) = [min(dd_),max(dd_), Nuclear_Dipole(prox_,dist_), crit_dd_];
    ClusterH.Frequency_Pair(icluster,:) = [min(freq_),max(freq_), Frequency_Pair(prox_,dist_), crit_freq_];
  end
  
  ClusterH.DeltaHyperfine_over_Nuclear_Dipole = ClusterH.DeltaHyperfine./ClusterH.Nuclear_Dipole;
  ClusterH.GaussianRMSD = getGaussianRMSD(ClusterH.ModulationDepth,ClusterH.Frequency_Pair,TM_powder);
end 