function uncertainty = getClusterError(matfile, ...
  Nuclei, System, nOrientations, Clusters, AuxiliarySignal,Method, gridWeight, Order_n_SignalMean, Order_n_Signals)


% Load data file.
if nargin == 1
  
  disp(matfile);
  indata = load(matfile);
  
  Nuclei             = indata.Nuclei;
  System             = indata.System;
  nOrientations      = indata.nOrientations;
  Clusters           = indata.Clusters;
  AuxiliarySignal    = indata.AuxiliarySignal;
  Method             = indata.Method;
  gridWeight         = indata.gridWeight;
  Order_n_SignalMean = indata.Order_n_SignalMean;
  Order_n_Signals    = indata.Order_n_Signals;
end

% number of bath spin
numberClusters = Nuclei.numberClusters;

% number of time points
nt = System.timepoints;

% matrix of nucleus-nucleus deistances 
DistanceMatrix = Nuclei.DistanceMatrix;

% CCE order
Method_order = Method.order;

% CCE lower bound
Method_order_lower_bound = Method.order_lower_bound; 


% -------------------------------------------------------------------------
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Get Cartesian coordinates.
Coordinates.xyz = Nuclei.Coordinates;

% Get radial distances,
Coordinates.r = vecnorm(Nuclei.Coordinates,2,2);

% Get list of n, from n-CCE calculated.
orderrange = Method_order_lower_bound:Method_order;

% Determine minimum system size needed to hold each cluster.
ClusterGeo = getClusterGeoStats(Clusters,Coordinates, DistanceMatrix, numberClusters,Method_order); 


% -------------------------------------------------------------------------
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Initialize martrices.
Nuclear_Dipole = abs(Nuclei.Statistics{1}.Nuclear_Dipole_perpendicular); % = matrix(N);
bAmax = Nuclei.Statistics{1}.bAmax; % = matrix(N);
Adjacency = Nuclei.Adjacency; % = matrix(N);

% ENUM
MIN = 1; MAX = 2; EDGE = 3; CRIT = 4;
ClusterH = cell(1,Method_order);
ClusterOriH = cell(nOrientations,Method_order);

isize = 1;
ClusterH{isize}.ENUM = ['SELF']; 

ClustersToSearch = cell(1,Method_order);

for isize = 2:Method_order
  
  
  ClusterH{isize} = getClusterH_Stats(Nuclear_Dipole,bAmax,Adjacency,Coordinates,Clusters,ClusterGeo,numberClusters,isize);

  b_min = min(ClusterH{isize}.Nuclear_Dipole(:,4));
  is_b_min = ClusterH{isize}.Nuclear_Dipole(:,4) == b_min;
  b_thr = max(ClusterH{isize}.Nuclear_Dipole(is_b_min,4));
  b_thrIndex = find(ClusterH{isize}.Nuclear_Dipole(:,4) == b_thr);
  
  bAmax_min = min(ClusterH{isize}.bAmax(:,4));
  is_bAmax_min = ClusterH{isize}.bAmax(:,4) == bAmax_min;
  bAmax_thr = max(ClusterH{isize}.bAmax(is_bAmax_min,4));
  bAmax_thrIndex = find(ClusterH{isize}.bAmax(:,4) == bAmax_thr);
  
  
  ClustersToSearch{isize} = [b_thrIndex(1),bAmax_thrIndex(1)]; 
  
  
end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% -------------------------------------------------------------------------


sansClusterV = cell(1,Method_order);

% -------------------------------------------------------------------------

  % Loop though all cluster sizes available.
  for isize = orderrange
    
    nC_ = numel(ClustersToSearch{isize});
    sansClusterV{isize} = zeros(nC_,nt,nOrientations);
    
    for ii = 1:nC_
      icluster = ClustersToSearch{isize}(ii); 
    
      SuperClusterIndices = findSuperClusters(Clusters,isize,icluster);
      
      % Loop through orientations used in powder averaging.
      for iOri = 1:nOrientations
      
        sansClusterV{isize}(ii,:,iOri) =  Order_n_Signals{iOri}{isize};       
        
        for jsize = isize:Method_order
          superclusters = SuperClusterIndices(:,jsize);
          superclusters(superclusters==0) = [];
          v_ = prod(AuxiliarySignal{iOri}{jsize}(superclusters,:),1);
          sansClusterV{isize}(ii,:,iOri) = sansClusterV{isize}(ii,:,iOri)./v_;
        end
        
      end
      
    end
    
  end
  

% -------------------------------------------------------------------------
% Do powder averaging.
sansClusterV_powder = cell(1,Method_order);

sansClusterRMSD_powder = cell(1,Method_order);

for isize = orderrange
  
  nC_ = numel(ClustersToSearch{isize});
  sansClusterV_powder{isize} = sum(sansClusterV{isize},3);
  
  sansClusterRMSD_powder{isize} = zeros(1,nC_);
  
  sansClusterRMSD_powder{isize} = sqrt(mean(abs(sansClusterV_powder{isize}-Order_n_SignalMean{isize}).^2,2));
  
end

% -------------------------------------------------------------------------
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

numberPossibleClusters = NchooseK(numberClusters(1),1:Method_order );

nk = numberPossibleClusters - numberClusters;

v_notC = cell(1,Method_order);
err_C = cell(1,Method_order);
err_max = cell(1,Method_order);
err_unitPseudoGrad = cell(1,Method_order);
err_norm = ones(1,Method_order);
for isize = 2:Method_order

  v_notC{isize} = sum(sansClusterV{isize}(:,:,:),3); 
  err_C{isize} = sqrt(cumsum( abs(Order_n_SignalMean{isize}-v_notC{isize}).^2./(1:nt)) ).*nk(isize);
 
  err_max{isize} = max(err_C{isize}')';
  err_norm(isize) = norm(err_max{isize});
  err_unitPseudoGrad{isize} = err_max{isize}/err_norm(isize); 
end 

uncertainty.v_notC = v_notC;
uncertainty.err_C = err_C;
uncertainty.err_max = err_max;
uncertainty.err_unitPseudoGrad = err_unitPseudoGrad;
uncertainty.err_norm = err_norm;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% -------------------------------------------------------------------------

end



function ClusterH = getClusterH_Stats(Nuclear_Dipole,bAmax,Adjacency,Coordinates,Clusters,ClusterGeo,numberClusters,isize)
 
  % This is a trimmed version of getClusterHStats().
  
  ClusterH.ENUM = {'MIN', 'MAX', 'EDGE' ,'CRIT'};
  
  E_ = 1-eye(isize)>0;
  E_ = E_(:);
  
  nC_ = numberClusters(isize);
  ClusterH.Nuclear_Dipole = zeros(nC_,4);
  ClusterH.bAmax = zeros(nC_,4);
  
  
  for icluster = 1:nC_
    cluster_ = Clusters{isize}(icluster,:);
    adj_ = Adjacency(cluster_,cluster_) > 0;
    adj_ = adj_(:);
    adj_ = adj_(E_);
    
    prox_ = cluster_(find(Coordinates.r(cluster_) == ClusterGeo{isize}.cluster_proximity(icluster)));
    dist_ = cluster_(find(Coordinates.r(cluster_) == ClusterGeo{isize}.cluster_distance(icluster)));
    
    dd_ = Nuclear_Dipole(cluster_,cluster_);
    dd_ = dd_(:);
    dd_ = dd_(E_);
    
    bAmax_ = bAmax(cluster_,cluster_);
    bAmax_ = bAmax_(:);
    bAmax_ = bAmax_(E_);
    
    
    if any(adj_)
      crit_dd_ = min(dd_(adj_));
      crit_bAmax_ = min(bAmax_(adj_));
    else
      crit_dd_ = -1;
      crit_bAmax_ = -1;
    end
    
    
    ClusterH.Nuclear_Dipole(icluster,:) = [min(dd_),max(dd_), Nuclear_Dipole(prox_,dist_), crit_dd_];
    ClusterH.bAmax(icluster,:) = [min(bAmax_),max(bAmax_), Nuclear_Dipole(prox_,dist_), crit_bAmax_];
    
  end
end 
