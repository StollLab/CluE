function ana = getClusterContributions(matfile, ...
  Nuclei, System, nOrientations, Clusters, Signals, AuxiliarySignal,Method, t, gridWeight, TM_powder, Input, SignalMean, Order_n_SignalMean, Order_n_Signals)


% Load data file.
if nargin == 1
  
  disp(matfile);
  indata = load(matfile);
  
  Nuclei             = indata.Nuclei;
  System             = indata.System;
  nOrientations      = indata.nOrientations;
  Clusters           = indata.Clusters;
  Signals            = indata.Signals;
  AuxiliarySignal    = indata.AuxiliarySignal;
  Method             = indata.Method;
  t                  = indata.experiment_time;
  gridWeight         = indata.gridWeight;
  TM_powder          = indata.TM_powder;
  Input              = indata.Input;
  SignalMean         = indata.SignalMean;
  Order_n_SignalMean = indata.Order_n_SignalMean;
  Order_n_Signals    = indata.Order_n_Signals;
end

% number of bath spin
nSpins = Nuclei.number;
numberClusters = Nuclei.numberClusters;
nt = System.timepoints;
DistanceMatrix = Nuclei.DistanceMatrix;
Method_order = Method.order;
Method_order_lower_bound = Method.order_lower_bound;
% CluserArray(iCluster,:,clusterSize) = nuclear indices.
maxClusterSize = min(12,Method_order);
maxNumberClusters = max(numberClusters(1:maxClusterSize));
ClusterArray = zeros(maxNumberClusters,maxClusterSize,maxClusterSize);
for isize = 1:Method_order  
  for ii = 1:Nuclei.numberClusters(isize)
    ClusterArray(ii,1:isize,isize) = Clusters{isize}(ii,:);
  end
end

% -------------------------------------------------------------------------
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Get spherical coordinates.

Coordinates.xyz = Nuclei.Coordinates;
Coordinates.r = vecnorm(Nuclei.Coordinates,2,2);
Coordinates.cos_theta = Nuclei.Coordinates(:,3)./Coordinates.r;
Coordinates.theta = acos(Coordinates.cos_theta); 
Coordinates.rho = vecnorm(Nuclei.Coordinates(:,1:2),2,2);
Coordinates.cos_phi= Nuclei.Coordinates(:,1)./Coordinates.rho;
Coordinates.phi = acos(Coordinates.cos_phi) + (1 - sign(Nuclei.Coordinates(:,2)))*pi/2;


% Get list of n, from n-CCE calculated.
orderrange = Method_order_lower_bound:Method_order;


% Determine minimum system size needed to hold each cluster.
ClusterGeo = cell(1,Method_order);
sansClusterV = cell(1,Method_order);

for isize = Method_order:-1:1
  
  nC_ = numberClusters(isize);
  C_ = 1:isize;
  ClusterGeo{isize}.number = nC_;
  ClusterGeo{isize}.cluster_indices = 1:nC_;
  ClusterGeo{isize}.cluster_distance = max(Coordinates.r(Clusters{isize}),[],2); 
  

  ClusterGeo{isize}.cluster_proximity = min(Coordinates.r(Clusters{isize}),[],2);
  ClusterGeo{isize}.cluster_distance_indices = diag(Clusters{isize}(ClusterGeo{isize}.cluster_indices', sum(  (Coordinates.r(Clusters{isize}) == ClusterGeo{isize}.cluster_distance).*C_,2  )));
  ClusterGeo{isize}.cluster_proximity_indices = diag(Clusters{isize}(ClusterGeo{isize}.cluster_indices', sum(  (Coordinates.r(Clusters{isize}) == ClusterGeo{isize}.cluster_proximity).*C_ ,2 )));

  r1_ = ClusterGeo{isize}.cluster_distance;
  ind1_ = ClusterGeo{isize}.cluster_distance_indices;
  
  r2_ = ClusterGeo{isize}.cluster_proximity;
  ind2_ = ClusterGeo{isize}.cluster_proximity_indices;
  if isize==1
    continue;
  end
  
  ClusterGeo{isize}.cluster_separation = diag(DistanceMatrix(ind1_, ind2_));

  R_ = ClusterGeo{isize}.cluster_separation;
  % R^2 = r1^2 + r2^2 -2*r1*r2*cos(theta).
  % 2*r1*r2*cos(theta) = r1^2 + r2^2 -R^2.
  % cos(theta) = (r1^2 + r2^2 -R^2)/(2*r1*r2).
  ClusterGeo{isize}.cluster_cos = (r1_.^2 + r2_.^2 - R_.^2)./(2*r1_.*r2_);

  sansClusterV{isize} = zeros(nC_,nt,nOrientations);
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% -------------------------------------------------------------------------
%{
% Set up parallel computing.

% remove current pool if it exists.
delete(gcp('nocreate'));

% determine number of cores available
numCores = feature('numcores');

% create parallel pool
pool = parpool(numCores);
%}

% -------------------------------------------------------------------------

  % Loop though all cluster sizes available.
  for isize = orderrange
    
    for icluster = 1:numberClusters(isize)
    
      SuperClusterIndices = findSuperClusters(ClusterArray,isize,icluster);
      
      % Loop through orientations used in powder averaging.
      for iOri = 1:nOrientations
      
        sansClusterV{isize}(icluster,:,iOri) = Signals{iOri};       
        
        for jsize = isize:Method_order
          superclusters = SuperClusterIndices(:,jsize);
          superclusters(superclusters==0) = [];
          v_ = prod(AuxiliarySignal{iOri}{jsize}(superclusters,:),1);
          sansClusterV{isize}(icluster,:,iOri) = sansClusterV{isize}(icluster,:,iOri)./v_;
        end
        
      end
      
    end
    
  end
  

% -------------------------------------------------------------------------
% Do powder averaging.
sansClusterV_powder = cell(1,Method_order);

sansClusterTM_powder = cell(1,Method_order);
sansClusterDeltaTM_powder = cell(1,Method_order);
sansClusterRMSD_powder = cell(1,Method_order);
sansClusterRMSDsorted_powder = cell(1,Method_order);
sansClusterRMSDsortOrder_powder = cell(1,Method_order);
for isize = orderrange
sansClusterV_powder{isize} = sum(sansClusterV{isize},3);

sansClusterTM_powder{isize} = zeros(1,numberClusters(isize));
sansClusterDeltaTM_powder{isize} = zeros(1,numberClusters(isize));
sansClusterRMSD_powder{isize} = zeros(1,numberClusters(isize)); 
sansClusterRMSDsorted_powder{isize} = zeros(1,numberClusters(isize));
sansClusterRMSDsortOrder_powder{isize} = zeros(1,numberClusters(isize));


for icluster = 1:numberClusters(isize)
  sansClusterTM_powder{isize}(icluster) = getTM(t,sansClusterV_powder{isize}(icluster,:));
end

sansClusterRMSD_powder{isize} = sqrt(mean(abs(sansClusterV_powder{isize}-SignalMean).^2,2));

[sansClusterRMSDsorted_powder{isize},sansClusterRMSDsortOrder_powder{isize}] = sort(sansClusterRMSD_powder{isize},'descend');
sansClusterDeltaTM_powder{isize} = sansClusterTM_powder{isize} - TM_powder;
end

% -------------------------------------------------------------------------
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

RMSDclu = cell(1,Method_order);
Vclu = cell(1,Method_order); 
for isize = 1:Method_order
  RMSDclu{isize} = zeros(1,numberClusters(isize)); 
  Vclu{isize}.found_10pc = false;
  Vclu{isize}.found_5pc = false;
  Vclu{isize}.found_2pc = false;
  Vclu{isize}.found_1pc = false;
  Vclu{isize}.found_5pm = false;
  Vclu{isize}.found_2pm = false;
  Vclu{isize}.found_1pm = false;
  
  for icluster = 1:numberClusters(isize)
    v_ = zeros(1,nt);
    theseClusters = sansClusterRMSDsortOrder_powder{isize}(1:icluster);
   
    for iOri = 1:nOrientations
      
      if isize ==1
        v0_ = ones(1,nt);
      else
        v0_ = Order_n_Signals{iOri}{isize -1};
        v0_ = v0_./v0_(1);
      end
      
      v_ = v_ + gridWeight(iOri)*prod(AuxiliarySignal{iOri}{isize}(theseClusters,:),1).*v0_;
      
    end
    RMSDclu{isize}(icluster)  = sqrt(mean(abs(v_ - Order_n_SignalMean{isize}).^2,2));
    
    if (~Vclu{isize}.found_10pc) && RMSDclu{isize}(icluster) <= 10e-2
      Vclu{isize}.found_10pc = true;
      Vclu{isize}.thr_10pc = v_;
      Vclu{isize}.icluster_5pc = icluster;
    end
    if (~Vclu{isize}.found_5pc) && RMSDclu{isize}(icluster) <= 5e-2
      Vclu{isize}.found_5pc = true;
      Vclu{isize}.thr_5pc = v_;
      Vclu{isize}.icluster_5pc = icluster;
    end
    if (~Vclu{isize}.found_2pc) && RMSDclu{isize}(icluster) <= 2e-2
      Vclu{isize}.found_2pc = true;
      Vclu{isize}.thr_2pc = v_;
      Vclu{isize}.icluster_2pc = icluster;
    end
    if (~Vclu{isize}.found_1pc) && RMSDclu{isize}(icluster) <= 1e-2
      Vclu{isize}.found_1pc = true;
      Vclu{isize}.thr_1pc = v_;
      Vclu{isize}.icluster_1pc = icluster;
    end
    if (~Vclu{isize}.found_5pm) && RMSDclu{isize}(icluster) <= 5e-3
      Vclu{isize}.found_5pm = true;
      Vclu{isize}.thr_5pm = v_;
      Vclu{isize}.icluster_5pm = icluster;
    end
    if (~Vclu{isize}.found_2pm) && RMSDclu{isize}(icluster) <= 2e-3
      Vclu{isize}.found_2pm = true;
      Vclu{isize}.thr_2pm = v_;
      Vclu{isize}.icluster_2pm = icluster;
    end
    if (~Vclu{isize}.found_1pm) && RMSDclu{isize}(icluster) <= 1e-3
      Vclu{isize}.found_1pm = true;
      Vclu{isize}.thr_1pm = v_;
      Vclu{isize}.icluster_1pm = icluster;
    end
  end
end
ana.RMSDclu = RMSDclu;
ana.Vclu = Vclu;
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ModulationDepth = Nuclei.ModulationDepth; % = matrix(N);
Hyperfine = abs(Nuclei.Hyperfine); % = matrix(N,1);
Nuclear_Dipole = abs(Nuclei.Nuclear_Dipole); % = matrix(N);
Frequency_Pair = Nuclei.Frequency_Pair; % = matrix(N);
DeltaHyperfine = abs(Nuclei.DeltaHyperfine); % = matrix(N);
Adjacency = Nuclei.ValidPair; % = matrix(N);

% ENUM
MIN = 1; MAX = 2; EDGE = 3; CRIT = 4;
ClusterH = cell(1,Method_order);

isize = 1;
ClusterH{isize}.ENUM = ['SELF']; 
ClusterH{isize}.Hyperfine  = Hyperfine;


for isize = 2:Method_order
  ClusterH{isize}.ENUM = {'MIN', 'MAX', 'EDGE' ,'CRIT'};
  
  E_ = 1-eye(isize)>0;
  E_ = E_(:);
  
  nC_ = numberClusters(isize);
  ClusterH{isize}.ModulationDepth = zeros(nC_,4);
  ClusterH{isize}.Hyperfine = zeros(nC_,2);
  ClusterH{isize}.DeltaHyperfine = zeros(nC_,4);
  ClusterH{isize}.Nuclear_Dipole = zeros(nC_,4);
  ClusterH{isize}.Frequency_Pair = zeros(nC_,4);
  
  for icluster = 1:nC_
    cluster_ = ClusterArray(icluster,1:isize,isize);
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
    
    
    ClusterH{isize}.ModulationDepth(icluster,:) = [min(moddepth_),max(moddepth_), ModulationDepth(prox_,dist_), min(moddepth_(adj_))];
    ClusterH{isize}.Hyperfine(icluster,:) = [min(abs(hyperfine_)),max(hyperfine_)];
    ClusterH{isize}.DeltaHyperfine(icluster,:) = [min(deltahyperfine_),max(deltahyperfine_), DeltaHyperfine(prox_,dist_), min(deltahyperfine_(adj_))];
    ClusterH{isize}.Nuclear_Dipole(icluster,:) = [min(dd_),max(dd_), Nuclear_Dipole(prox_,dist_), min(dd_(adj_))];
    ClusterH{isize}.Frequency_Pair(icluster,:) = [min(freq_),max(freq_), Frequency_Pair(prox_,dist_), min(freq_(adj_))];
  end
  
  ClusterH{isize}.DeltaHyperfine_over_Nuclear_Dipole = ClusterH{isize}.DeltaHyperfine./ClusterH{isize}.Nuclear_Dipole;
  
end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% -------------------------------------------------------------------------


% Pack output variable together. ------------------------------------------

% total evolution time
sim.time = t;
% set of nth order CCE sims
sim.Order_n_SignalMean = Order_n_SignalMean;
% original simulation input
sim.Input = Input;
% original simulation TM of the full signal
sim.TM = TM_powder;
% list of atom identities
sim.Type = Nuclei.Type;
% full signal
sim.V = SignalMean;

ana.sim = sim;

% electron-nucleus coordinates
ana.Coordinates = Coordinates;
ana.Clusters = Clusters;
ana.ClusterArray = ClusterArray;
ana.Adjacency = Adjacency;
ana.DistanceMatrix = DistanceMatrix;

% sansCluster
sansCluster.sansClusterV = sansClusterV;
sansCluster.sansClusterTM_powder = sansClusterTM_powder;
sansCluster.sansClusterDeltaTM_powder = sansClusterDeltaTM_powder;
sansCluster.sansClusterRMSD_powder = sansClusterRMSD_powder;
ana.sansCluster = sansCluster;

ana.ClusterGeo  = ClusterGeo;
 
ana.ClusterH = ClusterH;

% Save.
save( [matfile(1:end-4),'_analysis.mat'] ,'ana');

% Close parallel pool.
% delete(pool);
end
