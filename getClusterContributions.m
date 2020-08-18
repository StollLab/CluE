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
numberClusters = Nuclei.numberClusters;
nt = System.timepoints;
DistanceMatrix = Nuclei.DistanceMatrix;
Method_order = Method.order;
Method_order_lower_bound = Method.order_lower_bound; 


% -------------------------------------------------------------------------
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Get coordinates.
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
ClusterGeo = getClusterGeoStats(Clusters,Coordinates, DistanceMatrix, numberClusters,Method_order); 
sansClusterV = cell(1,Method_order);

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
    
    nC_ = numberClusters(isize);
    sansClusterV{isize} = zeros(nC_,nt,nOrientations);
    
    for icluster = 1:numberClusters(isize)
    
      SuperClusterIndices = findSuperClusters(Clusters,isize,icluster);
      
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
numberPossibleClusters = NchooseK(numberClusters(1),1:Method_order );

nk = numberPossibleClusters - numberClusters;
for isize = 1:Method_order
  minIndexClusters = sansClusterRMSDsortOrder_powder{isize}(end);
  Vclu{isize}.error = nk(isize)*sqrt(cumsum(abs(sansClusterV_powder{isize}(minIndexClusters,:)-Order_n_SignalMean{isize}).^2)./(1:nt));
  Vclu{isize}.lowerBound = Order_n_SignalMean{isize} - 2*numberClusters(isize)*sqrt(cumsum(abs(sansClusterV_powder{isize}(minIndexClusters,:)-Order_n_SignalMean{isize}).^2)./(1:nt));
  
  RMSDclu{isize} = zeros(1,numberClusters(isize)); 
  Vclu{isize}.found_10pc   = false;
  Vclu{isize}.found_5pc    = false;
  Vclu{isize}.found_2pc    = false;
  Vclu{isize}.found_1pc    = false;
  Vclu{isize}.found_5pm    = false;
  Vclu{isize}.found_2pm    = false;
  Vclu{isize}.found_1pm    = false;
  Vclu{isize}.found_500ppm = false;
  Vclu{isize}.found_200ppm = false;
  Vclu{isize}.found_100ppm = false;
  Vclu{isize}.found_50ppm  = false;
  Vclu{isize}.found_20ppm  = false;
  Vclu{isize}.found_10ppm  = false;
  Vclu{isize}.found_5ppm   = false;
  Vclu{isize}.found_2ppm   = false;
  Vclu{isize}.found_1ppm   = false;
  
  v_ori_ = ones(nOrientations,nt);
  for iOri = 1:nOrientations
    if isize ==1
      v0_ = ones(1,nt);
    else
      v0_ = Order_n_Signals{iOri}{isize -1};
      v0_ = v0_./v0_(1);
    end
    v_ori_(iOri,:) = gridWeight(iOri)*v0_;
  end
  v_ = ones(1,nt);
  v__ = ones(1,nt);
  for icluster = 1:numberClusters(isize)
    indexClusters = sansClusterRMSDsortOrder_powder{isize}(icluster);
   
    for iOri = 1:nOrientations
      v_ori_(iOri,:) = v_ori_(iOri,:).*prod(AuxiliarySignal{iOri}{isize}(indexClusters,:),1); 
    end
    v__ = v_;
    v_ = sum(v_ori_,1);
    
    RMSDclu{isize}(icluster)  = sqrt(mean(abs(v_ - Order_n_SignalMean{isize}).^2,2));
    
    
    if (~Vclu{isize}.found_10pc) && RMSDclu{isize}(icluster) <= 10e-2
      Vclu{isize}.found_10pc = true;
      Vclu{isize}.thr_10pc = v_;
      Vclu{isize}.icluster_5pc = icluster;
      
      Vclu{isize}.error_10pc = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt)); 
 
    end
    if (~Vclu{isize}.found_5pc) && RMSDclu{isize}(icluster) <= 5e-2
      Vclu{isize}.found_5pc = true;
      Vclu{isize}.thr_5pc = v_;
      Vclu{isize}.icluster_5pc = icluster;
      
      Vclu{isize}.error_5pc = (numberPossibleClusters(isize) - icluster)*sqrt(abs(v__-v_).^2./(1:nt)); 
      
    end
    if (~Vclu{isize}.found_2pc) && RMSDclu{isize}(icluster) <= 2e-2
      Vclu{isize}.found_2pc = true;
      Vclu{isize}.thr_2pc = v_;
      Vclu{isize}.icluster_2pc = icluster;
       
      Vclu{isize}.error_2pc = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt)); 
      
    end
    if (~Vclu{isize}.found_1pc) && RMSDclu{isize}(icluster) <= 1e-2
      Vclu{isize}.found_1pc = true;
      Vclu{isize}.thr_1pc = v_;
      Vclu{isize}.icluster_1pc = icluster;
      
      Vclu{isize}.error_1pc = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt)); 
    end
    if (~Vclu{isize}.found_5pm) && RMSDclu{isize}(icluster) <= 5e-3
      Vclu{isize}.found_5pm = true;
      Vclu{isize}.thr_5pm = v_;
      Vclu{isize}.icluster_5pm = icluster;
      
      Vclu{isize}.error_5pm = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt)); 
      
    end
    if (~Vclu{isize}.found_2pm) && RMSDclu{isize}(icluster) <= 2e-3
      Vclu{isize}.found_2pm = true;
      Vclu{isize}.thr_2pm = v_;
      Vclu{isize}.icluster_2pm = icluster;
      
      
      Vclu{isize}.error_2pm = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt)); 
    end
    if (~Vclu{isize}.found_1pm) && RMSDclu{isize}(icluster) <= 1e-3
      Vclu{isize}.found_1pm = true;
      Vclu{isize}.thr_1pm = v_;
      Vclu{isize}.icluster_1pm = icluster;
      
      
      Vclu{isize}.error_1pm = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt)); 
    end
    if (~Vclu{isize}.found_500ppm) && RMSDclu{isize}(icluster) <= 5e-4
      Vclu{isize}.found_500ppm = true;
      Vclu{isize}.thr_500ppm = v_;
      Vclu{isize}.icluster_500ppm = icluster;
      
      
      Vclu{isize}.error_500ppm = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt));
    end
    if (~Vclu{isize}.found_200ppm) && RMSDclu{isize}(icluster) <= 2e-4
      Vclu{isize}.found_200ppm = true;
      Vclu{isize}.thr_200ppm = v_;
      Vclu{isize}.icluster_200ppm = icluster;
      
      
      Vclu{isize}.error_200ppm = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt));
    end
    if (~Vclu{isize}.found_100ppm) && RMSDclu{isize}(icluster) <= 1e-4
      Vclu{isize}.found_100ppm = true;
      Vclu{isize}.thr_100ppm = v_;
      Vclu{isize}.icluster_100ppm = icluster;
      
      
      Vclu{isize}.error_100ppm = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt));
    end
    if (~Vclu{isize}.found_50ppm) && RMSDclu{isize}(icluster) <= 5e-5
      Vclu{isize}.found_50ppm = true;
      Vclu{isize}.thr_50ppm = v_;
      Vclu{isize}.icluster_50ppm = icluster;
      
      
      Vclu{isize}.error_50ppm = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt));
    end
    if (~Vclu{isize}.found_20ppm) && RMSDclu{isize}(icluster) <= 2e-5
      Vclu{isize}.found_20ppm = true;
      Vclu{isize}.thr_20ppm = v_;
      Vclu{isize}.icluster_20ppm = icluster;
      
      
      Vclu{isize}.error_20ppm = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt));
    end
    if (~Vclu{isize}.found_10ppm) && RMSDclu{isize}(icluster) <= 1e-5
      Vclu{isize}.found_10ppm = true;
      Vclu{isize}.thr_10ppm = v_;
      Vclu{isize}.icluster_10ppm = icluster;
      
      
      Vclu{isize}.error_10ppm = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt));
    end
    if (~Vclu{isize}.found_5ppm) && RMSDclu{isize}(icluster) <= 5e-6
      Vclu{isize}.found_5ppm = true;
      Vclu{isize}.thr_5ppm = v_;
      Vclu{isize}.icluster_5ppm = icluster;
      
      
      Vclu{isize}.error_5ppm = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt));
    end
    if (~Vclu{isize}.found_2ppm) && RMSDclu{isize}(icluster) <= 2e-6
      Vclu{isize}.found_2ppm = true;
      Vclu{isize}.thr_2ppm = v_;
      Vclu{isize}.icluster_2ppm = icluster;
      
      
      Vclu{isize}.error_2ppm = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt));
    end
    if (~Vclu{isize}.found_1ppm) && RMSDclu{isize}(icluster) <= 1e-6
      Vclu{isize}.found_1ppm = true;
      Vclu{isize}.thr_1ppm = v_;
      Vclu{isize}.icluster_1ppm = icluster;
      
      
      Vclu{isize}.error_1ppm = (numberPossibleClusters(isize) - icluster)*sqrt(cumsum(abs(v__-v_).^2)./(1:nt));
    end
    
  end
end
ana.RMSDclu = RMSDclu;
ana.Vclu = Vclu;
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ModulationDepth = Nuclei.Statistics{1}.Modulation_Depth_p; % = matrix(N);
Hyperfine = abs(Nuclei.Statistics{1}.Hyperfine_perpendicular); % = matrix(N,1);
Nuclear_Dipole = abs(Nuclei.Statistics{1}.Nuclear_Dipole_perpendicular); % = matrix(N);
Frequency_Pair = Nuclei.Statistics{1}.Frequency_Pair_p; % = matrix(N);
DeltaHyperfine = abs(Nuclei.Statistics{1}.DeltaHyperfine_perpendicular); % = matrix(N);
bAmax = abs(Nuclei.Statistics{1}.bAmax); % = matrix(N);
Adjacency = Nuclei.Adjacency; % = matrix(N);

% ENUM
MIN = 1; MAX = 2; EDGE = 3; CRIT = 4;
ClusterH = cell(1,Method_order);
ClusterOriH = cell(nOrientations,Method_order);

isize = 1;
ClusterH{isize}.ENUM = ['SELF']; 
ClusterH{isize}.Hyperfine  = Hyperfine;


for isize = 2:Method_order
  ClusterH{isize} = getClusterHStats(Hyperfine,DeltaHyperfine, ...
    Nuclear_Dipole,bAmax, ModulationDepth,Frequency_Pair,Adjacency, ...
    Coordinates,Clusters,ClusterGeo,numberClusters,isize,TM_powder);
  
  if nOrientations==1
    for iOri = 1:nOrientations
      ModulationDepth_ori = Nuclei.Statistics{iOri}.Modulation_Depth; % = matrix(N);
      Hyperfine_ori = abs(Nuclei.Statistics{iOri}.Hyperfine); % = matrix(N,1);
      Nuclear_Dipole_ori = abs(Nuclei.Statistics{iOri}.Nuclear_Dipole); % = matrix(N);
      Frequency_Pair_ori = Nuclei.Statistics{iOri}.Frequency_Pair; % = matrix(N);
      DeltaHyperfine_ori = abs(Nuclei.Statistics{iOri}.DeltaHyperfine); % = matrix(N);
      
      
      ClusterOriH{iOri,isize} = getClusterHStats(Hyperfine_ori,DeltaHyperfine_ori,...
        Nuclear_Dipole_ori,bAmax, ModulationDepth_ori,Frequency_Pair_ori,...
        Adjacency,Coordinates,Clusters,ClusterGeo,numberClusters,isize,TM_powder);
    end
  end
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
ana.ClusterOriH = ClusterOriH;

% Save.
save( [matfile(1:end-4),'_analysis.mat'] ,'ana','-v7.3');

% Close parallel pool.
% delete(pool);
end


