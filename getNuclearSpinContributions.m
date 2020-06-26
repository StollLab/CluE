function out = getNuclearSpinContributions(matfile, ...
  Nuclei, System, nOrientations, Clusters, Signals, AuxiliarySignal,Method, experiment_time, gridWeight, TM_powder, Input, SignalMean, Order_n_SignalMean)

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
  experiment_time    = indata.experiment_time;
  gridWeight         = indata.gridWeight;
  TM_powder          = indata.TM_powder;
  Input              = indata.Input;
  SignalMean         = indata.SignalMean;
  Order_n_SignalMean = indata.Order_n_SignalMean;
end

% number of bath spin
nSpins = Nuclei.number;

% Get list electron-nucleus separations.
r = vecnorm(Nuclei.Coordinates,2,2);

% Initialize other output variables.
sansSpinV = zeros(nSpins,System.timepoints,nOrientations);
V_of_R = zeros(nSpins,System.timepoints,nOrientations);

% Get list of n, from n-CCE calculated.
orderrange = Method.order_lower_bound:Method.order;

% Determine minimum system size needed to hold each cluster.
cluster_distance = cell(1,Method.order);
for isize = orderrange
  cluster_distance{isize} = max(r(Clusters{isize}),[],2);
end

% remove current pool if it exists.
delete(gcp('nocreate'));

% determine number of cores available
numCores = feature('numcores');

% create parallel pool
pool = parpool(numCores);

% Loop thourgh orientations used in powder averaging.
parfor iOri = 1:nOrientations
  fprintf('Orientations: %d/%d.\n', iOri,nOrientations);
  
  % Loop through all bath spins.
  for ispin = 1:nSpins
    
    % Initialize to full signal.
    sansSpinV(ispin,:,iOri) = Signals{iOri};
    V_of_R(ispin,:,iOri) = gridWeight(iOri)*ones(size(experiment_time));
    
    % Loop though all cluster sizes available.
    for isize = orderrange
      
      % Calculate signal for system without spins farther away than ispin
      rclusters = cluster_distance{isize} < r(ispin);
      v_ = prod(AuxiliarySignal{iOri}{isize}(rclusters,:),1);
      V_of_R(ispin,:,iOri) = V_of_R(ispin,:,iOri).*v_;
      
      % Get list of all clusters of size isize that contain spin ispin.
      clusters = any(Clusters{isize}==ispin,2);      
      v_ = prod(AuxiliarySignal{iOri}{isize}(clusters,:),1);
      sansSpinV(ispin,:,iOri) = sansSpinV(ispin,:,iOri)./v_;
        
    end
  end
  
end

% Do powder averaging.
sansSpinVpowder = sum(sansSpinV,3);
V_of_R_powder = sum(V_of_R,3);

% Initialize TMs for systems missing spin n.
sansSpinTM = zeros(1,nSpins);
TM_of_R = zeros(1,nSpins);
for ispin = 1:nSpins
  % Determine 1/e decay times 
  sansSpinTM(ispin) = getTM(experiment_time,sansSpinVpowder(ispin,:));
  TM_of_R(ispin) = getTM(experiment_time,V_of_R_powder(ispin,:));
end

% Initialize coherence of signals evaluated at the TM of the full signal.
sansSpinVofTM = zeros(1,nSpins);

% Loop through all spins
for ispin = 1:nSpins
  % Interpolate between time grid points.
  sansSpinVofTM(ispin) = interp1(experiment_time,abs(sansSpinVpowder(ispin,:)),TM_powder);
end

% Pack output variable together. ------------------------------------------

% total evolution time
out.time = experiment_time;

% spin indices
out.N_of_R = 1:nSpins;

% number of bath spins
out.numberNuclei = nSpins;

% set of nth order CCE sims
out.Order_n_SignalMean = Order_n_SignalMean;

% electron-nucleus distances
out.R = r;

% orientation signals with spin n removed
out.sansSpinV = sansSpinV;

% orientationally averaged signals with spin n removed
out.sansSpinVpowder = sansSpinVpowder;

% 1/e times for each sansSpinVpowder
out.sansSpinTM = sansSpinTM;

% coherence of sansSpinVpowder evaluated at the TM of the full signal
out.sansSpinVofTM = sansSpinVofTM;

% original simulation input
out.SimulationInput = Input;

% original simulation TM of the full signal
out.TM = TM_powder;

% list TMs when only spin within distance R of the electron are included
out.TM_of_R = TM_of_R;

% list of atom identities
out.Type = Nuclei.Type;

% change in TM from removing spin n 
out.DeltaTM = sansSpinTM - out.TM;

% full signal
out.V = SignalMean;

% signals from keeping only spins n, where r_{n} <= R
out.V_of_R = V_of_R;

% powder signals from keeping only spins n, where r_{n} <= R
out.V_of_R_powder = V_of_R_powder;

% out file name
outfilename = [matfile(1:end-4),'TM.pdb'];
disp(outfilename)

% Save.
save( [matfile(1:end-4),'_out.mat'] ,'out');

% DeltaTM in ns 
DeltaTM_ns = out.DeltaTM*1e9; % s -> ns;

% Write DeltaTMs into the system pdb. 
writeSpinPDB(Nuclei,DeltaTM_ns,outfilename);

% Close parallel pool.
delete(pool);
end
