function out = getNuclearSpinContributions(matfile)

% Load data file.
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

nSpins = Nuclei.number;

% Get list electron-nucleus separations.
r = vecnorm(Nuclei.Coordinates,2,2);

% Initialize other output variables.
sansSpinV = zeros(nSpins,System.timepoints,nOrientations);
V_of_R = zeros(nSpins,System.timepoints,nOrientations);

orderrange = Method.order_lower_bound:Method.order;

% Determine minimum system size needed to hold each cluster.
cluster_distance = cell(1,Method.order);
for isize = orderrange
  cluster_distance{isize} = max(r(Clusters{isize}),[],2);
end

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
      v_ = prod(AuxiliarySignal{iOri}{isize}(rclusters,:));
      V_of_R(ispin,:,iOri) = V_of_R(ispin,:,iOri).*v_;
      
      % Get list of all clusters of size isize that contain spin ispin.
      clusters = any(Clusters{isize}==ispin,2);      
      v_ = prod(AuxiliarySignal{iOri}{isize}(clusters,:));
      sansSpinV(ispin,:,iOri) = sansSpinV(ispin,:,iOri)./v_;
        
    end
  end
  
end

% Do powder averaging
sansSpinVpowder = sum(sansSpinV,3);
V_of_R_powder = sum(V_of_R,3);

% Determine 1/e decay times
sansSpinTM = zeros(1,nSpins);
TM_of_R = zeros(1,nSpins);
for ispin = 1:nSpins
  sansSpinTM(ispin) = getTM(experiment_time,sansSpinVpowder(ispin,:));
  TM_of_R(ispin) = getTM(experiment_time,V_of_R_powder(ispin,:));
end

sansSpinVofTM = zeros(1,nSpins);
for ispin = 1:nSpins
  sansSpinVofTM(ispin) = interp1(experiment_time,abs(sansSpinVpowder(ispin,:)),TM_powder);
end

out.time = experiment_time;
out.N_of_R = 1:nSpins;
out.numberNuclei = nSpins;
out.Order_n_SignalMean = Order_n_SignalMean;

out.R = r;
out.sansSpinV = sansSpinV;
out.sansSpinVpowder = sansSpinVpowder;
out.sansSpinTM = sansSpinTM;
out.sansSpinVofTM = sansSpinVofTM;

out.SimulationInput = Input;

out.TM = TM_powder;
out.TM_of_R = TM_of_R;

out.DeltaTM = sansSpinTM - out.TM;

out.V = SignalMean;
out.V_of_R = V_of_R;
out.V_of_R_powder = V_of_R_powder;

outfilename = [matfile(1:end-4),'TM.pdb']
writeSpinPDB(Nuclei,out.DeltaTM,outfilename);

end
