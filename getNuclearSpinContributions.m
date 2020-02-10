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

N = Nuclei.number;

% Get sorted list electron-nucleus separations.
r = sqrt(sum(Nuclei.Coordinates.^2,2))*1e10;
[R,ind] = sort(r);

% Initialize other output variables.
sansSpinV = zeros(N,System.timepoints,nOrientations);
sansSpinVpowder = zeros(N,System.timepoints);
sansSpinTM = zeros(1,N);
sansSpinVofTM = zeros(1,N);
V_of_R = zeros(N,System.timepoints,nOrientations);
V_of_R_powder = zeros(N,System.timepoints);
TM_of_R = zeros(1,N);


time = experiment_time*1e6; % us.

Vtest = cell(1,nOrientations);
VtestPow = zeros(size(time));

% Determine minimum system size needed to hold each cluster.
cluster_distance = cell(1,Method.order);
for isize = Method.order_lower_bound:Method.order
  cluster_distance{isize} = max(r(Clusters{isize})');
end
cluster_distance{1}=r(Clusters{1})';

% Loop thourgh orientations used in powder averaging.
parfor ii = 1:nOrientations
  fprintf('Orientations: %d/%d.\n', ii,nOrientations);
  
  Vtest{ii} = ones(size(time));
  
  % Loop through all bath spins.
  for ispin = 1:N
    % Initialize to full signal.
    sansSpinV(ispin,:,ii) = Signals{ii};
    V_of_R(ispin,:,ii) = ones(size(time));
    % Loop though all cluster sizes available.
    for isize = Method.order_lower_bound:Method.order
      
      rclusters = cluster_distance{isize}<R(ispin);
      V_of_R(ispin,:,ii) =V_of_R(ispin,:,ii).*prod(AuxiliarySignal{ii}{isize}(rclusters,:));
      
      
      % Get list of all clusters of size isize that contain spin ispin.
      clusters = find(sum(Clusters{isize}==ispin,2)>0)';
      
      % Loop through clusters.
      for icluster = clusters
        
        % Get auxiliary signal for cluster icluster.
        v_ = AuxiliarySignal{ii}{isize}(icluster,:);
        Vtest{ii} = Vtest{ii}.*v_.^(1/isize);
        sansSpinV(ispin,:,ii) = sansSpinV(ispin,:,ii)./v_;
        
        
      end
    end
  end
  
end

for ii = 1:nOrientations
  VtestPow = VtestPow + Vtest{ii}*gridWeight(ii);
end

for ispin = 1:N
  for ii = 1:nOrientations
    sansSpinVpowder(ispin,:) = sansSpinVpowder(ispin,:) + sansSpinV(ispin,:,ii);
    V_of_R_powder(ispin,:) = V_of_R_powder(ispin,:) + V_of_R(ispin,:,ii);
  end
  sansSpinTM(ispin) = getTM(time,sansSpinVpowder(ispin,:));
  TM_of_R(ispin) = getTM(time,V_of_R_powder(ispin,:));
  
  sansSpinVofTM(ispin) = linearEval(experiment_time,abs( sansSpinVpowder(ispin,:) ),TM_powder);
end

out.index = ind;
out.time = time;

out.N_of_R = 1:N;
out.numberNuclei = N;

out.Order_n_SignalMean = Order_n_SignalMean;

out.R = R;
out.sansSpinV = sansSpinV(ind,:,:);
out.sansSpinVpowder = sansSpinVpowder(ind,:);
out.sansSpinTM = sansSpinTM(ind);
out.sansSpinVofTM = sansSpinVofTM(ind);

out.SimulationInput = Input;

out.TM = TM_powder*1e6;
out.TM_of_R = TM_of_R;

out.DeltaTM = sansSpinTM - out.TM;


out.V = SignalMean;
out.V_of_R = V_of_R(ind,:,:);
out.V_of_R_powder = V_of_R_powder(ind,:);
outfilename = [matfile(1:end-4),'TM.pdb']
writeSpinPDB(Nuclei,out.DeltaTM,outfilename)
end

function y0 = linearEval(x,y,x0)
N = length(x);

n = sum(x<x0);

if n == 0

  dx = x0 - x(1);
  Dx = x(2) - x(1);
  dy = y(2) - y(1);
  y0 = y(1) + dy/Dx*dx;

elseif n ==N
  
  dx = x0 - x(n);
  Dx = x(N) - x(N-1);
  dy = y(N) - y(N-1);
  y0 = y(N) + dy/Dx*dx;
  
else
  
  dx = x0 - x(n);
  Dx = x(n+1) - x(n);

  p = 1 - dx/Dx;
  
  y0 = p*y(n) + (1-p)*y(n+1);
end

end

