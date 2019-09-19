function [Parameters, Eta] = getConvergenceParameters(pdb_file,twotau_max, Electron_Coor,options)

% initializing ==========================================================
Data.InputData =  pdb_file;

if isfield(options,'System')
  System = options.System;
end
if isfield(options,'Method')
  Method = options.Method;
end

timepoints = 101;
System.dt = twotau_max/(timepoints-1)/2; % s.
System.timepoints = timepoints; % s.

System.Electron.Coordinates = Electron_Coor;
angstrom = 1e-10; % m.
System.gridSize = 1;
if isfield(options,'path')
  oldpath = path;
  path(options.path,oldpath);
end

% =======================================================================
% r0, radius, method, order, orientation, D2O
% =======================================================================
Method.r0 = options.r0; % m.
System.radius = options.radius; % m.
Method.order = 2;

% Set lab magnetic field.
if isfield(options,'magneticField')
  System.magneticField = options.magneticField;
end

% Set number of orientations to average over.
if isfield(options,'n_orientations')
  n_orientations = options.n_orientations;
else
  n_orientations = 6;
end

% Toggle verbose mode.
if isfield(options,'verbose')
  Method.verbose = options.verbose;
end

if isfield(options,'conserveMemory')
  Method.conserveMemory = options.conserveMemory;
end
% Orientation size.
System.averaging = 'none';

% Turn all water molecules into heavy water.
if isfield(options,'D2O')
  System.D2O = options.D2O;
else
  System.D2O = false;
end

% Turn hydrogen into deuterium.
if isfield(options,'deuterateAll')
  System.deuterateAll = options.deuterateAll;
else
  System.deuterateAll = false;
end

% Ignore everything save for water.
if isfield(options,'solventOnly')
  System.solventOnly = options.solventOnly;
else
  System.solventOnly = false;
end

% Set the fraction of protium in the sample in the protein an/or the solvent.
if isfield(options,'protiumFractionProtein')
  System.protiumFractionProtein = options.protiumFractionProtein;
end
if isfield(options,'protiumFractionSolvent')
  System.protiumFractionSolvent = options.protiumFractionSolvent;
end

% Select the clustering method.
if isfield(options,'method')
  Method.method = options.method;
elseif (~System.deuterateAll && ~System.D2O)
  Method.method = 'rCCE';
else
  Method.method = 'CCE';
end

% Select the highest cluster size.
if isfield(options,'maxOrder')
  maxOrder = options.maxOrder;
else
  maxOrder = 2;
end

% Look at the n+1 cluster size for error,
% but do not calculate more than one orientation.
if isfield(options,'orderPlus1')
  checkOrderPlus1 = options.orderPlus1;
else
  checkOrderPlus1 = false;
end
Method.parallelComputing = false;
Method.precalculateHamiltonian = false;

% Set convergence threshold.
if isfield(options,'threshold')
  threshold = options.threshold;
else
  threshold = 1e-3;
end

% Use a different threshold on the time consuming steps.
if isfield(options,'order_threshold')
  order_threshold = options.order_threshold;
else
  order_threshold = 1e-2;
end
if isfield(options,'orientation_threshold')
  orientation_threshold = options.orientation_threshold;
else
  orientation_threshold = 1e-2;
end
if isfield(options,'maxOrientations')
  maxOrientations = options.maxOrientations;
else
  maxOrientations = 5810;
end
if isfield(options,'max_R')
  max_R = options.max_R;
else
  max_R = inf;
end
if isfield(options,'max_r0')
  max_r0 = options.max_r0;
else
  max_r0 = inf;
end

% system search =========================================================
search_system = true;


savefile_name = ['SIM_' Data.InputData];
savefile_name = savefile_name(1:end-4);
savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) 'A_averaging_' System.averaging '_' num2str(System.gridSize) '.mat'];
Data.OutputData = savefile_name;
disp(savefile_name);
computeSignal = true;

% check for pre-existing file
if exist(savefile_name,'file') == 2
  load(savefile_name,'SignalMean','TM_powder','Progress')
  
  % do not repeat the simulation needlessly
  if Progress.complete
    computeSignal = false;
  end
end

% calculate signal
if computeSignal
  [SignalMean, ~, ~,~, TM_powder, Progress] = nuclear_spin_diffusion(System,Method,Data);
end

% initialize error variable
eta_sys = [];

while search_system
  
  % save the initial results for comparison
  SignalMean_sys = SignalMean;
  TM_sys = TM_powder;
  
  % radius search =========================================================
  search_R = true;
  R_ = System.radius;
  eta_R = [];
  
  % search system size
  while search_R
    
    % save the the previous itteration for comparison
    SignalMean_ = SignalMean; TM_ = TM_powder;
    R_ = System.radius;
    
    % increase system size
    System.radius = System.radius + 1*angstrom;
    
    if System.radius > max_R
      fprintf('Caution radius not converged by the maximum allowed value of %d m.\n', R_);
      break;
    end
    
    % save
    savefile_name = ['SIM_' Data.InputData];
    savefile_name = savefile_name(1:end-4);
    savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) 'A_averaging_' System.averaging '_' num2str(System.gridSize) '.mat'];
    Data.OutputData = savefile_name;
    disp(savefile_name);
    
    computeSignal = true;
    
    % check for pre-existing file
    if exist(savefile_name,'file') == 2
      load(savefile_name,'SignalMean','TM_powder','Progress')
      
      % do not repeat the simulation needlessly
      if Progress.complete
        computeSignal = false;
      end
      
    end
    
    if computeSignal
      [SignalMean, ~, ~,~, TM_powder, Progress] = nuclear_spin_diffusion(System,Method,Data);
    end
    
    % calculate fractional difference
    eta = 2*abs(TM_powder-TM_)/(TM_powder+TM_); %max( abs( SignalMean - SignalMean_ ) );
    
    % display step results
    fprintf('eta(R = %d A) = %f.\n',R_*1e10,eta);
    fprintf('TM(R = %d A) = %f us.\n',R_*1e10, TM_*1e6);
    fprintf('TM(R = %d A) = %f us.\n',System.radius*1e10, TM_powder*1e6);
    
    % save error erm 
    eta_R = [eta_R,eta];
    
    % check for convergence
    if eta < threshold
      search_R = false;
    end
    
  end
  
  % set system size to converged value
  System.radius = R_;
  
  % r0 search =============================================================
  search_r0 = true;
  r0_ = System.radius;
  eta_r0 = [];
  
  savefile_name = ['SIM_' Data.InputData];
  savefile_name = savefile_name(1:end-4);
  savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) 'A_averaging_' System.averaging '_' num2str(System.gridSize) '.mat'];
  Data.OutputData = savefile_name;
  
  disp(savefile_name);
  
  computeSignal = true;
  
  % check for pre-existing file
  if exist(savefile_name,'file') == 2
    load(savefile_name,'SignalMean','TM_powder','Progress')
 
    % do not needlessly repeat computations
    if Progress.complete
      computeSignal = false;
    end
    
  end
  
  % simulate coherence decay
  if computeSignal
    [SignalMean, ~, ~,~, TM_powder, Progress] = nuclear_spin_diffusion(System,Method,Data);
  end
  
  % find a converged r0
  while search_r0
    % save previous loop values
    SignalMean_ = SignalMean; TM_ = TM_powder;
    r0_ = Method.r0;
    
    % increase r0
    Method.r0 = Method.r0 + 1*angstrom;
    
    if Method.r0 > max_r0
      fprintf('Caution r0 not converged by maximum allowed value of %d m.\n', R_);
      break;
    end
    
    % save
    savefile_name = ['SIM_' Data.InputData];
    savefile_name = savefile_name(1:end-4);
    savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) 'A_averaging_' System.averaging '_' num2str(System.gridSize) '.mat'];
    Data.OutputData = savefile_name;
    
    disp(savefile_name);
    
    computeSignal = true;
    
    % check for existing file
    if exist(savefile_name,'file') == 2
      load(savefile_name,'SignalMean','TM_powder','Progress')
      
      % do not needlessly repeat computations
      if Progress.complete
        computeSignal = false;
      end
      
    end
    
    % compute decay curve
    if computeSignal
      [SignalMean, ~, ~,~, TM_powder, Progress] = nuclear_spin_diffusion(System,Method,Data);
    end
    
    % fractional difference
    eta = 2*abs(TM_powder-TM_)/(TM_powder+TM_); %max( abs( SignalMean - SignalMean_ ) );
    
    % display loop info
    fprintf('eta(r0 = %d A) = %f.\n',r0_*1e10, eta);
    fprintf('TM(r0 = %d A) = %f us.\n',r0_*1e10, TM_*1e6);
    fprintf('TM(r0 = %d A) = %f us.\n',Method.r0*1e10, TM_powder*1e6);
    
    % save error term
    eta_r0 = [eta_r0,eta];
    
    % check for convergence
    if eta < threshold
      search_r0 = false;
    end
    
  end
  
  % set r0 to converged value
  Method.r0 = r0_;
  
  % order search ==========================================================
  
  search_order = true;
  didReachMaxOrder = false;
  
  % switch to the appropriate higher order method
  if strcmp(Method.method, 'rCCE')
    Method.method = 'CCE';
  elseif strcmp(Method.method, 'rCE')
    Method.method = 'CE';
  end
  
  
  order_ = Method.order;
  eta_order = [];
  savefile_name = ['SIM_' Data.InputData];
  savefile_name = savefile_name(1:end-4);
  savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) 'A_averaging_' System.averaging '_' num2str(System.gridSize) '.mat'];
  Data.OutputData = savefile_name;
  
  disp(savefile_name);
  
  computeSignal = true;
  
  % check for pre-computed file
  if exist(savefile_name,'file') == 2
    load(savefile_name,'SignalMean','TM_powder','Progress')
  
    % do not needlessly repeat computations
    if Progress.complete
      computeSignal = false;
    end
    
  end
  
  % simulate coherence decay
  if computeSignal
    [SignalMean, ~, ~,~, TM_powder, Progress] = nuclear_spin_diffusion(System,Method,Data);
  end
  
  % odd orders are sometimes not stable
  isSignalNonsense = false;
  
  % find minimum needed order
  while search_order
    
    SignalMean_ = SignalMean; TM_ = TM_powder;
    
    % odd orders are sometimes not stable
    if ~isSignalNonsense
      order_ = Method.order;
    else
      isSignalNonsense = false;
    end
    
    Method.order = Method.order + 1;
    
    if Method.order>maxOrder
      disp('Caution order not converged by maximum order.');
      didReachMaxOrder = true;
      break;
    end
    
    savefile_name = ['SIM_' Data.InputData];
    savefile_name = savefile_name(1:end-4);
    savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) 'A_averaging_' System.averaging '_' num2str(System.gridSize) '.mat'];
    Data.OutputData = savefile_name;
    
    disp(savefile_name);
    
    computeSignal = true;
    
    % check for pre-computed file
    if exist(savefile_name,'file') == 2
      load(savefile_name,'SignalMean','TM_powder','Progress')
      
      % do not needlessly repeat computations
      if Progress.complete
        computeSignal = false;
      end
      
    end
    
    % simulate decay 
    if computeSignal
      [SignalMean, ~, ~,~, TM_powder, Progress] = nuclear_spin_diffusion(System,Method,Data);
    end
    
    % fractional difference in TM
    eta = 2*abs(TM_powder-TM_)/(TM_powder+TM_); %max( abs( SignalMean - SignalMean_ ) );
    
    % display loop info
    fprintf('eta(|C| = %d) = %f.\n',order_, eta);
    fprintf('TM(|C| = %d) = %f us.\n',order_, TM_*1e6);
    fprintf('TM(|C| = %d) = %f us.\n',Method.order, TM_powder*1e6);
    
    % save errro
    eta_order = [eta_order,eta];
 
    % check for nonsense simulation: v(t) <= 1, for all t 
    if abs(SignalMean(1) - max(SignalMean))>1e-9
      fprintf('Signal of order %d is nonsensical.  Referencing convergence to order %d.\n',Method.order, order_);
      SignalMean = SignalMean_;
      isSignalNonsense = true;
    end
    
    % check for convergence
    if eta < order_threshold
      search_order = false;
    end
    
  end
  Method.order = order_;
  
  % order plus 1 ========================================================
  % calculate the error term for the next order, but do not change
  % convergence parameters 
  
  if checkOrderPlus1 && didReachMaxOrder
    Method.order = Method.order + 1;
    savefile_name = ['SIM_' Data.InputData];
    savefile_name = savefile_name(1:end-4);
    savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) 'A_averaging_' System.averaging '_' num2str(System.gridSize) '.mat'];
    Data.OutputData = savefile_name;
    disp(savefile_name);
    computeSignal = true;
    if exist(savefile_name,'file') == 2
      load(savefile_name,'SignalMean','TM_powder','Progress')
      if Progress.complete
        computeSignal = false;
      end
    end
    if computeSignal
      [SignalMean, ~, ~,~, TM_powder, Progress] = nuclear_spin_diffusion(System,Method,Data);
    end
    
    eta = 2*abs(TM_powder-TM_)/(TM_powder+TM_); %max( abs( SignalMean - SignalMean_ ) );
    fprintf('eta(|C| = %d) = %f.\n',order_, eta);
    fprintf('TM(|C| = %d) = %f us.\n',order_, TM_*1e6);
    fprintf('TM(|C| = %d) = %f us.\n',Method.order, TM_powder*1e6);
    eta_order = [eta_order,eta];
  end
  Method.order = order_;
  % system check=========================================================
 
  savefile_name = ['SIM_' Data.InputData];
  savefile_name = savefile_name(1:end-4);
  savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) 'A_averaging_' System.averaging '_' num2str(System.gridSize) '.mat'];
  Data.OutputData = savefile_name;
  disp(savefile_name);
  computeSignal = true;
  if exist(savefile_name,'file') == 2
    load(savefile_name,'SignalMean','TM_powder','Progress')
    if Progress.complete
      computeSignal = false;
    end
  end
  if computeSignal
    [SignalMean, ~, ~,~, TM_powder, Progress] = nuclear_spin_diffusion(System,Method,Data);
  end
  
  eta = 2*abs(TM_powder-TM_sys)/(TM_powder+TM_sys);
  fprintf('eta(sys) = %f.\n',eta);
  fprintf('TM(sys0) = %f us.\n',TM_sys*1e6);
  fprintf('TM(sys) = %f us.\n',TM_powder*1e6);
  eta_sys = [eta_sys,eta];
  if eta < threshold
    search_system = false;
  end
  
  
  % method check ==========================================================
  % check if rCE/rCCE may be used instead of CE/CCE

  if order_==2 && (~System.deuterateAll && ~System.D2O)
    method_ = Method.method;
    
    % set save file name
    savefile_name = ['SIM_' Data.InputData];
    savefile_name = savefile_name(1:end-4);
    savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) 'A_averaging_' System.averaging '_' num2str(System.gridSize) '.mat'];
    Data.OutputData = savefile_name;
    disp(savefile_name);
    
    computeSignal = true;
    
    
    % check for pre-computed file
    if exist(savefile_name,'file') == 2
      load(savefile_name,'SignalMean','TM_powder','Progress')
      
      % do not needlessly repeat computations
      if Progress.complete
        computeSignal = false;
      end
      
    end
    
    % simulate coherence decay
    if computeSignal
      [SignalMean, ~, ~,~, TM_powder, Progress] = nuclear_spin_diffusion(System,Method,Data);
    end
    
    SignalMean_ = SignalMean; TM_ = TM_powder;
    
    if strcmp(Method.method, 'CCE')
      Method.method = 'rCCE';
    elseif (~System.deuterateAll && ~System.D2O) && strcmp(Method.method, 'CE')
      Method.method = 'rCCE';
    end
    
    savefile_name = ['SIM_' Data.InputData];
    savefile_name = savefile_name(1:end-4);
    savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) 'A_averaging_' System.averaging '_' num2str(System.gridSize) '.mat'];
    Data.OutputData = savefile_name;
    
    disp(savefile_name);
    
    computeSignal = true;
    
    % check for pre-computed file
    if exist(savefile_name,'file') == 2
      load(savefile_name,'SignalMean','TM_powder','Progress');
      
      % do not needlessly repeat computations
      if Progress.complete
        computeSignal = false; 
      end
      
    end
    
    % compute signal
    if computeSignal
      [SignalMean, ~, ~,~, TM_powder, Progress] = nuclear_spin_diffusion(System,Method,Data);
    end
    
    
    % fractional difference
    eta = 2*abs(TM_powder-TM_)/(TM_powder+TM_); %max( abs( SignalMean - SignalMean_ ) );
    
    % display loop info
    fprintf('eta(rCCE) = %f.\n',eta);
    fprintf('TM(rCCE) = %f us.\n', TM_*1e6);
    fprintf('TM(CCE) = %f us.\n', TM_powder*1e6);
    
    % check for convergence
    if eta >= threshold
      Method.method = method_;
    end
    
  end
end

% orienttion search =====================================================
% powder average
if maxOrientations==1
  fprintf('Using only the given orientation.\n');
  eta_orientation = nan;
else
  
  % use parallel computing
  Method.parallelComputing = true;
  
  % lebedev grid options
  gridopt= [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074,...
    3470, 3890, 4334, 4802, 5294, 5810];
  
  % grid option index
  igrid = find(gridopt == n_orientations);
  
  % start at the begining if no valid index is found
  if isempty(igrid)
    igrid = 1;
  end
  
  System.averaging = 'powder';
  search_orientation = true;
  eta_orientation = [];
  
  System.gridSize = gridopt(igrid);
  
  % set save file name
  savefile_name = ['SIM_' Data.InputData];
  savefile_name = savefile_name(1:end-4);
  savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) 'A_averaging_' System.averaging '_' num2str(System.gridSize) '.mat'];
  Data.OutputData = savefile_name;
  
  disp(savefile_name);
  
  computeSignal = true;
  
  % check for pre-computed file
  if exist(savefile_name,'file') == 2
    load(savefile_name,'SignalMean','TM_powder','Progress');
    
     % do not needlessly repeat computations
    if Progress.complete
      computeSignal = false;
    end
  end
  
  if computeSignal
    [SignalMean, ~, ~,~, TM_powder, Progress] = nuclear_spin_diffusion(System,Method,Data);
  end
  
  didReachMaxOrientation = false;
  
  % loop over orientation
  while search_orientation
    
    SignalMean_ = SignalMean; TM_ = TM_powder;
    
    orientations_ = igrid;
    
    igrid = igrid+1;
    
    if gridopt(igrid) > maxOrientations
      disp('Caution orientations not converged by maximum allowed orientation.');
      didReachMaxOrientation = true;
      break;
    end
    
    System.gridSize = gridopt(igrid);
    
    savefile_name = ['SIM_' Data.InputData];
    savefile_name = savefile_name(1:end-4);
    savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) 'A_averaging_' System.averaging '_' num2str(System.gridSize) '.mat'];
    Data.OutputData = savefile_name;
    
    disp(savefile_name);
    
    computeSignal = true;
    
    % check for pre-computed file
    if exist(savefile_name,'file') == 2
      load(savefile_name,'SignalMean','TM_powder','Progress');
      
      if Progress.complete
        computeSignal = false;
      end
      
    end
    
    if computeSignal
      [SignalMean, ~, ~,~, TM_powder, Progress] = nuclear_spin_diffusion(System,Method,Data);
    end
    
    eta = 2*abs(TM_powder-TM_)/(TM_powder+TM_); %max( abs( SignalMean - SignalMean_ ) );
    
    ngrid = 1;
    
    if orientations_>0
      ngrid = gridopt(orientations_);
    end
    
    fprintf('eta(n = %d) = %f.\n',ngrid,eta);
    fprintf('TM(n = %d) = %f us.\n',ngrid, TM_*1e6);
    fprintf('TM(n = %d) = %f us.\n',System.gridSize, TM_powder*1e6);
    
    eta_orientation = [eta_orientation,eta];
    
    if eta < orientation_threshold
      search_orientation = false;
    end
  end
  
  if orientations_>0
    System.gridSize = gridopt(orientations_);
  else
    System.gridSize = 1;
  end
  
end

% collect all error vectors together
Eta = {eta_R, eta_r0, eta_order, eta_sys, eta_orientation};

% save convergence parameters
Parameters.R = R_;
Parameters.r0 = r0_;
Parameters.order = order_;
Parameters.orientations = System.gridSize;


% calculate the final signal with higher resolution
if isfield(options,'timepoints_final')
  timepoints = options.timepoints_final;
else 
  timepoints = 1001;
end

System.dt = twotau_max/(timepoints-1)/2; % s.
System.timepoints = timepoints; % s.

timestr = ['nt_', num2str(timepoints), 'dt_', num2str(System.dt*1e9), '_ns'];

% set save file name
savefile_name = ['SIM_' Data.InputData];
savefile_name = savefile_name(1:end-4);
savefile_name = [savefile_name '_' num2str(Method.order) Method.method '_R_' ...
                 num2str(System.radius*1e10) 'A_r0_' num2str(Method.r0*1e10) ...
                 'A_av_' System.averaging '_' num2str(System.gridSize), ...
                 timestr '.mat'];
               
Data.OutputData = savefile_name;

% calculate blame factors
Method.getNuclearContributions = true;

computeSignal = true;

if exist(savefile_name,'file') == 2
  load(savefile_name,'SignalMean','TM_powder','Progress');
  
  if Progress.complete
    computeSignal = false;
  end
  
end

if computeSignal
  [~, ~, ~,~, ~, ~] = nuclear_spin_diffusion(System,Method,Data);
end

fprintf('=============================================================================\n');
fprintf('----------------------------Convergence Prameters----------------------------\n');
fprintf('R = %d m.\n',R_);
fprintf('r0 = %d m.\n',r0_);
fprintf('order = %d.\n',order_);
fprintf('orientations = %d.\n',System.gridSize);
fprintf('----------------------------------Save File----------------------------------\n');
disp(savefile_name);
fprintf('------------------------------------ TM -------------------------------------\n');
fprintf('TM = %d us.\n',TM_*1e6);
fprintf('=============================================================================\n');
save(['convergenceResults_' pdb_file(1:end-4) '.mat'],'Parameters','Eta');

end
