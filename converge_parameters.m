function [parameters,System,Method,Data]  = converge_parameters(System,Method,Data, options)
fileID = fopen('convergence_log.txt','w');

% ENUM
RADIUS  =1;  %DIPOLE = 2; BAMAX = 3; POWDER = 4;

if ~isfield(Method.cutoff,'dipole')
  Method.cutoff.dipole = -inf;
end
if ~isfield(Method.cutoff,'bAmax')
  Method.cutoff.bAmax = -inf;
end

options = setDefaults(options);
metric = options.metric;

if options.pruneClusters
  Method.getClusterContributions = true;
  Method.getClusterContributions = true;
end

if options.usePseudoGrad
  Method.getUncertainty = true;
end

System.gridSize = 1;
ID = 1;
ID_Ref = ID;
hline = '-----------------------------------------------------------------\n';
 
logPrint(fileID,hline);
logPrint(fileID,'Starting.\n');

Data0 = Data;
Data.OutputData = [Data0.OutputData,'_ID_0'];
calculate_signal = true;

if isfile([Data.OutputData,'.mat'])
  try
    load([Data.OutputData,'.mat'],'SignalMean','experiment_time','TM_powder','Progress','uncertainty');
    if Progress.complete
      calculate_signal = false;
      TM_powder_ = TM_powder;
    end
  catch
  end
end
if calculate_signal
  [SignalMean, ~, TM_powder_,~,~,uncertainty] = CluE(System,Method,Data);
end
Npreload = 16;
v = zeros(Npreload,size(SignalMean,2)); 

logPrint(fileID,'TM2 = %d s.\n',TM_powder_);
logPrint(fileID,hline);
 

% Initialize error vector.
Eta = inf(1,4);

% Set error on unused cutoffs to -inf, so they will not be looked at.
Eta(~options.useCutoffs) = -inf;

% Set powder error to zero, so that it will be ignored until the end.
usePowder = options.converge.powder;
if usePowder
  Eta(end) = -inf;
end

n_ = length(options.useCutoffs);
isThisCutoffConverged = false(1,n_);

useCutoffs = options.useCutoffs;
isThisCutoffConverged(~useCutoffs) = true;

is_singleOri_converged = false;
is_powder_converged = false;
is_converged = false;
cutoff.delta = options.Delta; 
cutoff.Max = options.Max; 
cutoff.Min = options.Min; 
nextCutoff = options.firstCutoff;

use_radiusfirst = true;
radiusfirst = use_radiusfirst;

[ParameterLog, ~ , ~, NameLog] = updateParameterLog([], [],System,Method, Data0.OutputData ,cutoff,ID);

% Run until all cutoffs are converged.
while ~is_converged
  
  % Perturb cutoff. 
  [System,Method, cutoff]= adjustCutoff('relax',System, Method,cutoff,nextCutoff,uncertainty);
  
  % Print info. 
  logPrint(fileID,[cutoff.name, ' = %d ',cutoff.units,'.\n'],cutoff.value);
  logPrint(fileID, 'r = %d A. b = %d Hz. bAmax = %d. nOri = %d. \n',System.radius,...
    Method.cutoff.dipole(1),Method.cutoff.bAmax(1), System.gridSize);
  
  % Check for out of bounds parameters.
  if false %cutoff.value < cutoff.Min(cutoff.ID) || cutoff.value > cutoff.Max(cutoff.ID)
    [System, Method, cutoff]= adjustCutoff('tighten',System, Method,cutoff,nextCutoff,uncertainty);
    logPrint(fileID,[cutoff.name, ' did not converge within set bounds.\n']);
    Eta(cutoff.ID) = -inf;
    if useCutoffs(cutoff.ID)
      continue;
    end
    error('An unused parameter is being tested.');
  end
  
  % Remember parameter sets for data rerieval.
  [ParameterLog, ID , OutputData, NameLog] = updateParameterLog(ParameterLog, NameLog ,System,Method, Data0.OutputData ,cutoff,ID);
  Data.OutputData = OutputData;
  
  % Either load or calculate the appropriate simulation. 
  calculate_signal = true;
  uncertainty_ = uncertainty;
  if isfile([Data.OutputData,'.mat'])
    try
      load([Data.OutputData,'.mat'],'SignalMean','experiment_time','TM_powder','Progress','uncertainty');
      if Progress.complete
        calculate_signal = false;
      end
    catch
    end
  end
  if calculate_signal
    [SignalMean, experiment_time, TM_powder,~,~,uncertainty] = CluE(System,Method,Data);
  end
  
  % Store data.
  v(ID,:) = SignalMean;
  
  % Get error metric.
  if strcmp(nextCutoff,'pseudograd')
    eta = uncertainty.err_norm(Method.order);
  else
    eta = getErrorMetric(v(ID_Ref,:),v(ID,:),metric,experiment_time,experiment_time,options);
  end
  Eta(cutoff.ID) = eta;
  
  % Print info.
  logPrint(fileID,'TM2 = %d s.\n',TM_powder);
  logPrint(fileID,'eta = %d.\n',eta);
  logPrint(fileID,hline);
   
  % Check if data should be plotted.
  if options.doPlot
    plot(experiment_time*1e6,abs(SignalMean));
    xlabel('time of echo (\mus)');
    ylabel('coherence');
    drawnow;
  end
  
  % Check if the single orientation or powder avergage is being converged.
  if is_singleOri_converged
    
    % Check if the powder average converged.
    is_powder_converged = strcmp(nextCutoff,'powder') && eta < options.Threshold(cutoff.ID);
    
    % Update ID.
    ID_Ref = ~is_powder_converged*ID + is_powder_converged*ID_Ref;
    
    if is_powder_converged
      [System, Method, cutoff]= adjustCutoff('tighten',System, Method,cutoff,nextCutoff,uncertainty_);
    end
  else
    
    % Check if the current cutoff parameter is converged in isolation.
    isThisCutoffConverged(cutoff.ID) = eta < options.Threshold(cutoff.ID);
    
    
    if isThisCutoffConverged(cutoff.ID)
      % Revert to the previous state.
      [System, Method, cutoff]= adjustCutoff('tighten',System, Method,cutoff,nextCutoff,uncertainty_);
      if cutoff.ID == RADIUS
        radiusfirst = false;
      end
      
      % Decide if the other cutoffs need to be rechecked.
      isCutoffUnchanged = strcmp(nextCutoff_,nextCutoff);
      if isCutoffUnchanged
        isThisCutoffConverged(useCutoffs) = false;  
        isThisCutoffConverged(cutoff.ID) = true;
      end
      
    else     
      
      % Update ID.
      ID_Ref = ID;
      
      % Update cutoff convegence statuses. 
      isThisCutoffConverged(:) = false;
      if cutoff.ID == RADIUS
        radiusfirst = use_radiusfirst;
      end        
    end
  
    % Check if the single orientation system is converged.
    is_singleOri_converged =  all(Eta(useCutoffs) < options.Threshold(useCutoffs))...
      &&  all( isThisCutoffConverged( useCutoffs(1:end-1) ) );
  end
  

  % Record this cutoff.
  nextCutoff_ = nextCutoff;
  
  % Only use this black when converging a single orientation.
  if ~strcmp(nextCutoff,'powder')
    
    % Check if the orientations should be converged.
    if is_singleOri_converged
      nextCutoff = 'powder';
      Method.parallelComputing = options.parpow;
      is_powder_converged = ~usePowder;
    
    % Do not switch parameters until system size is re-converged.  
    elseif ~radiusfirst
      
      % Find the least converged cutoff.
      newMaxEta_ = max(Eta(useCutoffs & ~isThisCutoffConverged));
      
      if ~isempty(newMaxEta_)
        % Find the cutoff that is likely to produce the largest change.
        nextCutoff_ID = find(Eta == newMaxEta_);
        
        nextCutoff = options.cutoffNames{nextCutoff_ID(1)};
      end
    end
    
    % Error check.
    if isempty(Eta(useCutoffs & ~isThisCutoffConverged)) && ~is_singleOri_converged
      error('Single orientation is both converged and not converged.');
    end
    
    % Error check.
    if strcmp(nextCutoff_,nextCutoff) && eta < options.Threshold(cutoff.ID)
      error('Cutoff cannot increment.')
    end
  end
  
  % Check for convergence.
  is_converged  = is_singleOri_converged  && is_powder_converged;
  
  
end
  


parameters.radius = System.radius;
if options.converge.dipole || options.usePseudoGrad
parameters.dipole = Method.cutoff.dipole;
end
if options.converge.bAmax || options.usePseudoGrad
parameters.bAmax = Method.cutoff.bAmax;
end

parameters.gridSize = System.gridSize;

fclose(fileID);
end

function options = setDefaults(options)

if ~isfield(options,'metric')
  options.metric = 'rms';
end
if ~isfield(options,'pruneClusters')
  options.pruneClusters = false;
end
if ~isfield(options,'vmin')
  options.vmin = exp(-5);
end

options.cutoffNames = {'radius','dipole','bAmax','powder'};
% ENUM
RADIUS  =1;  DIPOLE = 2; BAMAX = 3; POWDER = 4; 

options.useCutoffs = false(1,POWDER-1);

% system radius
if ~isfield(options.converge,'radius')
  options.converge.radius = true;
end
options.useCutoffs(RADIUS) = options.converge.radius;

if ~isfield(options.threshold,'radius')
  options.threshold.radius = 1e-3;
end
if ~isfield(options.delta,'radius')
  options.delta.radius = 1e-10; % m
end
if ~isfield(options.limit,'radius')
  options.limit.radius = 20e-10; % m
end

% neighbor cutoff

% dipole
if ~isfield(options.converge,'dipole')
  options.converge.dipole = false;
end
options.useCutoffs(DIPOLE) = options.converge.dipole;

if ~isfield(options.threshold,'dipole')
  options.threshold.dipole = 1e-3;
end
if ~isfield(options.delta,'dipole')
  options.delta.dipole = -2; 
end
if ~isfield(options.limit,'dipole')
  options.limit.dipole = -15; 
end

if ~isfield(options,'usePseudoGrad')
  options.usePseudoGrad = false;
end
if ~isfield(options.threshold,'pgrad')
  options.threshold.pgrad = 1e-3;
end
if ~isfield(options.delta,'pgrad')
  options.delta.pgrad = -1000;
end
if ~isfield(options.limit,'pgrad')
  options.limit.pgrad = -15; 
end

if options.usePseudoGrad
  options.cutoffNames = {'radius','pseudograd','NA','powder'};
  options.converge.bAmax = false;
  options.threshold.dipole = options.threshold.pgrad;
  options.delta.dipole = options.delta.pgrad;
  options.limit.dipole = options.limit.pgrad;
end
% bAmax
if ~isfield(options.converge,'bAmax')
  options.converge.bAmax = false;
end
options.useCutoffs(BAMAX) = options.converge.bAmax;
if ~isfield(options.threshold,'bAmax')
  options.threshold.bAmax = 1e-3;
end
if ~isfield(options.delta,'bAmax')
  options.delta.bAmax = -4; 
end
if ~isfield(options.limit,'bAmax')
  options.limit.bAmax = -15; 
end


% powder
if ~isfield(options.converge,'powder')
  options.converge.powder = true;
end
options.useCutoffs(POWDER) = options.converge.powder;

if ~isfield(options.threshold,'powder')
  options.threshold.powder = 1e-2;
end
if ~isfield(options.limit,'grid_points')
  options.limit.grid_points = 5810;
end
if ~isfield(options,'parpow')
  options.parpow = true;
end
options.Threshold = [options.threshold.radius, options.threshold.dipole, options.threshold.bAmax,options.threshold.powder];
options.Delta = [options.delta.radius, options.delta.dipole, options.delta.bAmax,1];
options.Min = [0, 10^options.limit.dipole, 10^options.limit.bAmax,1];
options.Max = [options.limit.radius, inf, inf,options.limit.grid_points];
if options.converge.radius 
  options.firstCutoff = 'radius';
elseif options.converge.dipole
  options.firstCutoff = 'dipole';
elseif options.converge.bAmax
  options.firstCutoff = 'bAmax';
else
  options.firstCutoff = 'powder';
end
 

if ~isfield(options,'doPlot')
  options.doPlot = false;
end
end

function [ParameterLog_out, ID_out , OutputData, NameLog_out] = updateParameterLog(ParameterLog, NameLog ,System,Method, OutputData0 ,cutoff,ID)

% Initialize log variables.
ParameterLog_out = ParameterLog;
NameLog_out = NameLog;

% Update parameter log.
if isempty(ParameterLog)
  ParameterLog_out = [System.radius, Method.cutoff.dipole, Method.cutoff.bAmax, System.gridSize];
  OutputData = OutputData0;
  NameLog_out{1} = OutputData;
  ID_out = ID;
  return;
end

ParameterLog_out(ID+1, :) = [System.radius, Method.cutoff.dipole, Method.cutoff.bAmax, System.gridSize];

% Look for previous landings on these parameters.
doParametersExist = all(ParameterLog_out == ParameterLog_out(ID+1, :),2);


ID_out = ID + 1;

% Output previous data or update the parameters
if sum(doParametersExist) > 1
  ID_ref = find(doParametersExist);
  ID_ref = ID_ref(1);
  ParameterLog_out = ParameterLog;
  OutputData = NameLog{ID_ref};
else
  OutputData = [OutputData0,'_ID_',num2str(ID), '_',cutoff.shortname,'_',cutoff.value_str, cutoff.units];
  NameLog_out{ID_out} = OutputData;
end

end

function [System, Method, cutoff]= adjustCutoff(direction_str,System0, Method0,cutoff0,nextCutoff,uncertainty)
System = System0;
Method = Method0;
cutoff = cutoff0;
order = Method.order;
% ENUM
RADIUS  =1;  DIPOLE = 2; BAMAX = 3; POWDER = 4; 
             DIPOLE_BAMAX = 2;
switch direction_str  
  case 'relax'
    reltig = 1;
  case 'tighten'
    reltig = -1;
end

cutoff.name = nextCutoff;

switch cutoff.name
  case 'radius'
    cutoff.ID = RADIUS;
    System.radius = System.radius + reltig*cutoff.delta(cutoff.ID);
    cutoff.value = System.radius;
    cutoff.value_str = num2str(cutoff.value*1e10);
    cutoff.shortname = 'r';
    cutoff.units = 'A';
  case 'dipole'
    cutoff.ID = DIPOLE;
    Method.cutoff.dipole = Method.cutoff.dipole*10^(reltig*cutoff.delta(cutoff.ID));
    cutoff.value = Method.cutoff.dipole(1);
    cutoff.value_str = num2str(round(cutoff.value,0));
    cutoff.shortname = 'b';
    cutoff.units = 'Hz';
    
  case 'bAmax'
    cutoff.ID = BAMAX;
    Method.cutoff.bAmax = Method.cutoff.bAmax*10^(reltig*cutoff.delta(cutoff.ID) );
    cutoff.value = Method.cutoff.bAmax(1);
    cutoff.value_str = num2str(round(cutoff.value,0));
    cutoff.shortname = 'bAmax';
    cutoff.units = 'Hz';
    
  case 'pseudograd'
    cutoff.ID = DIPOLE_BAMAX;
    cutoff.shortname = 'pgrad';
    cutoff.units = 'Hz';
    
    switch Method.pseudogradType
      case 'lin_pgrad'
      % linear pseudorad steps
      Method.cutoff.dipole = Method.cutoff.dipole + reltig*cutoff.delta(cutoff.ID)*uncertainty.err_max{order}(1);
      Method.cutoff.bAmax = Method.cutoff.bAmax + reltig*cutoff.delta(cutoff.ID)*uncertainty.err_max{order}(2);
      
      case 'log_off_pgrad'
        cutoff.ID = DIPOLE_BAMAX;
        cutoff.shortname = 'pgrad';
        cutoff.units = 'Hz';
        
        % log off-pseudorad steps
        Method.cutoff.dipole = Method.cutoff.dipole*10^(log(10)*reltig*cutoff.delta(cutoff.ID)*uncertainty.err_max{order}(1));
        Method.cutoff.bAmax = Method.cutoff.bAmax*10^(log(10)*reltig*cutoff.delta(cutoff.ID)*uncertainty.err_max{order}(2));
        
      case 'log_pgrad'
        cutoff.ID = DIPOLE_BAMAX;
        cutoff.shortname = 'pgrad';
        cutoff.units = 'Hz';
        
        % log pseudorad steps
        Method.cutoff.dipole = Method.cutoff.dipole*10^(reltig*(1 + cutoff.delta(cutoff.ID)/Method.cutoff.dipole * uncertainty.err_max{order}(1)));
        Method.cutoff.bAmax = Method.cutoff.bAmax*10^(reltig*(1 + cutoff.delta(cutoff.ID)/Method.cutoff.bAmax * uncertainty.err_max{order}(2)));
        
      case 'lin_varStep'
        
        cutoff.ID = DIPOLE_BAMAX;
        cutoff.shortname = 'pgrad';
        cutoff.units = 'Hz';
        
        % linear pseudorad with variable steps
        c = reltig*cutoff.delta(cutoff.ID); % *uncertainty.err_norm(order)/min(uncertainty.err_unitPseudoGrad{order});
        if abs(c) >= 0.9
          c = sign(c)*0.9;
        end
        if abs(c) <= 0.1
          c = sign(c)*0.1;
        end
        c = c*min(Method.cutoff.dipole/uncertainty.err_unitPseudoGrad{order}(1),Method.cutoff.bAmax/uncertainty.err_unitPseudoGrad{order}(2));
        
        Method.cutoff.dipole = Method.cutoff.dipole + c*uncertainty.err_unitPseudoGrad{order}(1);
        Method.cutoff.bAmax = Method.cutoff.bAmax + c*uncertainty.err_unitPseudoGrad{order}(2);
        
        cutoff.value = abs(c);
        cutoff.value_str = num2str(round(c,0));
    end
  case 'powder'
    cutoff.ID = POWDER;
    grid_options = [1,6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, ...
      266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, ...
      3074, 3470, 3890, 4334, 4802, 5294, 5810];
    
    grid_point = find(grid_options == System.gridSize);
    System.gridSize = grid_options(grid_point + reltig*cutoff.delta(cutoff.ID));
    cutoff.value = System.gridSize;
    cutoff.value_str = num2str(cutoff.value);
    cutoff.unit = '';
    
    
end
end

