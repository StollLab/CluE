function [parameters,System,Method,Data]  = converge_parameters(System,Method,Data, options)

options = setDefaults(options,Method);
System.gridSize = 1;
ID = 0;
hline = '-----------------------------------------------------------------';
if options.verbose
  disp(hline);
  fprintf('Starting.\n');
end
Data0 = Data;
Data.OutputData = [Data0.OutputData,'_ID_', num2str(ID)];
calculate_signal = true;
if isfile([Data.OutputData,'.mat'])
  try
    load([Data.OutputData,'.mat'],'SignalMean','experiment_time','TM_powder','Progress');
    if Progress.complete
      calculate_signal = false;
      SignalMean_ = SignalMean;
      experiment_time_ = experiment_time;
      TM_powder_ = TM_powder;
    end
  catch
  end
end
if calculate_signal
[SignalMean_, experiment_time_, TM_powder_] = nuclear_spin_diffusion(System,Method,Data);
end
if options.verbose
  fprintf('TM2 = %d s.\n',TM_powder_);
  disp(hline);
end
% System size, R ----------------------------------------------------------
are_R_r0_converged = false;
while ~are_R_r0_converged
  are_R_r0_converged = true;
  if options.converge.radius
    is_R_converged = false;
    while ~is_R_converged
      System.radius = System.radius + options.delta.radius;
      if options.verbose
        fprintf('R = %d nm.\n',System.radius*1e9);
      end
      if System.radius > options.limit.radius
        disp('System radius could did not converge within set bounds.')
        break;
      end
      
      ID = ID + 1; 
      Data.OutputData = [Data0.OutputData,'_ID_',num2str(ID), 'R_', num2str(System.radius*1e10)];
      
      calculate_signal = true;
      if isfile([Data.OutputData,'.mat'])
        try
          load([Data.OutputData,'.mat'],'SignalMean','experiment_time','TM_powder','Progress');
          if Progress.complete
            calculate_signal = false;
          end
        catch
        end
      end
      if calculate_signal
        [SignalMean, experiment_time, TM_powder] = nuclear_spin_diffusion(System,Method,Data);
      end
      
      eta = abs( (TM_powder-TM_powder_)/TM_powder);
      if options.verbose
        fprintf('TM2 = %d s.\n',TM_powder);
        fprintf('|TM2-TM1|/TM2 = %d.\n',eta);
        disp(hline);
      end
      
      if eta < options.threshold.radius
        is_R_converged = true;
        System.radius = System.radius - options.delta.radius;
      else
        are_R_r0_converged = false;
        SignalMean_ = SignalMean;
        experiment_time_ = experiment_time;
        TM_powder_ = TM_powder;
      end
    end
  end
  
  % Neighbor cutoff r0 ----------------------------------------------------
  
  if options.converge.r0
    is_neighbor_converged = false;
    
    while ~is_neighbor_converged
      Method.r0 = Method.r0 + options.delta.r0;
      if options.verbose
        fprintf('r0 = %d nm.\n',Method.r0*1e9);
      end
      
      if Method.r0 > options.limit.r0
        disp('Neighbor cutoff r0 did not converge within set bounds.')
        break;
      end
      ID = ID + 1;
      Data.OutputData = [Data0.OutputData,'_ID_',num2str(ID), '_r0_', num2str(Method.r0*1e10)];
      
      calculate_signal = true;
      if isfile([Data.OutputData,'.mat'])
        try
          load([Data.OutputData,'.mat'],'SignalMean','experiment_time','TM_powder','Progress');
          if Progress.complete
            calculate_signal = false;
          end
        catch
        end
      end
      if calculate_signal
        [SignalMean, experiment_time, TM_powder] = nuclear_spin_diffusion(System,Method,Data);
      end
      
      eta = abs( (TM_powder-TM_powder_)/TM_powder);
      if options.verbose
        fprintf('TM2 = %d s.\n',TM_powder);
        fprintf('|TM2-TM1|/TM2 = %d.\n',eta);
        disp(hline);
      end
      
      if eta < options.threshold.r0
        is_neighbor_converged = true;
        Method.r0 = Method.r0 - options.delta.r0;
      else
        are_R_r0_converged = false;
        SignalMean_ = SignalMean;
        experiment_time_ = experiment_time;
        TM_powder_ = TM_powder;
      end
    end
  end
end

% Neighbor cutoff modulation ----------------------------------------------
  
  if options.converge.modulation
    is_neighbor_converged = false;
    
    while ~is_neighbor_converged
      Method.cutoff.modulation = Method.cutoff.modulation*10^(options.delta.modulation);
      
      if options.verbose
        fprintf('modulation depth = %d.\n',Method.cutoff.modulation);
      end
      
      if Method.cutoff.modulation < 10^(options.limit.modulation)
        disp('Neighbor cutoff modulation did not converge within set bounds.')
        break;
      end
      ID = ID + 1;
      Data.OutputData = [Data0.OutputData,'_ID_',num2str(ID), '_mod_', num2str(num2str(log(Method.cutoff.modulation)/log(10)))];
      
      calculate_signal = true;
      if isfile([Data.OutputData,'.mat'])
        try
          load([Data.OutputData,'.mat'],'SignalMean','experiment_time','TM_powder','Progress');
          if Progress.complete
            calculate_signal = false;
          end
        catch
        end
      end
      if calculate_signal
        [SignalMean, experiment_time, TM_powder] = nuclear_spin_diffusion(System,Method,Data);
      end
      
      eta = abs( (TM_powder-TM_powder_)/TM_powder);
      if options.verbose
        fprintf('TM2 = %d s.\n',TM_powder);
        fprintf('|TM2-TM1|/TM2 = %d.\n',eta);
        disp(hline);
      end
      
      if eta < options.threshold.modulation
        is_neighbor_converged = true;
        Method.cutoff.modulation = Method.cutoff.modulation*10^(-options.delta.modulation);
      else
        are_R_r0_converged = false;
        SignalMean_ = SignalMean;
        experiment_time_ = experiment_time;
        TM_powder_ = TM_powder;
      end
    end
  end

  % dipole coupling  ------------------------------------------------------
  
  if options.converge.dipole
    is_neighbor_converged = false;
    
    while ~is_neighbor_converged
      Method.cutoff.dipole = Method.cutoff.dipole*10^(options.delta.dipole);
      
      if options.verbose
        fprintf('dipole coupling = %d Hz.\n',Method.cutoff.dipole);
      end
      
      if Method.cutoff.dipole < 10^(options.limit.dipole)
        disp('Neighbor cutoff dipole coupling did not converge within set bounds.')
        break;
      end
      ID = ID + 1;
      Data.OutputData = [Data0.OutputData,'_ID_',num2str(ID), '_dip_', num2str(num2str(log(Method.cutoff.dipole)/log(10)))];
      
      calculate_signal = true;
      if isfile([Data.OutputData,'.mat'])
        try
          load([Data.OutputData,'.mat'],'SignalMean','experiment_time','TM_powder','Progress');
          if Progress.complete
            calculate_signal = false;
          end
        catch
        end
      end
      if calculate_signal
        [SignalMean, experiment_time, TM_powder] = nuclear_spin_diffusion(System,Method,Data);
      end
      
      eta = abs( (TM_powder-TM_powder_)/TM_powder);
      if options.verbose
        fprintf('TM2 = %d s.\n',TM_powder);
        fprintf('|TM2-TM1|/TM2 = %d.\n',eta);
        disp(hline);
      end
      
      if eta < options.threshold.dipole
        is_neighbor_converged = true;
        Method.cutoff.dipole = Method.cutoff.dipole*10^(-options.delta.dipole);
      else
        are_R_r0_converged = false;
        SignalMean_ = SignalMean;
        experiment_time_ = experiment_time;
        TM_powder_ = TM_powder;
      end
    end
  end
  
  % hyperfine coupling infimum --------------------------------------------
  
  if options.converge.hyperfine && options.delta.hyperfine < 0
    is_neighbor_converged = false;
    
    while ~is_neighbor_converged
      
      Method.cutoff.hyperfine_inf = Method.cutoff.hyperfine_inf*10^(options.delta.hyperfine);
      
      if options.verbose
        fprintf('hyperfine coupling = %d Hz.\n',Method.cutoff.hyperfine_inf);
      end
      
      if Method.cutoff.hyperfine_inf < 10^(options.limit.hyperfine)
        disp('Neighbor cutoff hyperfine coupling did not converge within set bounds.')
        break;
      end
      ID = ID + 1;
      Data.OutputData = [Data0.OutputData,'_ID_',num2str(ID), '_hf_', num2str(num2str(log(Method.cutoff.hyperfine_inf)/log(10)))];
      
      calculate_signal = true;
      if isfile([Data.OutputData,'.mat'])
        try
          load([Data.OutputData,'.mat'],'SignalMean','experiment_time','TM_powder','Progress');
          if Progress.complete
            calculate_signal = false;
          end
        catch
        end
      end
      if calculate_signal
        [SignalMean, experiment_time, TM_powder] = nuclear_spin_diffusion(System,Method,Data);
      end
      
      eta = abs( (TM_powder-TM_powder_)/TM_powder);
      if options.verbose
        fprintf('TM2 = %d s.\n',TM_powder);
        fprintf('|TM2-TM1|/TM2 = %d.\n',eta);
        disp(hline);
      end
      
      if eta < options.threshold.hyperfine
        is_neighbor_converged = true;
        Method.cutoff.hyperfine_inf = Method.cutoff.hyperfine_inf*10^(-options.delta.hyperfine);
      else
        are_R_r0_converged = false;
        SignalMean_ = SignalMean;
        experiment_time_ = experiment_time;
        TM_powder_ = TM_powder;
      end
    end
  end 
  
   % hyperfine coupling supremum ------------------------------------------
  
  if options.converge.hyperfine && options.delta.hyperfine > 0
    is_neighbor_converged = false;
    
    while ~is_neighbor_converged
      
      Method.cutoff.hyperfine_sup = Method.cutoff.hyperfine_sup*10^(options.delta.hyperfine);
      
      if options.verbose
        fprintf('hyperfine coupling = %d Hz.\n',Method.cutoff.hyperfine_sup);
      end
      
      if Method.cutoff.hyperfine_sup > 10^(options.limit.hyperfine)
        disp('Neighbor cutoff hyperfine coupling did not converge within set bounds.')
        break;
      end
      ID = ID + 1;
      Data.OutputData = [Data0.OutputData,'_ID_',num2str(ID), '_hf_', num2str(num2str(log(Method.cutoff.hyperfine_sup)/log(10)))];
      
      calculate_signal = true;
      if isfile([Data.OutputData,'.mat'])
        try
          load([Data.OutputData,'.mat'],'SignalMean','experiment_time','TM_powder','Progress');
          if Progress.complete
            calculate_signal = false;
          end
        catch
        end
      end
      if calculate_signal
        [SignalMean, experiment_time, TM_powder] = nuclear_spin_diffusion(System,Method,Data);
      end
      
      eta = abs( (TM_powder-TM_powder_)/TM_powder);
      if options.verbose
        fprintf('TM2 = %d s.\n',TM_powder);
        fprintf('|TM2-TM1|/TM2 = %d.\n',eta);
        disp(hline);
      end
      
      if eta < options.threshold.hyperfine
        is_neighbor_converged = true;
        Method.cutoff.hyperfine_sup = Method.cutoff.hyperfine_sup*10^(-options.delta.hyperfine);
      else
        are_R_r0_converged = false;
        SignalMean_ = SignalMean;
        experiment_time_ = experiment_time;
        TM_powder_ = TM_powder;
      end
    end
  end 
  
  % minimum frequency -----------------------------------------------------
  if options.converge.minimum_frequency
    is_neighbor_converged = false;
    
    while ~is_neighbor_converged
      Method.cutoff.minimum_frequency = Method.cutoff.minimum_frequency*10^(options.delta.minimum_frequency);
      
      if options.verbose
        fprintf('minimum frequency coupling = %d Hz.\n',Method.cutoff.minimum_frequency);
      end
      
      if Method.cutoff.minimum_frequency < 10^(options.limit.minimum_frequency)
        disp('Neighbor cutoff minimum frequency coupling did not converge within set bounds.')
        break;
      end
      ID = ID + 1;
      Data.OutputData = [Data0.OutputData,'_ID_',num2str(ID), '_freq_', num2str(num2str(log(Method.cutoff.minimum_frequency)/log(10)))];
      
      calculate_signal = true;
      if isfile([Data.OutputData,'.mat'])
        try
          load([Data.OutputData,'.mat'],'SignalMean','experiment_time','TM_powder','Progress');
          if Progress.complete
            calculate_signal = false;
          end
        catch
        end
      end
      if calculate_signal
        [SignalMean, experiment_time, TM_powder] = nuclear_spin_diffusion(System,Method,Data);
      end
      
      eta = abs( (TM_powder-TM_powder_)/TM_powder);
      if options.verbose
        fprintf('TM2 = %d s.\n',TM_powder);
        fprintf('|TM2-TM1|/TM2 = %d.\n',eta);
        disp(hline);
      end
      
      if eta < options.threshold.minimum_frequency
        is_neighbor_converged = true;
        Method.cutoff.minimum_frequency = Method.cutoff.minimum_frequency*10^(-options.delta.minimum_frequency);
      else
        are_R_r0_converged = false;
        SignalMean_ = SignalMean;
        experiment_time_ = experiment_time;
        TM_powder_ = TM_powder;
      end
    end
  end
  

% Powder orientations
grid_options = [1,6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, ...
  266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, ...
  3074, 3470, 3890, 4334, 4802, 5294, 5810];

limit.grid_index = length(grid_options);
is_grid_converged = false;
grid_point = 1;
if options.converge.powder
  while ~is_grid_converged
    grid_point = grid_point + 1;
    
    if grid_point > limit.grid_index
      disp('Powder grid could not converge.')
      break;
    end
    
    System.gridSize = grid_options(grid_point);
    
    if options.verbose
      fprintf('grid = %d points.\n',System.gridSize);
    end
    
    if System.gridSize > options.limit.grid_points
      disp('Powder grid did not converge within set bounds.')
      break;
    end
    
    
    ID = ID + 1;
    
    Data.OutputData = [Data0.OutputData,'_ID_',num2str(ID), '_powder_', num2str(System.gridSize)];

    calculate_signal = true;
      if isfile([Data.OutputData,'.mat'])
        try
          load([Data.OutputData,'.mat'],'SignalMean','experiment_time','TM_powder','Progress');
          if Progress.complete
            calculate_signal = false;
          end
        catch
        end
      end
      if calculate_signal
        [SignalMean, experiment_time, TM_powder] = nuclear_spin_diffusion(System,Method,Data);
      end
    
    eta = abs( (TM_powder-TM_powder_)/TM_powder);
    if options.verbose
      
      fprintf('TM2 = %d s.\n',TM_powder);
      fprintf('|TM2-TM1|/TM2 = %d.\n',eta);
      disp(hline);
      
    end
    
    if eta < options.threshold.powder
      is_grid_converged = true;
      grid_point = grid_point - 1;
      System.gridSize = grid_options(grid_point);
    else
      SignalMean_ = SignalMean;
      experiment_time_ = experiment_time;
      TM_powder_ = TM_powder;
    end
    
  end
end

parameters = [];
if options.converge.radius
  parameters.radius = System.radius;
end
if options.converge.r0
parameters.r0 = Method.r0;
end
if options.converge.modulation
  parameters.modulation = Method.cutoff.modulation;
end
if options.converge.dipole
parameters.dipole = Method.cutoff.dipole;
end
if options.converge.powder
  parameters.gridSize = System.gridSize;
end
end

function options = setDefaults(options,Method)
% system radius
if ~isfield(options.converge,'radius')
  options.converge.radius = true;
end
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

if false % isfield(Method,'Criteria')
num_criteria = numel(Method.Criteria);  
  for ii = 1:num_criteria
    if strcmp(Method.Criteria{ii},'neighbor')
      options.converge.r0 = true;
    elseif strcmp(Method.Criteria{ii},'modulation')
      options.converge.modulation = true;
    elseif strcmp(Method.Criteria{ii},'psuedo-secular')
      options.converge.dipole = true;
    elseif strcmp(Method.Criteria{ii},'minimum-frequency')
      options.converge.minimum_frequency = true;
    end
  end
end

% distance
if ~isfield(options.converge,'r0')
  options.converge.r0 = false;
end
if ~isfield(options.threshold,'r0')
  options.threshold.r0 = 1e-3;
end
if ~isfield(options.delta,'r0')
  options.delta.r0 = 1e-10; % m
end
if ~isfield(options.limit,'r0')
  options.limit.r0 = 20e-10; % m
end

% modulation depth
if ~isfield(options.converge,'modulation')
  options.converge.modulation = false;
end
if ~isfield(options.threshold,'modulation')
  options.threshold.modulation = 1e-3;
end
if ~isfield(options.delta,'modulation')
  options.delta.modulation = -1; 
end
if ~isfield(options.limit,'modulation')
  options.limit.modulation = -15; 
end

% dipole
if ~isfield(options.converge,'dipole')
  options.converge.dipole = false;
end
if ~isfield(options.threshold,'dipole')
  options.threshold.dipole = 1e-3;
end
if ~isfield(options.delta,'dipole')
  options.delta.dipole = -1; 
end
if ~isfield(options.limit,'dipole')
  options.limit.dipole = -15; 
end

% hyperfine
if ~isfield(options.converge,'hyperfine')
  options.converge.hyperfine = false;
end
if ~isfield(options.threshold,'hyperfine')
  options.threshold.hyperfine = 1e-3;
end
if ~isfield(options.delta,'hyperfine')
  options.delta.hyperfine = -1; 
end
if ~isfield(options.limit,'hyperfine')
  options.limit.hyperfine = -15; 
end

% minimum_frequency
if ~isfield(options.converge,'minimum_frequency')
  options.converge.minimum_frequency = false;
end
if ~isfield(options.threshold,'r0')
  options.threshold.minimum_frequency = 1e-3;
end
if ~isfield(options.delta,'r0')
  options.delta.minimum_frequency = -1;
end
if ~isfield(options.limit,'r0')
  options.limit.minimum_frequency = -15;
end


% powder
if ~isfield(options.converge,'powder')
  options.converge.powder = true;
end
if ~isfield(options.threshold,'powder')
  options.threshold.powder = 1e-2;
end
if ~isfield(options.limit,'grid_points')
  options.limit.grid_points = 5810;
end
% verbosity
if ~isfield(options,'verbose')
  options.verbose = true;
end

end





