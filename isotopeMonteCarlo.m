function [signals,TM,statistics,twotau] = isotopeMonteCarlo(System,Method,Data, savefile,N,dN,threshold,options)

if (dN < 1)
  error('The parameter dN >= 1.')
end
if (N < 1)
  error('The parameter N >= 1.')
end


if ~isfield(options,'saveEveryN')
  options.saveEveryN = dN;
end
if ~isfield(options,'maxN')
  options.maxN = inf;
end
if ~isfield(options,'N_ave')
  options.N_ave = 1000;
end

if ~isfield(options,'newProgress')
  options.newProgress = [];
end
if ~isfield(options,'metric')
  options.metric = 'rms';
end
if ~isfield(options,'vmin')
  options.vmin = exp(-5);
end
if ~isfield(options,'seed')
  options.seed = 42;
end
if ~isempty(options.seed)
  fprintf('Using rng seed to options.seed = %d.\n',options.seed)
  rng(options.seed);
end
if ~isfield(options,'parallelComputing')
  options.parallelComputing = false;
end
if isfield(Method,'parallelComputing') && Method.parallelComputing && options.parallelComputing
  error('Method.parallelComputing & options.parallelComputing cannot both be true.')
end

if ~isfield(options,'numCores')
  options.numCores = feature('numcores');
end
if options.parallelComputing
  delete(gcp('nocreate'));
  pool = parpool(options.numCores);
end

if ~isfield(options,'conserveMemory')
  options.conserveMemory = true;
end
if ~isfield(options,'concentration_1H')
  options.concentration_1H = true;
end

% Set RNG to ensure that repeats on the same initial conditions
% give identical results. 
iMax = 2^20;
try 
  rng(options.seed)
catch 
  options.seed = 42;
  rng(options.seed)
end
nextSeed = randi(iMax);



if ~isfield(System,'randomOrientation')
  System.randomOrientation = true;
end

saveEveryN = options.saveEveryN;
maxN = options.maxN;
N_ave = options.N_ave;
newProgress = options.newProgress;

saveCounter = 0;

doReset = true;

% Try to continue from canceled run.
if isfile(savefile) && ~(isfield(Data,'overwriteLevel') && Data.overwriteLevel ==2 )
  try
    load(savefile,'signals','progress','TM','N','twotau','statistics');
    disp('Save data loaded.');
    
    doReset = false;
  catch
  end
end

% Start from scratch.
if doReset
  disp('Initializing MonteCarlo.');
  progress = false(1,2);
  signals = zeros(2*N,System.timepoints);
  TM = zeros(1,2*N);
  statistics = cell(1,2*N);
 end

% ENUM
INITIAL_TRIALS = 1;  CONVERGENCE_TRIALS = 2;

% Set progress flags.
if ~isempty(newProgress)
  progress = newProgress;
end

% Decide if initial trials already been loaded.
if ~progress(INITIAL_TRIALS)
  
  % Determine number of trials to run.
  N = min(N,maxN);
  
  % Run initial set of trials.
  iiStart = find( signals(:,1)==0,1);
  if isempty(iiStart)
    iiStart = N + 1;
  end
  
  rng(nextSeed(1));
  Method.seed = nextSeed(1);
  nextSeed = randi(iMax,1,N);
  indices = iiStart:N;
  N_ = length(indices);
  
  if options.parallelComputing
    % Skip loaded trials.
    if ~all(signals(indices,1)==1)
      temp_signals = signals(indices,:);
      temp_twotau = -1*ones(size(signals(indices,:)));
      temp_TM = TM(indices);
      temp_statistics = {statistics{indices}};
      parfor ii=1:N_
        
        Method_ = Method;
        Method_.seed = nextSeed(ii);
        
        fprintf('Setting Method.seed to %d.\n',Method_.seed)
        
        fprintf('Running initial trial %d/%d.\n', indices(ii),N);
        [temp_signals(ii,:),temp_twotau(ii,:),temp_TM(ii),~,~,temp_statistics{ii}] = CluE(System,Method_,Data);
        
      end
      
      signals(indices,:) = temp_signals;
      TM(indices) = temp_TM;
      for ii=1:dN
        statistics{indices(ii)} = temp_statistics{ii};
        if temp_twotau(ii,1)~= -1
          twotau = temp_twotau(ii,:);
        end
      end
    end
  else
    for ind_ = 1:N_
      ii = indices(ind_);
      
      % Skip loaded trials.
      if signals(ii,1)==0
        
        Method.seed = nextSeed(ind_);
        
        fprintf('Setting Method.seed to %d.\n',Method.seed)
        
        fprintf('Running initial trial %d/%d.\n', ii,N);
        [signals(ii,:),twotau,TM(ii),~,~,statistics{ii}] = CluE(System,Method,Data);
        
        saveCounter = saveCounter + 1;
        
        % Save in case run is canceled.
        if saveCounter >= saveEveryN
          saveCounter = 0;
          save(savefile,'-v7.3');
        end
        
      end
    end
    % Find the overall mean signal, and TM.
    v3 = mean(signals(1:ii,:),1);
    TM(ii) = getTM(twotau,v3);
    
    fprintf('TM  = %d us.\n',TM(ii)*1e6);
    
    
  end
  
  % Update progress flag.
  progress(INITIAL_TRIALS) = true;
  
  % Save.
  save(savefile,'-v7.3');
end

saveCounter = 0;
conNum = 1;

% Check if convergence runs are needed.
if ~progress(CONVERGENCE_TRIALS)
  
  isConverged = false;
  
  % Run convergence loop.
  while ~isConverged
    
    
    rng(nextSeed(1));
    Method.seed = nextSeed(1);
    nextSeed = randi(iMax,1,dN);
    
    % Loop over new trials.
    indices = N+1:(N+dN);
    
    if options.parallelComputing
      % Skip loaded trials.
      if ~all(signals(indices,1)==1)
        
        
        temp_signals = signals(indices,:);
        temp_TM = TM(indices);
        temp_statistics = {statistics{indices}};
        parfor ii=1:dN
          
          Method_ = Method;
          Method_.seed = nextSeed(ii);
          
          fprintf('Setting Method.seed to %d.\n',Method_.seed)
          
          
          fprintf('Running convergene trial %d: %d/%d.\n',conNum, ii,N+dN);
          [temp_signals(ii,:),~,temp_TM(ii),~,~,temp_statistics{ii}] = CluE(System,Method_,Data);
          
        end
        
        signals(indices,:) = temp_signals;
        TM(indices) = temp_TM;
        for ii=1:dN
          statistics{indices(ii)} = temp_statistics{ii};
        end
      end
    else
      for ind_ = 1:dN
        ii = indices(ind_);
        
        % Skip loaded trials
        if signals(ii,1)==0
          
          Method.seed = nextSeed(ind_);
          
          fprintf('Setting Method.seed to %d.\n',Method.seed)
          
          
          fprintf('Running convergene trial %d: %d/%d.\n',conNum, ii,N+dN);
          [signals(ii,:),~,TM(ii),~,~,statistics{ii}] = CluE(System,Method,Data);
          
          saveCounter = saveCounter + 1;
          
          % Save in case run is canceled.
          if saveCounter >= saveEveryN
            saveCounter = 0;
            save(savefile,'-v7.3');
          end
        end
      end
    end
    % Find the overall mean signal, and TM.
    v3 = mean(signals(1:ii,:),1);
    TM(ii) = getTM(twotau,v3);
    
    fprintf('TM  = %d us.\n',TM(ii)*1e6);
    save(savefile,'-v7.3');
    % Partition trials into two qual parts.
    N_ = (N+dN);
    
    
    % Get measure of difference.
    eta = 0;
    for ii = 1:N_ave
      perm = randperm(N_);
      
      % Find the mean signal of each partition.
      v1 = mean(  signals(  perm( 1:ceil(N_/2)    ) ,: ) , 1);
      v2 = mean(  signals(  perm( ( ceil(N_/2) + 1 ) :N_ ) ,: ) , 1);
      
      eta = eta + getErrorMetric(v1,v2,options.metric,twotau,twotau,options)/N_ave;
    end
    fprintf('eta  = %d.\n',eta);
    
    % Compare measure of difference to the set threshold.
    if eta < threshold
      isConverged = true;
      progress(CONVERGENCE_TRIALS) = true;
    else
      
      % Check that another round will not exceed limit..
      if N + 2*dN > maxN
        save(savefile,'-v7.3');
        return;
      end
      
      % Update N.
      N = N + dN;
      
      % Initialize memory.
      signals(N+dN,:) = 0;
      TM(N+dN) = 0;
      statistics{N+dN} = {};
      
      conNum = conNum + 1;
    end
    
    % Save.
    save(savefile,'-v7.3');
  end
end
if exist('v1')
  clear('v1');
end
if exist('v2')
  clear('v2');
end
if exist('v3')
  clear('v3');
end
if options.concentration_1H
  number_samples = 0;
  number_1H_exchangeable = 0;
  number_1H_nonExchangeable = 0;
  number_2H_exchangeable = 0;
  number_2H_nonExchangeable = 0;
  for ii = 1:numel(statistics)
    if isfield(statistics{ii},'Statistics')
      for jj = 1:numel(statistics{ii}.Statistics)
        if isfield(statistics{ii}.Statistics{jj},'isotopologueStatistics')
          iStat = statistics{ii}.Statistics{jj}.isotopologueStatistics;
          number_samples = number_samples + 1;
          number_1H_exchangeable = number_1H_exchangeable ...
            + iStat.number_1H_exchangeable;
          number_1H_nonExchangeable = number_1H_nonExchangeable ...
            + iStat.number_1H_nonExchangeable;
          number_2H_exchangeable = number_2H_exchangeable ...
            + iStat.number_2H_exchangeable;
          number_2H_nonExchangeable = number_2H_nonExchangeable ...
            + iStat.number_2H_nonExchangeable;

        end
        if options.conserveMemory
          if isfield(statistics{ii}.Statistics{jj},'Nuclear_Dipole')
            statistics{ii}.Statistics{jj} = ...
              rmfield(statistics{ii}.Statistics{jj},'Nuclear_Dipole');
          end
          if isfield(statistics{ii}.Statistics{jj},'Nuclear_Dipole_x_iy_Z')
            statistics{ii}.Statistics{jj} = ...
              rmfield(statistics{ii}.Statistics{jj},'Nuclear_Dipole_x_iy_Z');
          end
        end
      end
    end
  end

  mean_number_1H = (number_1H_exchangeable + number_1H_nonExchangeable)...
    /number_samples;
  mean_number_2H = (number_2H_exchangeable + number_2H_nonExchangeable)...
    /number_samples;
  volume = 4*pi/3*System.radius^3;
  concentration_1H_number_per_m3 = mean_number_1H/volume;
  avogadro= 6.022140857e23;
  concentration_1H_Molar = concentration_1H_number_per_m3/avogadro/1000; 
  %   mean_number_H = number_1H + number_2H;
  %   mean_proton_fraction = number_1H/number_H;
  %   mean_deuteron_fraction = number_2H/number_H;
end
SignalMean = mean(signals,1);
SignalMean = SignalMean/SignalMean(1);
if options.parallelComputing
  delete(pool)
end
save(savefile,'-v7.3');
end