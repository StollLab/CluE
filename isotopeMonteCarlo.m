function [signals,TM] = isotopeMonteCarlo(System,Method,Data, savefile,N,dN,threshold,options)

if (mod(dN,2) ~= 0) || (dN < 0)
  error('The parameter dN must be an even natural number.')
end
if (mod(N,2) ~= 0) || (N < 0)
  error('The parameter N must be an even natural number.')
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
primeRange = [2^10,2^20];



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
if isfile(savefile)
  try
    load(savefile,'signals','progress','TM','N','twotau');
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
  for ii=1:N
    
    % Skip loaded trials.
    if signals(ii,1)==0
      if ~isempty(options.seed)
        Method.seed = nthprime(randi(primeRange));
        fprintf('Setting Method.seed to %d.\n',Method.seed)
      end
      
      fprintf('Running inital trial %d/%d.\n', ii,N);
      [signals(ii,:),twotau,TM(ii)] = CluE(System,Method,Data);
      
      saveCounter = saveCounter + 1;
      
      % Save in case run is canceled.
      if saveCounter >= saveEveryN
        saveCounter = 0;
        save(savefile);
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
  save(savefile);
end

saveCounter = 0;
conNum = 1;

% Check if convergence runs are needed.
if ~progress(CONVERGENCE_TRIALS)
  
  isConverged = false;
  
  % Run convergence loop.
  while ~isConverged
    
    % Loop over new trials.
    for ii=N+1:(N+dN)
      
      % Skip loaded trials
      if signals(ii,1)==0
        
        if ~isempty(options.seed)
          Method.seed = nthprime(randi(primeRange));
          fprintf('Setting Method.seed to %d.\n',Method.seed)
        end
        
        fprintf('Running convergene trial %d: %d/%d.\n',conNum, ii,N+dN);
        [signals(ii,:),twotau,TM(ii)] = CluE(System,Method,Data);
        
        saveCounter = saveCounter + 1;
        
        % Save in case run is canceled.
        if saveCounter >= saveEveryN
          saveCounter = 0;
          save(savefile);
        end
      end
      
      % Find the overall mean signal, and TM.
      v3 = mean(signals(1:ii,:),1);
      TM(ii) = getTM(twotau,v3);
      
      fprintf('TM  = %d us.\n',TM(ii)*1e6);
      
    end
    
    % Partition trials into two qual parts.
    N_ = (N+dN);
    
    
    % Get measure of difference.
    eta = 0;
    for ii = 1:N_ave
      perm = randperm(N_);
      
      % Find the mean signal of each partition.
      v1 = mean(  signals(  perm( 1:N_/2    ) ,: ) , 1);
      v2 = mean(  signals(  perm( N_/2+1:N_ ) ,: ) , 1);
      
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
        save(savefile);
        return;
      end
      
      % Update N.
      N = N + dN;
      
      % Initialize memory.
      signals(N+dN,:) = zeros(1,System.timepoints);
      TM(N+dN) = 0;
      
      conNum = conNum + 1;
    end
    
    % Save.
    save(savefile);
  end
end
end