function [signals,TM] = isotopeMonteCarlo(System,Method,Data, savefile,N,threshold,options)

if ~isfield(options,'saveEveryN')
  options.saveEveryN = 16;
end
if ~isfield(options,'maxN')
  options.maxN = inf;
end
if ~isfield(options,'newProgress')
  options.newProgress = [];
end

saveEveryN = options.saveEveryN;
maxN = options.maxN;
newProgress = options.newProgress;

saveCounter = 0;

doReset = true; 
if isfile(savefile)
  try
    load(savefile,'signals','progress','TM','N','twotau');   
    disp('Save data loaded.');
    
    doReset = false;
  catch
  end
end

if doReset
  disp('Initializing MonteCarlo.');
  progress = false(1,2);
  signals = zeros(2*N,System.timepoints);
  TM = zeros(1,2*N);
end

INITIAL_TRIALS = 1;
CONVERGENCE_TRIALS = 2;

if ~isempty(newProgress)
  progress = newProgress;
end

if ~progress(INITIAL_TRIALS)
  N = min(N,maxN);
  for ii=1:N
    
    if signals(ii,1)==0
      fprintf('Running inital trail %d/%d.\n', ii,N);
      [signals(ii,:),twotau,TM(ii)] = CluE(System,Method,Data);
      saveCounter = saveCounter + 1;
      if saveCounter >= saveEveryN
        saveCounter = 0;
        save(savefile);
      end
    end
    
  end
  
  progress(INITIAL_TRIALS) = true;
  save(savefile);
end

saveCounter = 0;
conNum = 1;
if ~progress(CONVERGENCE_TRIALS)
  isConverged = false;
  while ~isConverged
    for ii=N+1:2*N
      if signals(ii,1)==0
        fprintf('Running convergene trail %d: %d/%d.\n',conNum, ii,2*N);
        [signals(ii,:),twotau,TM(ii)] = CluE(System,Method,Data);
        saveCounter = saveCounter + 1;
        if saveCounter >= saveEveryN
          saveCounter = 0;
          save(savefile);
        end
      end
      
    end
    v1 = mean(signals(1:N,:),1);
    v2 = mean(signals(N+1:2*N,:),1);
    v3 = mean(signals,1);
    TM1 = getTM(twotau,v1);
    TM2 = getTM(twotau,v2);
    TM3 = getTM(twotau,v3);
    eta = abs((TM2-TM1)/TM3);
    fprintf('TM  = %d us.\n',TM3*1e6);
    fprintf('eta  = %d.\n',eta);
      
    if eta < threshold
      isConverged = true;
      progress(CONVERGENCE_TRIALS) = true;
    else
      
      if 2*N > maxN
        save(savefile);
        return;
      end
      
      N = 2*N;
      signals(2*N,:) = zeros(1,System.timepoints);
      TM(2*N) = 0;
      conNum = conNum + 1;
    end
    save(savefile);
  end
end
end