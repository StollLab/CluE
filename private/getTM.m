% ========================================================================
% Find decay time (where signal drops to 1/e of initial amplitude)
% ========================================================================
function TM = getTM(t,Signal)

% normalize signal maximum to 1
Signal = abs(Signal/Signal(1));

if any(isnan(Signal))
  TM = inf;
  return
end

e_inv = exp(-1);
if min(Signal) > e_inv
  TM = inf;
  return
end

% translate signal to put TM at the x-intercept
searchVec = abs(Signal) - e_inv;
negatives = find(searchVec<=0);
if isempty(negatives)
  TM = inf;
  return
end

% initial search index
iTM = negatives(1);
if length(iTM)>1
  iTM = iTM(1);
end
% number of time points
N = length(t);

% initial TM guess
TM = t(iTM);

% search
searching = true;
while searching
  
  % get TM guess
  v1 = Signal(iTM);
  i_mismatch = 0;
  
  % examine guess
  if v1 == e_inv
    
    % TM found
    searching = false;
    break;
    
  elseif v1 > e_inv
    
    % TM guess too soon
    % make second guess
    jTM = iTM + 1;
    i_mismatch = +1;
    
  elseif v1 < e_inv
    
    % TM guess too late
    % make second guess
    jTM = iTM - 1;
    i_mismatch = -1;
    
  end
  
  % checkguess validity
  if jTM<1
    
    % guess outside range of data
    TM = NaN;
    return;
    
  end
  if jTM>N
    
    % guess outside range of data
    TM = inf;
    return;
    
  end
  
  % new guess
  v2 = Signal(jTM);
  jmismatch = 0;
  
  % examine new guess
  if v2 == e_inv
    
    % TM found
    TM = t(jTM);
    iTM = jTM;
    v1 = v2;
    searching = false;
    break;
  elseif v1 > e_inv
    
    % TM guess too soon
    jmismatch = +1;
    
  elseif v2 < e_inv
    
    % TM guess too late
    jmismatch = -1;
    
  end
  
  % error parameter
  mismatch = i_mismatch*jmismatch;
  
  %decide how to make another guess
  if mismatch > 0
    
    % shift over one index
    TM = t(jTM);
    iTM = jTM;
    v1=v2;
    
  else
    
    % one guess is too late, the other too soon
    searching = false;
    
    % interpolate TM between the two guesses
    d2tau = t(jTM) - t(iTM);
    dv   =  v2-v1;
    TM = (d2tau/dv)*(e_inv - v1 + t(iTM)*dv/d2tau);
    
  end
  
end
end