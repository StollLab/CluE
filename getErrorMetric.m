function eta = getErrorMetric(v1,v2,metric,t1,t2,options)

if nargin < 3
  metric = 'rms';
end
if nargin < 6
  options = [];
end

if isfield(options,'tmax')
  v1 = v1( t1 <= options.tmax);
  t1 = t1( t1 <= options.tmax);
  v2 = v2( t2 <= options.tmax);
  t2 = t2( t2 <= options.tmax);
end

if isfield(options,'vmin')
  v_=max(abs(v1),abs(v2));
  to_remove = v_<options.vmin;
  v1(to_remove) = [];
  t1(to_remove) = [];
  v2(to_remove) = [];
  t2(to_remove) = [];
end

if any(isnan(v1)) || any(isnan(v1))
  disp('Warning, removing NaN entries');
  vnan = isnan(v1) | isnan(v2);
  v1(vnan) = [];
  v2(vnan) = [];
  t1(vnan) = [];
  t2(vnan) = [];
  if isempty(v1)
    disp('Warning, no non-NaN entries to compares.');
    eta = nan;
    return;
  else
    fprintf('Removed %d NaN entries.\n',sum(vnan));
    fprintf('Calculating error metric with the %d non-NaN entries.\n',sum(~vnan));
  end
  
end


switch metric
  
  case 'rms'
    eta = sqrt(mean(  abs(v1-v2).^2   ));
  case 'root-max'
    eta = sqrt(max(  abs(v1-v2).^2   ));
  case 'max-abs'
    eta = max(abs(v1-v2));
  case 'mean-abs'
    eta = mean(abs(v1-v2));
  case '1/e'
    TM1 = getTM(t1,v1);
    TM2 = getTM(t2,v2);
    eta = abs( 2*(TM1-TM2)/(TM1+TM2) );
  case '1/e 2'
    TM1 = getTM(t1,v1);
    TM2 = getTM(t2,v2);
    eta = abs( (TM1-TM2)/TM2);
  case '1/e 1'
    TM1 = getTM(t1,v1);
    TM2 = getTM(t2,v2);
    eta = abs( (TM1-TM2)/TM1);

  case 'max(rms,1/e)'   
    eta_rms = sqrt(mean(  abs(v1-v2).^2   ));
    TM1 = getTM(t1,v1);
    TM2 = getTM(t2,v2);
    eta_e = abs( 2*(TM1-TM2)/(TM1+TM2) );
    eta = max(eta_rms,eta_e);
end

end
