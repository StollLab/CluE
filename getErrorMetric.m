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

switch metric
  
  case 'rms'
    eta = sqrt(mean(  (v1-v2).^2   ));
  case 'root-max'
    eta = sqrt(max(  (v1-v2).^2   ));
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
    
end

end