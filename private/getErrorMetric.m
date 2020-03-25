function eta = getErrorMetric(v1,v2,metric,t1,t2)

if nargin < 3
  metric = 'rms';
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