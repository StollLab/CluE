% ========================================================================
% S- matrix
% ========================================================================
function Sminus = spinLower(spin)
Sminus = zeros(2*spin+1);
for ii = 1:(2*spin+1)-1
  Sminus(ii+1,ii) = sqrt(spin*(spin+1)-(spin - ii)*(spin+ 1 - ii));
end
end