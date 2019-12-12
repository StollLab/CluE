% ========================================================================
% S+ matrix
% ========================================================================
function Splus = spinRaise(spin)
Splus = zeros(2*spin+1);
for ii = 1:(2*spin+1)-1
  Splus(ii,ii+1) = sqrt(spin*(spin+1)-(spin - ii)*(spin+ 1 - ii));
end
end