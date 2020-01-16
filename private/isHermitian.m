function [ishermitian,nonHermiticity] = isHermitian(H,threshold)
if numel(H) == 1 && H==0
  ishermitian = true;
  nonHermiticity = nan;
  return;
end
mma = @(A)max(max(abs(A)));
nonHermiticity = mma(H-H')/mma(H);
ishermitian = nonHermiticity<=threshold;
end
