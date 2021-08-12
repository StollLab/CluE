function [ishermitian,nonHermiticity] = isHermitian(H,threshold)
if numel(H) == 1 && H==0
  ishermitian = true;
  nonHermiticity = 0;
  return;
end
mma = @(A)max(max(abs(A)));

% Check for the zero matrix.
if mma(H) < threshold
  nonHermiticity = mma(H);
else
  nonHermiticity = mma(H-H')/mma(H);
end

ishermitian = nonHermiticity<=threshold;
end
