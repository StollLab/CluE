function [outCoor, Nout] = removeOverlap(inCoor,refCoor,r0)

% Get number of reference points.
M = size(refCoor,1);

% Get number of input points.
N = size(inCoor,1);

% Initialize list of points to remove.
removeCoor  = false(N,1);

% Loop through all reference points.
for m = 1:M
  % Remove any input point that is within r0 of a reference point m.
  removeCoor  = removeCoor |  vecnorm(inCoor -refCoor(m,:),2,2) < r0(m);
end

% output all input points that were not removed.
outCoor = inCoor(~removeCoor,:);

Nout = size(outCoor,1);

end