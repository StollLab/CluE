% ========================================================================
% Rotation matrix around x axis by angle gamma
% ========================================================================

function Rx = rotateX(gamma)
Rx  =   [1,  0,            0; ...
  0,  cos(gamma),  -sin(gamma); ...
  0,  sin(gamma),   cos(gamma )];
end