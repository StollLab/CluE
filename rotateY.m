% ========================================================================
% Rotation matrix around y axis by angle theta
% ========================================================================

function Ry = rotateY(theta)
Ry =   [ cos(theta),  0,         sin(theta); ...
  0,           1,         0;          ...
  -sin(theta),  0,         cos(theta)];
end
