% ========================================================================
% Rotation matrix around z axis by angle phi
% ========================================================================

function Rz = rotateZ(phi)
Rz  =   [cos(phi ), -sin(phi ), 0; ...
  sin(phi ),  cos(phi ), 0; ...
  0,          0,         1];
end