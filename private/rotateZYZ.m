% ========================================================================
% Convert Euler angles to rotation matrix
% ========================================================================

function R = rotateZYZ(alpha,beta,gamma)
sa = sin(alpha);
ca = cos(alpha);
sb = sin(beta);
cb = cos(beta);
sc = sin(gamma);
cc = cos(gamma);
R = [ cc*cb*ca - sc*sa,  cc*cb*sa + sc*ca, -cc*sb; ...
     -sc*cb*ca - cc*sa, -sc*cb*sa + cc*ca,  sc*sb; ...
         sb*ca,             sb*sa,             cb];
end
