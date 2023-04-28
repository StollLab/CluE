function R = rotate_theta_phi(theta,phi)


cos_phi = cos(phi);
sin_phi = sin(phi);

cos_theta = cos(theta);
sin_theta = sin(theta);

R = [cos_phi*cos_theta, -sin_phi, cos_phi*sin_theta;...
     sin_phi*cos_theta,  cos_phi, sin_phi*sin_theta;...
             -sin_theta,         0,         cos_theta];

end