% Given a normalized x_axis and z_axis, R = alignCoordinates(x_axis,z_axis) is 
% an active rotation matrix R*x_axis = [1;0;0], and R*z_axis = [0;0;1].
function R = alignCoordinates(x_axis,z_axis)

assert( numel(x_axis)==3 );
assert( numel(z_axis)==3 );

assert(isreal(x_axis));
assert(isreal(z_axis));

assert(norm(x_axis)>0);
assert(norm(z_axis)>0);

normalized_z_axis = z_axis(:)/norm(z_axis);
normalized_x_axis = x_axis(:)/norm(x_axis);

assert( abs(normalized_x_axis'*normalized_z_axis) < 1)

y_axis = cross(normalized_z_axis,normalized_x_axis); 
normalized_y_axis = y_axis(:)/norm(y_axis);

R = [1;0;0]*normalized_x_axis' ...
  + [0;1;0]*normalized_y_axis'...
  + [0;0;1]*normalized_z_axis';

end