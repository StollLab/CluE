%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function tests = test_rotate()
tests = functiontests(localfunctions);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function test_rotate_theta_phi(testCase)

theta_list = linspace(0,pi,11);
phi_list = linspace(0,2*pi,11);
z = [0;0;1];

for theta = theta_list
  for phi= phi_list

    R = rotate_theta_phi(theta,phi);
    v = R*z;

    ref = [cos(phi)*sin(theta); sin(phi)*sin(theta); cos(theta)];

    assert(all(v==ref));

  end
end
end
%-------------------------------------------------------------------------------
function test_alignCoordinates(testCase)
  x = [1;0;0];
  z = [0;0;1];
  R = alignCoordinates(x,z);
  assert( all(all( R==eye(3) )) );

  x = [0;0;1];
  z = [1;0;0];
  R = alignCoordinates(x,z);
  ref = [0,0,1;...
         0,-1,0;...
         1,0,0];
  assert( max(max(abs( R - ref ))) < 1e-12 );


  A = rand(3);
  A = A + A';
  while abs(det(A))<1e-1
    A = rand(3);
    A = A + A';
  end
  A = A/det(A);

  [evec,~] = eig(A);
  evec = evec/det(evec);
  R = alignCoordinates(evec(:,1),evec(:,3));
  assert( max(max(abs( R*R' - eye(3) ))) < 1e-12 );
  assert( max(max(abs( R*evec - eye(3) ))) < 1e-12 );

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
