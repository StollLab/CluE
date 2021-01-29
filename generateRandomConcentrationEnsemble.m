% The function generateRandomConcentrationEnsemble(concentration,R,R0,r0) generates
% a set of coordinates with concentration in molarity and volume 4/3*pi*R^3, 
% such that all points are at least R0 [=] m from the origin, and 
% r0 [=] m from each other.
%
%
% concentration [=] M.  R,R0,r0 [=] m.
% coor,DistanceMatrix,coorDistance [=] m.  
function [coor,N,DistanceMatrix,coorDistance] = generateRandomConcentrationEnsemble(concentration,R,R0,d0)

% Define conversion factors.
NA = 6.02214076e23; % exact,2021-01-21, https://physics.nist.gov/cgi-bin/cuu/Value?na
litre_to_m3 = 1e-3; % m^3;

% number density
C = concentration*NA/litre_to_m3; % particle per m^3/

% volume
V = 4/3*pi*R^3;

% number of particles
N = round(C*V);

% excludedVolume = 4/3*pi*(R0^3 + N*r0^3); 
% disp(excludedVolume/V);

% Initialize coordinates in a cube with sides of length 2R,
% centered on the origin.
coor = R*(2*rand(N,3) -1);

% particle-particle distances
DistanceMatrix = zeros(N);
counter = 0;
listSize = 2^10;
n_list = zeros(1,listSize);
doCalc = true;
while doCalc
  
  % Get distances from the origen.
  coorDistance = vecnorm(coor,2,2);
  
  % Recalculate any points outside of the sphere of radius R.
  recalc = coorDistance > R;
  
  % Recalculate any points within the sphere of radius R0.
  recalc = recalc | coorDistance < R0;
  

  % Loop over all particles.
  for ii = 1:N
    % Loop over all particles with a greater index than the current index.
    for jj = ii+1:N
      
      % Find particle-particle separation.
      separation = norm(coor(ii,:)-coor(jj,:));
      
      % Build Distance Matrix.
      DistanceMatrix(ii,jj) = separation;
      DistanceMatrix(jj,ii) = separation;
      
      % Recalculate any points within the r0 of each other.
      if separation  < d0
        recalc(jj) = true;
      end
    end
  end
  
  % Get the number of points to recalculate.
  n_ = sum(recalc);
  
  n_list(1+mod(counter,listSize)) = n_;

  % Update counter.
  counter = counter + 1;
  
%   if mod(counter,listSize)==0 && n_ > 0 
%     newPerm_ = randperm(N);
%     recalc(newPerm_ (1:n_)) = true;
%     n_ = sum(recalc);
%   end
  % Get new coordinates for the selected particles.
  coor(recalc,:) = R*(2*rand(n_,3) -1);
  
  % If no particles needed to be recalculated then exit loop.
  if n_ ==0
    doCalc = false;
  end
end
counter
  disp(counter)
end




