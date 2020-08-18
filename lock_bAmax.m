function bAmax_lock = lock_bAmax(pwstat,Method)

notDiag = eye(size(pwstat.Nuclear_Dipole_perpendicular))==0;
% Find maximum |b|.
b = abs( pwstat.Nuclear_Dipole_perpendicular(notDiag) );
bmax = max( max( b ) );

% Find all edges within 10 percent of bmax.
selection_matrix = abs(pwstat.Nuclear_Dipole_perpendicular(notDiag)) >= 0.9*bmax; 

% Find the sqrt|bAmax| value for these edges.
bAmax = abs(pwstat.bAmax(notDiag));
bAmax_set = bAmax(selection_matrix);

% Find the maximum and minimum sqrt|bAmax| value 
% from the top 10 percent of b.
bAmax_min = min(bAmax_set);
bAmax_max = max(bAmax_set);

% Find the geometric mean of the max and min sqrt|bAmax|.
bAmax_1 = (bAmax_min^(1-Method.lockValue) )*bAmax_max^(Method.lockValue);


bmin = Method.cutoff.dipole(1);

bAmax_0 = sqrt(bmin/bmax)*bAmax_1; 

bAmax_lock = min(bAmax_0,bAmax_1);

if Method.lock_bAmax
  fprintf('bAmax_lock = %d Hz.\n',bAmax_lock);
end

end