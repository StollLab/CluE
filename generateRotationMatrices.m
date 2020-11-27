function rot = generateRotationMatrices(spinDim,numberMethyl)

% identity matrix for spin dof
identity = eye(spinDim);

% rotational transition matrix
K = [0,1,1; 1,0,1; 1,1,0];

% 3 by 3 identity matrix
E3 = eye(3);

% Initialize output.
rot = zeros(spinDim*3^numberMethyl);

% Loop over methyl rotors.
for ii= 1:numberMethyl

  % direct product over all methyl transition matrices.
  K_ = 1;
  for jj = 1:numberMethyl
    if ii==jj
      K_ = kron(K,K_);
    else
      K_ = kron(E3,K_);
    end
  end  
  rot =  rot + kron(K_,identity);
  
end

end