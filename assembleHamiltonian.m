function H_out = assembleHamiltonian(H_in,Cluster,System, eState,Nuclei,zeroIndex,clusterSize)


  Cluster0 = Cluster;
  for ispin = Cluster0
    if strcmp('CH3', Nuclei.Type{ispin})
      clusterSize = clusterSize +2;
      %       Cluster = [Cluster, Nuclei.Auxiliary_ID(ispin,:)];
      Cluster = [Cluster, [ispin+1,ispin+2,ispin+3]];
      Cluster(Cluster==ispin)=[];
    end
  end

  Cluster = sort(unique(Cluster));

if abs(double(zeroIndex) + 1 -double(Cluster(1)))>=1
  error('Cluster reference failure.');
end

n_cluster = length(Cluster);

if clusterSize ~= n_cluster
  error('Cluster reference failure.');
end

% electron manifold
ms = eState - 1  - System.Electron.spin;

% dimension of the electron manifold
dim = prod(Nuclei.NumberStates(Cluster));

% initialize output cluster Hamiltonian
H_out = eye(dim);

% electron Zeeman splitting
eZeeman = H_in{1,1}*eye(3); % Hz

% electron Zeeman frequency
H_out = ms*H_out*eZeeman(3,3); % Hz


% icluster is a index for the cluster
for icluster = 1:n_cluster
  
  % the ith nuclear spin in the system
  ispin = Cluster(icluster) ;
  
  % the ith nuclear spin in the input Hamiltonian
  inucleus = ispin - zeroIndex + 1; % since the electron is given position 1;
  
  % dimension of the identity matrix preceding the ith nucleus
  pre_dim = prod(Nuclei.NumberStates(Cluster(1:icluster-1)));
  
  % dimension of the identity matrix following the ith nucleus
  post_dim = prod(Nuclei.NumberStates(Cluster(icluster+1:end)));
  
  % dimension of the matrix of the ith nucleus
  dim_i =  prod(Nuclei.NumberStates(ispin));
  
  
  % nuclear Zeeman
  
  K = zeros(3);
  if isempty(H_in{inucleus,inucleus})
    continue;
  end
  K = K + H_in{inucleus,inucleus};
  
  % spin matrix
  %clear matrix;
  matrix = kron(eye(pre_dim), spinZ(Nuclei.Spin(ispin))) ;
  matrix = kron(matrix, eye(post_dim));
  
  if abs(trace(matrix))>1e-12
    error('Spin Matrix is not traceless.')
  end
  
  nZeeman_ = -K(3,3)*matrix; % e kron e ... kron Iz_inucleus kron ... kron e.
  H_out = H_out + nZeeman_;
  
  % secular hyperfine
  
  % K = H_in{1,inucleus}*eye(3) + H_in{inucleus,1}*eye(3);
  % but either H_in may be empty
  K = zeros(3);
  if ~isempty(H_in{1,inucleus})
    K = K +   H_in{1,inucleus};
  end
  if ~isempty(H_in{inucleus,1})
    K = K +   H_in{inucleus,1};
  end
  
  % hyperfine spin Hamiltonian
  hyperfine_ =  -K(3,3)*ms*matrix; % e kron e ... kron Iz_inucleus kron ... kron e.
  
  if System.full_Sz_Hyperfine
    % spin matrix
    %clear matrix;
    matrix = kron(eye(pre_dim), spinX(Nuclei.Spin(ispin))) ;
    matrix = kron(matrix, eye(post_dim));
    hyperfine_SzIx = -K(1,3)*ms*matrix;
    
    matrix = kron(eye(pre_dim), spinY(Nuclei.Spin(ispin))) ;
    matrix = kron(matrix, eye(post_dim));
    hyperfine_SzIy = -K(2,3)*ms*matrix;
    
    hyperfine_ =  hyperfine_  + hyperfine_SzIx + hyperfine_SzIy;
    
    %   hyperfine_SzIxy = +(K(1,3) + 1i*K(2,3))*ms*matrix/4;
    %    hyperfine_SzIxy =  hyperfine_SzIxy +  hyperfine_SzIxy';
    %   hyperfine_ =  hyperfine_  + hyperfine_SzIxy;
  end
  
  H_out = H_out + hyperfine_;
  
  
  % loop over all nuclei with index greater than the ith nucleus
  for jcluster = icluster+1:n_cluster
    if jcluster <= icluster
      continue;
    end
    
    % the jth nuclear spin in the system
    jspin = Cluster(jcluster);
    
    % the jth nuclear spin in the input Hamiltonian
    jnucleus = jspin - zeroIndex + 1; % since the electron is given position 1;
    
    % dimension of the identity matrix between the ith and jth nucleus
    mid_dim = prod(Nuclei.NumberStates(Cluster(icluster+1:jcluster-1)));
    
    % dimension of the identity matrix following the jth nucleus
    post_dim = prod(Nuclei.NumberStates(Cluster(jcluster+1:end)));
    
    % dimension of the matrix of the ith nucleus
    dim_j =  prod(Nuclei.NumberStates(jspin));
    
    % get H_ij
    K = zeros(3);
    if ~isempty(H_in{inucleus,jnucleus})
      K = K +   H_in{inucleus,jnucleus};
    end
    if ~isempty(H_in{jnucleus,inucleus})
      K = K +   H_in{jnucleus,inucleus};
    end
    
    % use full dipole-dipole Hamiltonian
    if System.fullDipoleTensor
      % i spin vector
      Spin_i = {spinX(Nuclei.Spin(ispin)),spinY(Nuclei.Spin(ispin)),spinZ(Nuclei.Spin(ispin))};
      
      % j spin vector
      Spin_j = {spinX(Nuclei.Spin(jspin)),spinY(Nuclei.Spin(jspin)),spinZ(Nuclei.Spin(jspin))};
      
      % spin dipolar Hamilltonian
      Hdd = zeros(dim);
      
      for ix = 1:3
        for jx = 1:3
          Gdd_left_ = kron(eye(pre_dim), eye(dim_i));
          Gdd_left_ = kron( Gdd_left_, eye(mid_dim) );
          Gdd_left_ = kron( Gdd_left_,Spin_j{jx}  );
          Gdd_left_ = kron(Gdd_left_,eye(post_dim) );
          
          Gdd_right_ = kron(eye(pre_dim), Spin_i{ix});
          Gdd_right_ = kron( Gdd_right_ , eye(mid_dim) );
          Gdd_right_ = kron( Gdd_right_ ,eye(dim_j)  );
          Gdd_right_ = kron(Gdd_right_ ,eye(post_dim) );
          
          %           Hdd = Hdd - kron( eye(dim_i),Spin_j{jx}  )'*K(jx,ix) * kron( Spin_i{ix},eye(dim_j) );
          Hdd = Hdd + Gdd_left_'*K(jx,ix) * Gdd_right_;
          
        end
      end
      H_out = H_out + Hdd;
      continue;
    end
    
    % or just the secular and pseudo-secular
    if System.nuclear_dipole_A
      % nucleus-nucleus secular
      
      %clear matrix;
      matrix = kron(eye(pre_dim), spinZ(Nuclei.Spin(ispin)) );
      matrix = kron(matrix, eye(mid_dim));
      matrix = kron(matrix, spinZ(Nuclei.Spin(jspin)) );
      matrix = kron(matrix, eye(post_dim));
      
      if abs(trace(matrix))>1e-12
        error('Spin Matrix is not traceless.')
      end
      
      % select Iz*Iz
      H_out = H_out + K(3,3)*matrix;
    end
    
    if System.nuclear_dipole_B
      % nucleus-nucleus pseudo-secular
      
      %clear matrix;
      matrix = kron(eye(pre_dim), spinRaise(Nuclei.Spin(ispin)) );
      matrix = kron(matrix, eye(mid_dim));
      matrix = kron(matrix, spinLower(Nuclei.Spin(jspin)) );
      matrix = kron(matrix, eye(post_dim));
      matrix = matrix + matrix';
      
      if abs(trace(matrix))>1e-12
        error('Spin Matrix is not traceless.')
      end
      
      
      bath = 0.25*(K(1,1) + K(2,2)); % e kron e ... kron Iz_inucleus kron ... kron e ... kron Iz_jnucleus ... kron e.
      
      H_out = H_out + bath*matrix;
      
    end
    
    if System.nuclear_dipole_CD
      %clear matrix;
      matrix = kron(eye(pre_dim), spinRaise(Nuclei.Spin(ispin)) );
      matrix = kron(matrix, eye(mid_dim));
      matrix = kron(matrix, spinZ(Nuclei.Spin(jspin)) );
      matrix = kron(matrix, eye(post_dim));      
      matrix_ = kron(eye(pre_dim), spinZ(Nuclei.Spin(ispin)) );
      matrix_ = kron(matrix_, eye(mid_dim));
      matrix_ = kron(matrix_, spinRaise(Nuclei.Spin(jspin)) );
      matrix_ = kron(matrix_, eye(post_dim));      
      
      matrix = 1/2*(K(1,3) - 1i*K(2,3))*(matrix + matrix_);
      H_out = H_out + matrix + matrix';
      
    end

    if System.nuclear_dipole_EF
      %clear matrix;
      matrix = kron(eye(pre_dim), spinRaise(Nuclei.Spin(ispin)) );
      matrix = kron(matrix, eye(mid_dim));
      matrix = kron(matrix, spinRaise(Nuclei.Spin(jspin)) );
      matrix = kron(matrix, eye(post_dim));
      
      matrix  = 1/4*(K(1,1)  - K(2,2) - 1i*K(1,2) - 1i*K(2,1))*matrix;
      H_out = H_out + matrix + matrix';
    end 
    
  end
end
if max(max(abs( H_out - H_out') ) )>1e-12
  error('Cluster Hamiltonian is not Hermitian.');
end


% Find all methyls groups. 
selectionRule = eye(size(H_out));
 for ispin = Cluster0
    if strcmp('CH3', Nuclei.Type{ispin})
      Methyl_Cluster = [ispin + 1,ispin + 2,ispin + 3];
      pre_dim = prod(Nuclei.NumberStates(Cluster(Cluster<ispin)));
      post_dim = prod(Nuclei.NumberStates(Cluster(Cluster>ispin+3)));
      H_methyl = assembleHamiltonian(H_in,Methyl_Cluster,System, eState,Nuclei,ispin,3);
      H_methyl = kron(eye(pre_dim),kron(H_methyl,eye(post_dim)));
      H_out = H_out -H_methyl + H_methyl.*selectionRule;
    end
 end

if max(max(abs( H_out - H_out') ) )>1e-12
  error('Cluster Hamiltonian is not Hermitian.');
end

end


