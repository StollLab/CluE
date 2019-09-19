function H_out = importHamiltonian(H_in,Cluster,System, eState,Nuclei)
  ms = eState - 1 - System.Electron.spin;
  dim = prod(Nuclei.NumberStates(Cluster));

  H_out = eye(dim);
  % electron Zeeman
  eZeeman = H_in{1,1}*eye(3); % e kron e kron ... kron e.
  H_out = H_out*eZeeman(3,3);
  
  for ii = Cluster
    inucleus = ii +1; % since the electron is given position 1;
    dim_ii = Nuclei.NumberStates(ii);
    pre_dim = prod(Nuclei.NumberStates(Cluster(1:ii-1)));
    post_dim = prod(Nuclei.NumberStates(Cluster(ii+1:end)));
    matrix = kron(eye(pre_dim), spinZ(Nuclei.Spin(ii))) ;
    matrix = kron(matrix, eye(post_dim));
    
    % nuclear Zeeman
    K = H_in{inucleus,inucleus}*eye(3);
    nZeeman_ = K(3,3)*matrix; % e kron e ... kron Iz_inucleus kron ... kron e.
    H_out = H_out + nZeeman_;

    % secular hyperfine
    K = H_in{1,inucleus}*eye(3);% + H_in{inucleus,1}*eye(3);
    hyperfine_ = K(3,3)*ms*matrix;% e kron e ... kron Iz_inucleus kron ... kron e.
    H_out = H_out + hyperfine_;
    
    for jj = Cluster
      if jj <= ii
        continue;
      end
      
      mid_dim = prod(Nuclei.NumberStates(Cluster(ii+1:jj-1)));
      post_dim = prod(Nuclei.NumberStates(Cluster(jj+1:end)));
      dim_jj = Nuclei.NumberStates(jj);
      % nucleus-nucleus secular
      
      matrix = kron(eye(pre_dim), spinZ(Nuclei.Spin(ii)) );
      matrix = kron(matrix, eye(mid_dim));
      matrix = kron(matrix, spinZ(Nuclei.Spin(jj)) );
      matrix = kron(matrix, eye(post_dim));
      
      jnucleus = jj + 1;
      K = H_in{inucleus,jnucleus}*eye(3);% + H_in{jnucleus,inucleus}*eye(3);
      H_out = H_out + K(3,3)*matrix; 
      
      % nucleus-nucleus pseudo-secular
      
      bath = 0.25*(K(1,1) + K(2,2)); % e kron e ... kron Iz_inucleus kron ... kron e ... kron Iz_jnucleus ... kron e.
      
      matrix = kron(eye(pre_dim), spinRaise(Nuclei.Spin(ii)) );
      matrix = kron(matrix, eye(mid_dim));
      matrix = kron(matrix, spinLower(Nuclei.Spin(jj)) );
      matrix = kron(matrix, eye(post_dim));
      matrix = matrix + matrix';
      H_out = H_out + bath*matrix;
    end
  end

end


function Iz = spinZ(spin)
Iz = eye(2*spin+1);
for ii = 1:(2*spin+1)
  Iz(ii,ii) = spin+ 1 - ii;
end
end

function Iminus = spinLower(spin)
Iminus = zeros(2*spin+1);
for ii = 1:(2*spin+1)-1
  Iminus(ii+1,ii) = sqrt(spin*(spin+1)-(spin - ii)*(spin+ 1 - ii));
end
end

function Iplus = spinRaise(spin)
Iplus = zeros(2*spin+1);
for ii = 1:(2*spin+1)-1
  Iplus(ii,ii+1) = sqrt(spin*(spin+1)-(spin - ii)*(spin+ 1 - ii));
end
end