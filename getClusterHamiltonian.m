% HC = sum_n( kron_i(G_ni) +sum_{m!~=n}(G_ni^m) ).
% G_ni = H_i for n==i, else G_ni = eye(size(H_i,1)).
% G_ni^m = H_i^m for i==n, H_i^n for i==m, else eye(size(H_i,1)).
function HC = getClusterHamiltonian(Hamiltonian_diagonal,Hamiltonian_offDiagonal, eState, Cluster, Nuclei,zeroIndex)

 
  dimension = prod(Nuclei.NumberStates(Cluster)); 
  HC = Hamiltonian_diagonal{eState,size(Hamiltonian_diagonal,2),size(Hamiltonian_diagonal,3)}*eye(dimension); % electron Zeeman.
  
  for  i_index_nucleus = Cluster
      inucleus = i_index_nucleus - zeroIndex;
    HC = HC + directProduct_diag(Hamiltonian_diagonal, eState, Cluster, Nuclei, inucleus, inucleus,zeroIndex);
    for  j_index_nucleus = Cluster
      jnucleus = j_index_nucleus - zeroIndex;
      if jnucleus >= inucleus
        break;
      end
       
      HC = HC + directProduct_diag(Hamiltonian_diagonal, eState, Cluster, Nuclei, jnucleus, inucleus,zeroIndex);
      H_temp = directProductOD(Hamiltonian_offDiagonal, eState, Cluster, Nuclei, jnucleus, inucleus,zeroIndex);
      HC = HC + H_temp + H_temp';
      
    end
  end
  if max( max( isnan(HC) ) ) || max( max( isinf(HC) ) )
    error('Cluster Hamiltonian could not be calculated.');
  end

  if max(max(abs(HC-HC'))) > 1e-9
    error('Hamiltonian is not Hermetian.');
  end
end

function Hdp = directProduct(Hamiltonian, eState, Cluster, Nuclei, inucleus, jnucleus, zeroIndex)

      Hdp = 1;
      
      for k_index_nucleus = Cluster
        knucleus = k_index_nucleus - zeroIndex;
    
          if knucleus == inucleus
              G = diag(Hamiltonian{eState,inucleus,jnucleus});
          elseif knucleus == jnucleus
              G = diag(Hamiltonian{eState,jnucleus,inucleus});
          else
              G = eye(Nuclei.NumberStates(k_index_nucleus));
          end
          
          if isempty(G)
              G = zeros(Nuclei.NumberStates(k_index_nucleus));
          end
          
          Hdp = kron(Hdp,G); 
        %  Hdp = kron_diag(diag(Hdp),diag(G));     
      end      
end
  
function Hdp = directProductOD(Hamiltonian, eState, Cluster, Nuclei, inucleus, jnucleus, zeroIndex)

      Hdp = 1;
      
      for  k_index_nucleus = Cluster
        knucleus = k_index_nucleus - zeroIndex;

          if knucleus == inucleus              
              if ((size(Hamiltonian,2)>=inucleus)&&(size(Hamiltonian,3)>=jnucleus))
                G = expandODHamiltonian(Hamiltonian{eState,inucleus,jnucleus});
              else 
                G = [];
              end
          elseif knucleus == jnucleus
              if ((size(Hamiltonian,2)>=jnucleus)&&(size(Hamiltonian,3)>=inucleus))
                G = expandODHamiltonian(Hamiltonian{eState,jnucleus,inucleus});
              else 
                G = [];
              end
          else
              G = eye(Nuclei.NumberStates(k_index_nucleus));
          end
          
          if isempty(G)
              G = zeros(Nuclei.NumberStates(k_index_nucleus));
          end
          
          Hdp = kron(Hdp,G);
               
      end      
end

function H = expandODHamiltonian(Hamiltonian)
  H = Hamiltonian(1:end-1);
  H = diag(H);
  if isempty(H)
    return;
  elseif Hamiltonian(end) > 0
    H = [zeros(size(H,1),1), H];
    H = [H ; zeros(1,size(H,2))];
  elseif Hamiltonian(end) < 0
    H = [H,zeros(size(H,1),1)];
    H = [zeros(1,size(H,2));H];     
  end
end

function Hdp = directProduct_diag(Hamiltonian, eState, Cluster, Nuclei, inucleus, jnucleus, zeroIndex)

      Hdp = 1;
      
      for k_index_nucleus = Cluster
        knucleus = k_index_nucleus - zeroIndex;
    
          if knucleus == inucleus
              G = Hamiltonian{eState,inucleus,jnucleus};
          elseif knucleus == jnucleus
              G = Hamiltonian{eState,jnucleus,inucleus};
          else
              G = ones(1,Nuclei.NumberStates(k_index_nucleus));
          end
          
          if isempty(G)
              G = zeros(1,Nuclei.NumberStates(k_index_nucleus));
          end
          
          %Hdp = kron(Hdp,G); 
          Hdp = kron_diag(Hdp,G);     
      end    
      Hdp = diag(Hdp);
end

function C = kron_diag(A,B)
  a = length(A);
  b = length(B);
  dim = a*b;
  C = zeros(dim,1);
  nB1 = 1;
  nB2 = 0;
  for n = 1:a 
    nB1= nB2+1;
    nB2 = n*b;
    C(nB1:nB2) = A(n)*B;
  end
  %C= diag(C);
end



