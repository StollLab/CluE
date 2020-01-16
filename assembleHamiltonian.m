function [H_alpha,H_beta] = assembleHamiltonian(tensors,Cluster,System,Nuclei,zeroIndex,clusterSize)

theory = System.theory;

useEZ       = theory(1);
useNZ       = theory(2);
useHF_SzIz  = theory(3);
useHF_SzIxy = theory(4);
useNucA     = theory(5);
useNucB     = theory(6);
useNucCD    = theory(7);
useNucEF    = theory(8);
useNQ       = theory(9);

%{
Cluster0 = Cluster;
for ispin = Cluster0
  if strcmp('CH3', Nuclei.Type{ispin})
    clusterSize = clusterSize +2;
    %       Cluster = [Cluster, Nuclei.Auxiliary_ID(ispin,:)];
    Cluster = [Cluster, [ispin+1,ispin+2,ispin+3]];
    Cluster(Cluster==ispin)=[];
  end
end
%}
Cluster = sort(unique(Cluster));

if abs(double(zeroIndex) + 1 -double(Cluster(1)))>=1
  error('Cluster reference failure.');
end

% dimension of the electron manifold
dim = prod(Nuclei.NumberStates(Cluster));

% initialize output cluster Hamiltonian
H_out = eye(dim);

Hnuc = 0;
Hhf = 0;

% icluster is a index for the cluster
for iSpin = 1:clusterSize
  
  % nuclear spin quantum number
  I = Nuclei.Spin(Cluster(iSpin));
  
  % the ith nuclear spin in the input Hamiltonian
  inucleus = Cluster(iSpin) - zeroIndex + 1; % since the electron is given position 1;
  
  % Hilbert subspace dimensions nuclei 1:i-1, i+1:end
  pre_dim = prod(Nuclei.NumberStates(Cluster(1:iSpin-1)));
  post_dim = prod(Nuclei.NumberStates(Cluster(iSpin+1:end)));

  if useNZ || HF
    Iz = fullop(pre_dim,spinZ(I),post_dim);
  end
  
  % Calculate nuclear Zeeman Hamiltonian
  if useNZ
    H_nuclear_Zeeman_Iz = -tensors{inucleus,inucleus}(3,3)*Iz;
  end
  
  % Calculate hyperfine Hamiltonian
  A = zeros(3);
  if ~isempty(tensors{1,inucleus})
    A = A + tensors{1,inucleus};
  end
  if ~isempty(tensors{inucleus,1})
    A = A + tensors{inucleus,1};
  end
  if useHF_SzIz
    H_hyperfine_SzIz = -A(3,3)*Iz;
  else
    H_hyperfine_SzIz = 0;
  end
  if useHF_SzIxy
    Ix = fullop(pre_dim,spinX(I),post_dim);
    Iy = fullop(pre_dim,spinY(I),post_dim);
    H_hyperfine_SzIx = -A(3,1)*Ix;
    H_hyperfine_SzIy = -A(3,2)*Iy;
  else
    H_hyperfine_SzIx = 0;
    H_hyperfine_SzIy = 0;
  end
    
  % Assemble single-nucleus terms in nuclear Hamiltonian
  Hnuc = Hnuc + H_nuclear_Zeeman_Iz + H_nuclear_quadrupole;
  Hhf = Hhf + H_hyperfine_SzIz + H_hyperfine_SzIx + H_hyperfine_SzIy;
  
  H_out = H_out + hyperfine_;
  
  % loop over all nuclei with index greater than the ith nucleus
  for jcluster = iSpin+1:clusterSize
    
    % the jth nuclear spin in the system
    jspin = Cluster(jcluster);
    
    % the jth nuclear spin in the input Hamiltonian
    jnucleus = jspin - zeroIndex + 1; % since the electron is given position 1;
    
    % dimension of the identity matrix between the ith and jth nucleus
    mid_dim = prod(Nuclei.NumberStates(Cluster(iSpin+1:jcluster-1)));
    
    % dimension of the identity matrix following the jth nucleus
    post_dim = prod(Nuclei.NumberStates(Cluster(jcluster+1:end)));
    
    % dimension of the matrix of the ith nucleus
    dim_j =  prod(Nuclei.NumberStates(jspin));
    
    % get dipolar coupling tensor
    dd = zeros(3);
    if ~isempty(tensors{inucleus,jnucleus})
      dd = dd +   tensors{inucleus,jnucleus};
    end
    if ~isempty(tensors{jnucleus,inucleus})
      dd = dd +   tensors{jnucleus,inucleus};
    end
    
    % nucleus-nucleus secular (A term)
    if useNucA
      IzJz = fullop2(pre_dim,spinZ(I),mid_dim,spinZ(Nuclei.Spin(jspin)),post_dim0);
      Hnn_A = dd(3,3)*IzJz;
    else
      Hnn_A = 0;
    end
    
    % nucleus-nucleus pseudo-secular (B term)      
    if useNucB
      IpJm = fullop2(pre_dim,spinRaise(I),mid_dim,spinLower(Nuclei.Spin(jspin)),post_dim);
      ImJp = IpJm';
      Hnn_B = 0.25*(dd(1,1) + dd(2,2))*(IpJm+ImJp);
    else
      Hnn_B = 0;  
    end
    
    % nucleus-nucleus dipolar C and D terms
    if useNucCD
      IzJp = fullop2(pre_dim,spinZ(I),mid_dim,spinRaise(Nuclei.Spin(jspin)),post_dim);      
      IpJz = fullop2(pre_dim,spinRaise(I),mid_dim,spinZ(Nuclei.Spin(jspin)),post_dim);      
      IzJm = IzJp';
      ImJz = IpJz';
      cd = 1/2*(dd(1,3) - 1i*dd(2,3));
      Hnn_CD = cd*(IzJp+IpJz) + cd'*(IzJm+ImJz);
    else
      Hnn_CD = 0;      
    end
    
    % nucleus-nucleus dipolar E and F terms
    if useNucEF
      IpJp = fullop2(pre_dim,spinRaise(I),mid_dim,spinRaise(Nuclei.Spin(jspin)),post_dim);
      ImJm = IpJp';
      ef = 1/4*(dd(1,1) - dd(2,2) - 1i*dd(1,2) - 1i*dd(2,1));
      Hnn_EF = ef*IpJp + ef'*ImJm;
    else
      Hnn_EF = 0;
    end
    
    Hnuc = Hnuc + Hnn_A + Hnn_B + Hnn_CD + Hnn_EF;
    
  end
end

% Calculate electron Zeeman Hamiltonian
if useEZ
  eZeeman = tensors{1,1}*eye(3); % Hz
  HEZ = eZeeman(3,3)*eye(size(Hnuc)); % Hz
else
  HEZ = 0;
end

% Calculate total nuclear Hamiltonians for alpha and beta electron manifolds
H_alpha = +1/2*(HEZ + Hhf) + Hnuc;
H_beta  = -1/2*(HEZ + Hhf) + Hnuc;

% Check Hermitianity
threshold = 1e-12;
[isHermA,nonHermiticityA] = isHermitian(H_alpha,threshold);
[isHermB,nonHermiticityB] = isHermitian(H_beta,threshold);
if ~isHermA || ~isHermB
  hline = '--------------------------------------------------------------';
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
  fprintf('Non-Hermiticity = {%d,%d}.\n',nonHermiticityA,nonHermiticityB);
  disp(hline);
  
  [isHerm,nonHermiticity] = isHermitian(HEZ,threshold);
  fprintf('HEZ non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('pass = %d.\n',isHerm); disp(hline);
  
  [isHerm,nonHermiticity] = isHermitian(Hhf,threshold);
  fprintf('Hhf non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('pass = %d.\n',isHerm); disp(hline);
  
  [isHerm,nonHermiticity] = isHermitian(H_hyperfine_SzIz,threshold);
  fprintf('  H_hyperfine_SzIz non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);

  [isHerm,nonHermiticity] = isHermitian(H_hyperfine_SzIy,threshold);
  fprintf('  H_hyperfine_SzIy non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);
  
  [isHerm,nonHermiticity] = isHermitian(H_hyperfine_SzIx,threshold);
  fprintf('  H_hyperfine_SzIx non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);

  [isHerm,nonHermiticity] = isHermitian(Hnuc,threshold);
  fprintf('Hnuc non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('pass = %d.\n',isHerm); disp(hline);
  
  [isHerm,nonHermiticity] = isHermitian(Hnn_A,threshold);
  fprintf('  Hnn_A non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);
  
  [isHerm,nonHermiticity] = isHermitian(Hnn_B,threshold);
  fprintf('  Hnn_B non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);
  
  [isHerm,nonHermiticity] = isHermitian(Hnn_CD,threshold);
  fprintf('  Hnn_CD non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);
  
  [isHerm,nonHermiticity] = isHermitian(Hnn_EF,threshold);
  fprintf('  Hnn_EF non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);
    
  [isHerm,nonHermiticity] = isHermitian(H_nuclear_quadrupole,threshold);
  fprintf('  H_nuclear_quadrupole non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);
  
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
  error('Cluster Hamiltonian is not Hermitian.');
end

% Hermitianize
H_alpha = (H_alpha+H_alpha')/2;
H_beta = (H_beta+H_beta')/2;

%{
% Find all methyls groups.
selectionRule = eye(size(H_out));
for ispin = Cluster0
  if strcmp('CH3', Nuclei.Type{ispin})
    Methyl_Cluster = [ispin + 1,ispin + 2,ispin + 3];
    pre_dim = prod(Nuclei.NumberStates(Cluster(Cluster<ispin)));
    post_dim = prod(Nuclei.NumberStates(Cluster(Cluster>ispin+3)));
    H_methyl = assembleHamiltonian(tensors,Methyl_Cluster,System, eState,Nuclei,ispin,3);
    H_methyl = kron(eye(pre_dim),kron(H_methyl,eye(post_dim)));
    H_out = H_out -H_methyl + H_methyl.*selectionRule;
  end
end
%}

end

function Op = fullop(n_pre,op,n_post)
Op = eye(n_pre);
Op = kron(Op,op);
Op = kron(Op,eye(n_post));
end

function Op = fullop2(n_pre,op1,n_mid,op2,n_post)
Op = eye(n_pre);
Op = kron(Op,op1);
Op = kron(Op,eye(n_mid));
Op = kron(Op,op2);
Op = kron(Op,eye(n_post));
end
