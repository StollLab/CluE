function [spinOp,IiIjop_HD,IiIjop_DH] = assembleMixedSpinOperator(SpinOp_1,SpinOp_2,doXiXjOps)

% ENUM
E = 1;  Z = 2; RAISE = 3; LOWER = 4;

% Initialize output.
nOp = 2;
dim = size(SpinOp_1,1)*size(SpinOp_2,1);
spinOp = zeros(dim,dim,4^nOp);
opNames = getOpNames(nOp);

for ii =1:4^nOp
  
  opstr = opNames{ii};
  
  % Identify single-spin operator;
  switch opstr(1)
    case 'E'
      op1_ = E;
    case 'Z'
      op1_ = Z;
    case 'R'
      op1_ = RAISE;
    case 'L'
      op1_ = LOWER;
    otherwise
      disp(['Unidentified reference ', opstr(1), '.']);
      error('Could not identify spin operator.');
  end
  
  switch opstr(2)
    case 'E'
      op2_ = E;
    case 'Z'
      op2_ = Z;
    case 'R'
      op2_ = RAISE;
    case 'L'
      op2_ = LOWER;
    otherwise
      disp(['Unidentified reference ', opstr(2), '.']);
      error('Could not identify spin operator.');
  end
  
  % Krocker the operator into the output.
  spinOp(:,:,ii) = kron(SpinOp_1(:,:,op1_),SpinOp_2(:,:,op2_));
end


E2 = eye(2);
S = 1;
Iop{1} = spinX(S);
Iop{2} = spinY(S);
Iop{3} = spinZ(S);

if ~doXiXjOps
  IiIjop_HD = [];
  IiIjop_DH = [];
else
  IiIjop_HD = zeros(6,6,9);
  IiIjop_DH = zeros(6,6,9);
  % Loop over the 9 IiIj operators.
  kk = 0;
  for ii = 1:3
    for jj = 1:3
      kk = kk + 1;
      IiIj = Iop{ii}*Iop{jj};
      IiIjop_HD(:,:,kk) = kron(E2,IiIj);
      IiIjop_DH(:,:,kk) = kron(IiIj,E2);
    end
  end
  
end
end