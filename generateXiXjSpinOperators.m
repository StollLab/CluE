function [SpinOp_1,SpinOp_2,SpinOp_3,SpinOp_4,SpinOp_5,SpinOp_6]= generateXiXjSpinOperators(spin)
% enum
X = 1;  Y = 2; Z = 3;
XX = 4; XY = 5;  XZ = 6;
YX = 7; YY = 8;  YZ = 9;
ZX = 10; ZY = 11;  ZZ = 12;

numOps = 12;
OO = [XX, XY, XZ; YX,YY,YZ; ZX,ZY,ZZ];
XYZ = [X,Y,Z];

% Set spin and spin multiplicity;
multiplicity = 2*spin+1; 

% 1-Cluster----------------------------------------------------------------
SpinOp_1 = zeros(multiplicity,multiplicity,numOps);


SpinOp_1(:,:,X) = spinX(spin);
SpinOp_1(:,:,Y) = spinY(spin);
SpinOp_1(:,:,Z) = spinZ(spin);

for ii = XYZ
  for jj = XYZ
    SpinOp_1(:,:,OO(ii,jj)) = SpinOp_1(:,:,ii)*SpinOp_1(:,:,jj);
  end
end

% 2-Cluster---------------------------------------------------------------
n = 2;
SpinOp_2 = zeros(multiplicity^n,multiplicity^n,2*numOps);
E1 = eye(multiplicity);
for ii = 1:numOps
SpinOp_2(:,:,ii) = kron(E1,SpinOp_1(:,:,ii));
SpinOp_2(:,:,ii+numOps) = kron(SpinOp_1(:,:,ii),E1);
end

% 3-Cluster----------------------------------------------------------------
n = 3;
SpinOp_3 = zeros(multiplicity^n,multiplicity^n,3*numOps);
E2 = eye(multiplicity^2);

for ii = 1:numOps
SpinOp_3(:,:,ii) = kron(E2,SpinOp_1(:,:,ii));
SpinOp_3(:,:,ii+numOps) = kron(SpinOp_2(:,:,ii),E1);
SpinOp_3(:,:,ii+2*numOps) = kron(SpinOp_1(:,:,ii),E2);
end

% 4-Cluster----------------------------------------------------------------
n = 4;
SpinOp_4 = zeros(multiplicity^n,multiplicity^n,4*numOps);
E3 = eye(multiplicity^3);

for ii = 1:numOps
SpinOp_4(:,:,ii)            = kron(              E3 , SpinOp_1(:,:,ii) );
SpinOp_4(:,:,ii + numOps)   = kron(SpinOp_3(:,:,ii) , E1               );
SpinOp_4(:,:,ii + 2*numOps) = kron(SpinOp_2(:,:,ii) , E2               );
SpinOp_4(:,:,ii + 3*numOps) = kron(SpinOp_1(:,:,ii) , E3               );
end

% 5-Cluster----------------------------------------------------------------
n = 5;
SpinOp_5 = zeros(multiplicity^n,multiplicity^n,5*numOps);
E4 = eye(multiplicity^4);

for ii = 1:numOps
SpinOp_5(:,:,ii)            = kron(              E4 , SpinOp_1(:,:,ii) );
SpinOp_5(:,:,ii + numOps)   = kron(SpinOp_4(:,:,ii) , E1               );
SpinOp_5(:,:,ii + 2*numOps) = kron(SpinOp_3(:,:,ii) , E2               );
SpinOp_5(:,:,ii + 3*numOps) = kron(SpinOp_2(:,:,ii) , E3               );
SpinOp_5(:,:,ii + 4*numOps) = kron(SpinOp_1(:,:,ii) , E4               );
end


% 6-Cluster----------------------------------------------------------------
n = 6;
SpinOp_6 = zeros(multiplicity^n,multiplicity^n,6*numOps);
E5 = eye(multiplicity^5);

for ii = 1:numOps
SpinOp_6(:,:,ii)            = kron(              E5 , SpinOp_1(:,:,ii) );
SpinOp_6(:,:,ii + numOps)   = kron(SpinOp_5(:,:,ii) , E1               );
SpinOp_6(:,:,ii + 2*numOps) = kron(SpinOp_4(:,:,ii) , E2               );
SpinOp_6(:,:,ii + 3*numOps) = kron(SpinOp_3(:,:,ii) , E3               );
SpinOp_6(:,:,ii + 4*numOps) = kron(SpinOp_2(:,:,ii) , E4               );
SpinOp_6(:,:,ii + 5*numOps) = kron(SpinOp_1(:,:,ii) , E5               );
end

