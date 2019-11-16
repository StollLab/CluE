function [P3,P4,P5,P6] = getMethylProjections(P,E)

numOp = size(P,3);
dimE = size(E,1);
dimP = size(P,1);

P3 = P;

P4 = zeros(dimP*dimE,dimP*dimE,2*numOp);
for iop = 1:numOp
  P4(:,:,iop) = kron(E,P(:,:,iop));
  P4(:,:,iop+numOp) = kron(P(:,:,iop),E);
end

P5 = zeros(dimP*dimE^2,dimP*dimE^2,3*numOp);
for iop = 1:numOp
  P5(:,:,iop) = kron(E,kron(E,P(:,:,iop)));
  P5(:,:,iop+numOp) = kron(E,kron(P(:,:,iop),E));
  P5(:,:,iop+2*numOp) = kron(P(:,:,iop),kron(E,E));
end

P6 = zeros(dimP*dimE^3,dimP*dimE^3,6*numOp);
for iop = 1:numOp
  P6(:,:,iop) = kron(E,kron(E,kron(E,P(:,:,iop))));
  P6(:,:,iop+numOp) = kron(E,kron(E,kron(P(:,:,iop),E)));
  P6(:,:,iop+2*numOp) = kron(E,kron(P(:,:,iop),kron(E,E)));
  P6(:,:,iop+3*numOp) = kron(kron(kron( P(:,:,iop),E),E),E);
  P6(:,:,iop + 4*numOp) = kron(P(:,:,iop),eye(dimP));
  P6(:,:,iop + 5*numOp) = kron(eye(dimP),P(:,:,iop));
end

end













