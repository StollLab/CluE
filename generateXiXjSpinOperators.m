function IiIjop = generateXiXjSpinOperators(S,maxClusterSize)

Iop{1} = spinX(S);
Iop{2} = spinY(S);
Iop{3} = spinZ(S);

for i = 1:3
  for j = 1:3
    IiIj{i,j} = Iop{i}*Iop{j};
  end
end

E = eye(2*S+1);
IiIjop = cell(maxClusterSize,1);

for clusterSize = 1:maxClusterSize
  N = (2*S+1)^clusterSize;
  IiIjop{clusterSize} = zeros(N,N,clusterSize*9);
  idx = 1;
  for iNuc = 1:clusterSize
    for iOp = 1:9
      
      Iop = 1;
      for j = 1:clusterSize
        if j==iNuc
          op_ = IiIj{iOp};
        else
          op_ = E;
        end
        Iop = kron(Iop,op_);
      end
      IiIjop{clusterSize}(:,:,idx) = Iop;
      idx = idx + 1;
    end
  end
  
end
