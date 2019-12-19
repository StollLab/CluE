% IiIjop = generateXiXjSpinOperators(S,maxClusterSize)
% S = spin.
% maxClusterSize = the maximum number of spins to be accounted for.
%
% IiIjop{clusterSize}}(:,:, x + 9n) is  the spin opererator that acts as IiIj
% for the nth spin in a cluster of size clusterSize, and as the identity for all
% other spins in the cluster; n in [0, maxClusterSize -1].  
%
% The number x gives the operator IiIj as shown below.
%
% ENUM for spin-operator indices.
% XX = 1;  XY = 4;  XZ = 7;
% YX = 2;  YY = 5;  YZ = 8;
% ZX = 3;  ZY = 6;  ZZ = 9;
function IiIjop = generateXiXjSpinOperators(S,maxClusterSize)

% Generate spin operators.
Iop{1} = spinX(S);
Iop{2} = spinY(S);
Iop{3} = spinZ(S);

% Form all 9 IiIj operators, in the order XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ
IiIj = cell(3,3);
for i = 1:3
  for j = 1:3
    IiIj{i,j} = Iop{i}*Iop{j};
  end
end

% Get the identity operator.
E = eye(2*S+1);

% Initialize output.
IiIjop = cell(maxClusterSize,1);

% Loop of cluster sizes.
for clusterSize = 1:maxClusterSize
  
  % Determine cluster Hilbert space dimensionality.
  N = (2*S+1)^clusterSize;
  
  % Initialize the set of cluster operators of order clusterSize. 
  IiIjop{clusterSize} = zeros(N,N,clusterSize*9);
  
  % Define an index.
  idx = 1;
  
  % Loop over nuclei within a the cluster.
  for iNuc = 1:clusterSize
    
    % Loop over the 9 IiIj operators.
    for iOp = 1:9
      
      % Initilize cluster spin-operator.
      Iop = 1;
      
      % Loop over cluster positions.
      for j = 1:clusterSize
        
        if j==iNuc % spin position == cluster position
          % Use IiIj operator.
          op_ = IiIj{iOp};
        else
          % Use identity.
          op_ = E;
        end
        
        % Build up cluster spin-operator with the appropriate spin operator. 
        Iop = kron(Iop,op_);
      end
      
      % Add cluster spin-operator to the output.
      IiIjop{clusterSize}(:,:,idx) = Iop;
      
      % Increment index.
      idx = idx + 1;
    end
  end
  
end
