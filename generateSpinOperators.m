function output = generateSpinOperators(spin,maxClusterSize)

% Initialize output.
output = cell(1,maxClusterSize);

% Set spin and spin multiplicity.
multiplicity = 2*spin+1;

% ENUM
E = 1;  Z = 2; RAISE = 3; LOWER = 4;

% Initialize 1-cluster operators.
SpinOp_1 = zeros(multiplicity,multiplicity,4);

% Set 1-cluster operators.
SpinOp_1(:,:,E) = eye(multiplicity);
SpinOp_1(:,:,Z) = spinZ(spin);
SpinOp_1(:,:,RAISE) = spinRaise(spin);
SpinOp_1(:,:,LOWER) = spinLower(spin);

% Set 1-cluster operators to output.
output{1} = SpinOp_1;
 
% Loop throurh cluster sizes > 1.
for clusterSize = 2:maxClusterSize
  
  % Get set of operator names.
  OpNames = getOpNames(clusterSize);
  
  % Get number of operators for the given cluster size.
  numOps_ = 1 + 3*clusterSize + 9*NchooseK(clusterSize,2);
  
  % Initialize output for the given cluster size.
  output{clusterSize} = ...
    zeros(multiplicity^clusterSize,multiplicity^clusterSize,numOps_);
  
  % Loop over operator indices. 
  for iop =1:numOps_
    
    % Get the name of the indexed operator.
    opstr = OpNames{iop};
    
    % Construct operator.
    output{clusterSize}(:,:,iop) = assembleSpinOperator(opstr,SpinOp_1);
  end
end
end

%==========================================================================
% Convert a string to a cluster spin-operator.
%==========================================================================
function spinOp = assembleSpinOperator(opstr,SpinOp_1)

% ENUM
E = 1;  Z = 2; RAISE = 3; LOWER = 4;

% Initialize output.
spinOp = 1;

% Loop of single-spin operators in the provided string.
for ichar = 1:numel(opstr)
  
  % Identify single-spin operator;
  switch opstr(ichar)
    case 'E'
      op_ = E;
    case 'Z'
      op_ = Z;
    case 'R'
      op_ = RAISE;
    case 'L'
      op_ = LOWER;
    otherwise
      disp(['Unidentified reference ', opstr(ichar), '.']);
      error('Could not identify spin operator.');
  end
  
  % Krocker the operator into the output.
  spinOp = kron(spinOp,SpinOp_1(:,:,op_));
end

end


