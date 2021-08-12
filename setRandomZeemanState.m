%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% This function generates an initial bath state.
function [ZeemanState, spinState] = setRandomZeemanState(Nuclei,nStates)
if nargin==1
  nStates = Nuclei.nStates;
end
maxnStates = max(nStates);
ZeemanState = zeros(maxnStates,Nuclei.number);
spinState = ZeemanState;

% Loop through all nuclei.
for iSpin = 1:Nuclei.number
  I = Nuclei.Spin(iSpin);
  ZeemanState(:,iSpin) =  randi(2*I+1,maxnStates,1); 
  spinState(:,iSpin) = ZeemanState(:,iSpin) -I - 1;
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>