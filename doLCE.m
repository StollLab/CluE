function [Signal,AuxiliarySignal,Order_n_Signal] = ...
  doLCE(System,Method, Nuclei, Clusters, verbose);

MethylID = [];

Signal = zeros(size(System.Time));
V2 = Signal;
Order_n_Signal{1} = Signal;
tau1 = System.Time;
tau2 = System.Time';

if Method.conserveMemory
  AuxiliarySignal = 'The auxiliary signals are not saved in memory conservation mode.';
end

for isize = 2
  for icluster = 1:Nuclei.numberClusters(isize)
    
    thisCluster = Clusters{isize}(icluster,:);
    
    [tensors,zeroIndex] = pairwisetensors_gpu(Nuclei.Nuclear_g, Nuclei.Coordinates,...
      thisCluster,Nuclei.Atensor, System.magneticField, System.ge, System.gMatrix(3,3), System.muB, System.muN, System.mu0, System.hbar,System.theory,MethylID);
    
    b = -2*pi*tensors(3,3,2,3)/4; % rad/s.
    omega = 2*pi*(tensors(3,3,1,2) -  tensors(3,3,1,3))/2; % rad/s.
    
    
    switch System.experiment
      case 'Hahn'
      case 'CPMG-2D'
        V2 = - 2*(b/omega*(cos(omega.*tau1) - cos(omega.*tau2)  )).^2;
        AuxiliarySignal_  = V2;
        Signal = Signal + AuxiliarySignal_;
    end
    % initialize the nth auxiliary signal
    if ~Method.conserveMemory
      AuxiliarySignal{icluster} = AuxiliarySignal_;
    end
  end
  

  
  Signal = exp(Signal);
  if System.dimensionality ==2
    Signal = reshape(Signal.',1,[]);
  end
  Order_n_Signal{1} = Signal;
  Order_n_Signal{isize} = Signal;
  
  
end

