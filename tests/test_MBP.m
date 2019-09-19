  clear;
  clf;
  clc;
  
  %InputData = '1n3w.pdb';
  %InputData = 'mbp_bound.pdb';
  %InputData = 'mbp_bound_solvated_10A.pdb';
  %InputData = '1lls_solv.pdb';
  %InputData = '1lls_gmx.pdb';
  InputData = '1lls_gmx_h_add.pdb';
  System.Time = linspace(0,10,50)*1e-6; % s.
  System.Electron.Coordinates = 2; % atom 12.
  %System.g = [2,2,2];
  System.averaging = 'none';
  
  Method.method = 'restrictedCE';
  Method.r0 = 10e-10; % m.
  Method.order = 2;
  Method.verbose = true;
  
  [SignalMean, Signals] = nuclear_spin_diffusion(System,Method,InputData);
  
  System.hydrogen = false;
  Method.method = 'CE';
  hold on;
  for isignal = 1:length(Signals)
    plot(2*System.Time*1e6,Signals{isignal},'color', [0, 0.0, 0.0,0.1]);
  end
  plot(2*System.Time*1e6,SignalMean,'-o', 'color', [0,0,0]);
  xlabel('2\tau (\mus)');
  ylabel('v/v_0');