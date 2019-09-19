  clear;
  clf;
  clc;
  %InputData = '1omp_Asp287R1_solvated.pdb';  
  InputData = '1omp_res287_to_R1_cryH2O_addH.pdb'; % D2O approximation
  %InputData = '1omp_res287_to_R1_cryH2O_addH_solvated.pdb'; % H2O 
  %InputData = '1omp_res287_to_R1_cryH2O_addH_40A_solvent_.pdb';
  System.Time = linspace(0,10,50)*1e-6; % s.
  R_O15 = [-0.127  17.489 -19.006]*1e-10; % m.
  R_N1 = [-0.020  16.336 -18.471]*1e-10; % m.
  System.Electron.Coordinates = (R_N1+R_O15)/2;
  System.radius = 10e-10; % m.
  System.averaging = 'none';
  %System.D2O = true;
  %System.deuterateAll = true;
  Method.method = 'CE';
  Method.r0 = 5e-10; % m.
  Method.order = 2;
  %Method.parallelComputing = false;
 %Method.verbose = true;
  
  %[SignalMean, Signals] = main_20180613(System,Method,InputData);
  [SignalMean, Signals] = nuclear_spin_diffusion(System,Method,InputData);
%   
%   System.hydrogen = false;
%   Method.method = 'CE';
%   [SignalMean2, Signals2] = nuclear_spin_diffusion(System,Method,InputData);
  hold on;
  for isignal = 1:length(Signals)
    plot(2*System.Time*1e6,Signals{isignal},'color', [0.4, 0.0, 0.0,0.1]);
  end
  plot(2*System.Time*1e6,SignalMean,'-', 'color', [1,0,0]);
  
%   plot(2*System.Time*1e6,SignalMean2,'-o', 'color', [1,0,0]);
  xlabel('2\tau (\mus)');
  ylabel('v/v_0');