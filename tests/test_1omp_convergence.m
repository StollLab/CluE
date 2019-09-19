function test_1omp_convergence()
  clear;
  clc;
  InputData_solvent = '1omp_res287_to_R1_cryH2O_addH_40A_solvent_.pdb';  
  %InputData = '1omp_res287_to_R1_cryH2O_addH.pdb';
  System.Time = linspace(0,10,50)*1e-6; % s.
  R_O15 = [-0.127  17.489 -19.006]*1e-10; % m.
  R_N1 = [-0.020  16.336 -18.471]*1e-10; % m.
  System.Electron.Coordinates = (R_N1+R_O15)/2;
  System.averaging = 'none';
  Method.method = 'restrictedCE';
  Method.r0 = 10e-10; % m.
  Method.order = 2;
  Method.verbose = true;lse;
  save('data_1omp_dueterated_distance.mat');
  % D D  
  exp = 'H D';
  System.radius = 30e-10; % m.
  System.D20 = true; 
  %System.deuterateProtein = true;
  InputData = InputData_solvent; 
  Method.method = 'CE';
  [SignalMean_HD_r30, Signals_HD_r30] = nuclear_spin_diffusion(System,Method,InputData);
  save('data_1omp_dueterated_distance.mat');
  % D D  
  exp = 'D D';
  System.radius = 20e-10; % m.
  System.D20 = true; 
  System.deuterateProtein = true;
  InputData = InputData_solvent; 
  Method.method = 'CE';
  [SignalMean_DD_r20, Signals_DD_r20] = nuclear_spin_diffusion(System,Method,InputData);
  save('data_1omp_dueterated_distance.mat');
  % D D  
  exp = 'D D';
  System.radius = 30e-10; % m.
  System.D20 = true; 
  System.deuterateProtein = true;
  InputData = InputData_solvent; 
  Method.method = 'CE';
  [SignalMean_DD_r30, Signals_DD_r30] = nuclear_spin_diffusion(System,Method,InputData);
  save('data_1omp_dueterated_distance.mat');
end