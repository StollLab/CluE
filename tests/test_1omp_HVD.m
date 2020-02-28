function test_1omp_HVD()
  clear;
  clf;
  clc;
  InputData_solvent = '1omp_res287_to_R1_cryH2O_addH_25A_solvent_.pdb';
  System.Time = linspace(0,10,50)*1e-6; % s.
  R_O15 = [-0.127  17.489 -19.006]*1e-10; % m.
  R_N1 = [-0.020  16.336 -18.471]*1e-10; % m.
  System.Electron.Coordinates = (R_N1+R_O15)/2;
  System.radius = 20e-10; % m.
  System.averaging = 'none';
  Method.method = 'restrictedCE';
  Method.r0 = 10e-10; % m.
  Method.order = 2;
  Method.verbose = true;
  % H H
  exp = 'H H';
  InputData = InputData_solvent;  
  Method.method = 'restrictedCE';
  [SignalMean_HH, Signals_HH] = CluE(System,Method,InputData);
  save('data_1omp.mat');
  % H D  
  exp = 'H D';
  System.D20 = true; 
  InputData = InputData_solvent; 
  Method.method = 'CE';
  [SignalMean_HD, Signals_HD] = CluE(System,Method,InputData);
  System.D20 = false;
  save('data_1omp.mat');
  % D H  
  exp = 'D H';
  System.deuterateProtein = true;
  InputData = InputData_solvent;  
  Method.method = 'CE';
  [SignalMean_DH, Signals_DH] = CluE(System,Method,InputData);
  System.deuterateProtein = false;
  save('data_1omp.mat');
  % D D  
  exp = 'D D';
  System.D20 = true; 
  System.deuterateProtein = true;
  InputData = InputData_solvent; 
  Method.method = 'CE';
  [SignalMean_DD, Signals_DD] = CluE(System,Method,InputData);
  System.D20 = false;
  System.deuterateProtein = false;
  save('data_1omp.mat');
  % H ~
  exp = 'H ~';
  InputData = '1omp_res287_to_R1_cryH2O_addH.pdb';
  Method.method = 'restrictedCE';
  [SignalMean_Hn, Signals_Hn] = CluE(System,Method,InputData);
  save('data_1omp.mat');
  
  % D ~
  
  exp = 'D ~';
  System.deuterateProtein = true;
  InputData = '1omp_res287_to_R1_cryH2O_addH.pdb';
  Method.method = 'CE';
  [SignalMean_Dn, Signals_Dn] = CluE(System,Method,InputData);
  System.deuterateProtein = false;
  save('data_1omp.mat');
  % ~ H
  
  exp = '~ H';
  System.solventOnly = true;
  InputData = InputData_solvent; 
  Method.method = 'restrictedCE';
  [SignalMean_nH, Signals_nH] = CluE(System,Method,InputData);
  System.solventOnly = false; 
  save('data_1omp.mat');
  
  % ~ D
  exp = '~ D';
  System.solventOnly = true;
  System.D20 = true;
  InputData = InputData_solvent; 
  Method.method = 'CE';
  [SignalMean_nD, Signals_nD] = CluE(System,Method,InputData);
  System.solventOnly = false; 
  System.D2O = false;
  save('data_1omp.mat');

  
end