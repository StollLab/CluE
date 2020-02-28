
  clear;
  clc;
  
  InputData = '116253_unitcell.pdb'; 
  % System ================================================================
  System.experiment = 'Hahn';
  System.Time = linspace(0,10,50)*1e-6; % s.
  System.Electron.Coordinates = 32; % atom #32. %[9.48360        2.99030        7.47299]*1e-10; % m.
  System.g = [2.051,2.051,2.258];
  System.radius = 10e-10; % m.
  System.Y = 46; % atom #46. %[10.506   3.320   9.057]*1e-10-System.Electron.Coordinates; %R_O46 - R_el.
  System.X = 37; % atom # 37. %[10.511   1.417   7.252]*1e-10-System.Electron.Coordinates; %R_O37 - R_el.
  %System.averaging = 'xy';
  % Method ================================================================
  Method.method = 'restrictedCE';
  Method.r0 = 10e-10; % m.
  Method.order = 2;
  Method.verbose = true;
  System.averaging = 'xy';
  
  
  
  % Simulation ============================================================
  [SignalMean, Signals] = CluE(System,Method,InputData);
  % Plot ==================================================================
  hold on;
  LenzData=csvread('cudbm_Lenz_fit_data.csv',1,0);
  LenzData(:,2) = LenzData(:,2) - min(LenzData(:,2));
  LenzData(:,2) = LenzData(:,2)/max(LenzData(:,2));  
  plot(LenzData(:,1),LenzData(:,2),'- red');
  Lsys.tau = System.Time*1e6;
  Lsys.TM = 9.18;%7.74; % us.
  Lsys.k = 2.68;
  Lenz_v = hahn(Lsys);
  plot(2*System.Time*1e6,Lenz_v, '-- red');
  for isignal = 1:length(Signals)
    plot(2*System.Time*1e6,Signals{isignal},'color', [0, 0.0, 0.0,0.1]);
  end
  plot(2*System.Time*1e6,SignalMean,'-o', 'color', [0,0,0]);
  xlabel('2\tau (\mus)');
  ylabel('v/v_0');
