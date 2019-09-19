%test_TEMPO
%clear;
%clf;
%clc;
InputData = 'TEMPO.pdb';
%InputData = 'TEMPO_solvated_.pdb';
%InputDat
% HETATM   28  N   LIG     1       3.347  10.998   7.746  1.00  0.00           N
R_N = [3.347  10.998   7.746]*1e-10; % m.
% HETATM   29  O   LIG     1       4.293  10.190   8.078  1.00  0.00           O
R_O = [4.293  10.190   8.078 ]*1e-10; % m.

System.Electron.Coordinates = (R_N+R_O)/2;
System.Time = linspace(0,50000,1000)*1e-6; % s.
System.radius = 20*1e-10; % m;
%System.averaging = 'none';
%System.D2O = true;
Method.method = 'restrictedCE';
Method.r0 = 5e-10; % m.
Method.order = 2;
Method.verbose = true;

[SignalMean, Signals] = nuclear_spin_diffusion(System,Method,InputData);

hold on;
%%{
for isignal = 1:length(Signals)
  plot(2*System.Time*1e6,Signals{isignal},'color', [0.4, 0.0, 0.0,0.1]);
end
%}
plot(2*System.Time*1e6,SignalMean,'-', 'color', [1,0,0]);
%semilogx(2*System.Time*1e6,SignalMean,'-', 'color', [1,0,0]);

xlabel('2\tau (\mus)');
ylabel('v/v_0');