% function test_centralSpinSystem_isotopologue()
clear;
clear centralSpinSystem;
clc;

Data.InputData = 'TEMPO_Gly_70A.pdb';

P = 1e-1;

System.radius = (P^(-1/3))*12e-10; % m.
System.Electron.Coordinates = {28,29};

System.timepoints = 2^7;
System.dt = 0.5e-6; % s
System.carbon = false;

System.particleOptions = {...
  'hydrogen','TEM', 'abundance', 0, ...
  'hydrogen','TEM', 'active', false, ...
  'hydrogen','SOL','switchParticle', 'void',...
  'hydrogen','SOL', 'abundance', P, ...
  'hydrogen','MGL','switchParticle', 'void', ...
  'hydrogen','MGL', 'abundance', P, ...
  };

System.deuterateProtein = true;
System.D2O = true;
System.deuteriumFraction = 1-P;

System.nitrogen = false;
System.spinHalfOnly = true;
Method.reparseNuclei = true;

Method.Criteria = {'dipole'};
Method.cutoff.dipole = 10^3; % Hz
Method.getNuclearContributions = true;

[System,Method, Data, ~] = setSystemDefaults(System,Method, Data);

numTrial = 10;
N0 = zeros(numTrial,1);
N = zeros(numTrial,1);
for ii=1:numTrial
  disp(['trial ', num2str(ii), '/', num2str(numTrial)]);
  if numTrial>1
    Nuclei0 = parseNuclei(System,Method,Data,Data.InputData);
    N0(ii) = Nuclei0.number;
  end
  Nuclei = centralSpinSystem(System,Method,Data);
  N(ii) = Nuclei.number;
end

%%
if numTrial>1
n0 = mean(N0);
N1 = N-19;
n = mean(N1);

sig0 = sqrt( sum( (N0-n0).^2 ) / (numTrial-1) );
sig  = sqrt( sum( (N1 -n ).^2 ) / (numTrial-1) );

Delta_n = (n0-n)
Delta_n/sig0

x = linspace( n0 - 5*sig0,n0 + 5*sig0,1001);
g0 = 1/sig0/sqrt(2*pi) * exp( -(x -n0).^2 /( 2* sig0^2 ) );
g  = 1/sig/sqrt(2*pi)  * exp( -(x -n ).^2 /( 2* sig^2  ) );

G0 = 1/sig0/sqrt(2*pi) * exp( -(N0-n0).^2 /( 2* sig0^2 ) );
G  = 1/sig/sqrt(2*pi)  * exp( -(N1 -n ).^2 /( 2* sig^2  ) );


clf
hold on;
color0 = [0,0,0];
color = [1,0,0];
plot(x,g0,'color',color0);
plot(x,g,'color',color);
xline(n0,'color',color0);
xline(n0-sig0,'--','color',color0);
xline(n0+sig0,'--','color',color0);
xline(n,'color',color);
xline(n-sig,'--','color',color);
xline(n+sig,'--','color',color);

plot(N0,exp(-1)*ones(size(N0))/sig0/sqrt(2*pi),'o','color',color0);
plot(N1,exp(-2)*ones(size(N1))/sig/sqrt(2*pi),'o','color',color);

xlabel('number')


end

