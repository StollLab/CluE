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
System.spinHalfOnly = ~true;
Method.reparseNuclei = true;

Method.Criteria = {'dipole'};
Method.cutoff.dipole = 10^3; % Hz
Method.getNuclearContributions = true;

[System,Method, Data, ~] = setDefaults(System,Method, Data);

numTrial = 1000;
N0 = zeros(numTrial,1);
N = zeros(numTrial,1);
NN0 = zeros(numTrial,1);
NN = zeros(numTrial,1);
NNex = zeros(numTrial,1);
NNnx = zeros(numTrial,1);
NNex0 = zeros(numTrial,1);
NNnx0 = zeros(numTrial,1);

pdb0 = parsePDB(Data.InputData,System);
pdb_ = parsePDBfile(Data.InputData, System.angstrom);
tic
parfor ii=1:numTrial
  disp(['trial ', num2str(ii), '/', num2str(numTrial)]);
  if numTrial>1
    Nuclei0 = parseNuclei(System,Method,Data,pdb0);
    Nuclei0 = newHydronIsotopologue(Nuclei0,System);
    NN0(ii) = Nuclei0.number;
    N0(ii) = sum(strcmp(Nuclei0.Type(:),'1H'));
    NNex0(ii) = Nuclei0.number_1H_exchangeable;
    NNnx0(ii) = Nuclei0.number_1H_nonExchangeable;
    
    NNN0 = Nuclei0.Isotopologue.TypeNumber-...
      Nuclei0.Isotopologue.Instance_2H_Number;
    
    N0(ii) = NNN0(1);
    NNex0(ii) = NNN0(2);
    NNnx0(ii) = NNN0(3);
    
  end
  Nuclei = centralSpinSystem(System,Method,Data,pdb_);
  NN(ii) = Nuclei.number;
  N(ii) = sum(strcmp(Nuclei.Type(:),'1H'));
  NNex(ii) = Nuclei.number_1H_exchangeable;
  NNnx(ii) = Nuclei.number_1H_nonExchangeable;
end
toc
%%
clf
doSaveFig = true
if doSaveFig
  figure();
  fontsize = 12;
  lw = 1;
  sz = 16;
  cbscale = 200;
else
  fontsize = 24;
  lw = 2;
  sz = 50;
  cbscale = 100;
end
for dataOption = 0:2
  subplot(1,3,dataOption+1);
switch dataOption
  case 0
    D1 = NNex; %-19;
    D0 = NNex0;
    name = 'exchangeable 1H';
  case 1
    D1 = NNnx;% -19;
    D0 = NNnx0;
    name = 'non-exchangeable 1H';
  case 2
    
    D1 = NN -19;
    D0 = N0;
    name = 'all 1H';
end

if numTrial>50

% 
% D1 = NNnx;% -19;
% D0 = NNnx0;



% D1 = N; % -19;
% D0 = N0;
% clf
hold on;
color0 = [0,0,0];
color = [1,0,0];

edges = 0:2:1000;
histogram(D0,edges)
histogram(D1,edges)

xlim([ min([D0;D1]), max([D0;D1])]);
title(name);
xline(mean(D0),'-', 'color',[0 0 0], 'linewidth',2*lw)
xline(mean(D0),'-', 'color',[0 0 0]+1, 'linewidth',lw)
xline(mean(D0),'--', 'color',[0 0.4470 0.7410], 'linewidth',lw)
xline(mean(D1),'-','color',[0 0 0], 'linewidth',2*lw)
xline(mean(D1),'-','color',[0 0 0]+1, 'linewidth',lw)
xline(mean(D1),'--','color',[0.8500 0.3250 0.0980], 'linewidth',lw)
legend('old','new')
xlabel('number of protons')
ylabel('count')
ax = gca;
ax.XAxis.Color = [0,0,0];
ax.YAxis.Color = [0,0,0];
set(gca,'fontsize',fontsize);
box on
elseif numTrial>1
% D1 = NNex; %-19;
% D0 = NNex0;
% 
% D1 = NNnx;% -19;
% D0 = NNnx0;

% D1 = NN -19;
% D0 = N0;

% D1 = N; % -19;
% D0 = N0;

n0 = mean(D0);
n = mean(D1);

sig0 = sqrt( sum( (D0-n0).^2 ) / (numTrial-1) );
sig  = sqrt( sum( (D1 -n ).^2 ) / (numTrial-1) );

Delta_n = (n0-n)
Delta_n/sig0

x = linspace( n0 - 3*sig,n0 + 3*sig,1001);
g0 = 1/sig0/sqrt(2*pi) * exp( -(x -n0).^2 /( 2* sig0^2 ) );
g  = 1/sig/sqrt(2*pi)  * exp( -(x -n ).^2 /( 2* sig^2  ) );

G0 = 1/sig0/sqrt(2*pi) * exp( -(D0-n0).^2 /( 2* sig0^2 ) );
G  = 1/sig/sqrt(2*pi)  * exp( -(D1 -n ).^2 /( 2* sig^2  ) );


clf
hold on;
% pd0 = makedist('Binomial','N',round(n0/P),'p',P);
% pd = makedist('Binomial','N',round(n/P),'p',P);
% pd5 = makedist('Binomial','N',round(n/P/5),'p',P);
% % h0 = histcounts(random(pd,1,round(n0/P)));
% % plot(h0/max(h0));
% counts0 = cell(10,1);
% bins0 = cell(10,1);
% counts = cell(10,1);
% bins = cell(10,1);
% counts5 = cell(10,1);
% bins5 = cell(10,1);
% for ii = 1:10
% h = histogram(random(pd0,1, round(n0/P) ) );
% counts0{ii} = h.Values;
% bins0{ii} = h.BinEdges;
% h = histogram(random(pd,1, round(n/P) ) );
% counts{ii} = h.Values;
% bins{ii} = h.BinEdges;
% 
% h = histogram(random(pd5,1, round(n/P) ) );
% counts5{ii} = h.Values;
% bins5{ii} = h.BinEdges;
% end
% clf
% hold on;
% for ii = 1:10
%   plot(bins0{ii}(1:end-1),counts0{ii}/max(counts0{ii})/sig0/sqrt(2*pi),...
%     'color', [1,1,1,1]/4 );
%   
%   plot(bins{ii}(1:end-1),counts{ii}/max(counts{ii})/sig/sqrt(2*pi),...
%     'color', [1,0,0,1/4] );
%   
%   plot(5*bins5{ii}(1:end-1),counts5{ii}/max(counts5{ii})/sig/sqrt(2*pi),...
%     'color', [0,0,1,1/4] );
% end

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

plot(D0,exp(-1)*ones(size(D0))/sig0/sqrt(2*pi),'o','color',color0);
plot(D1,exp(-2)*ones(size(D1))/sig/sqrt(2*pi),'o','color',color);

legend('old','new')

xlabel('number')
xlim([n0-3*sig,n0 + 3*sig])
title(name);
end
end
if doSaveFig
figName = 'fig_test_centralSpinSystem_isotopologue';
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 3*1.2 0.9]*9);
% set(gcf,'PaperUnits','points','PaperPosition',[0 0 1 3/4]*9/2.54);
print(figName,'-dpng','-r600');
% print(['SVG/',figName],'-dsvg','-r600');

 set(gcf,'WindowStyle','Docked')
end