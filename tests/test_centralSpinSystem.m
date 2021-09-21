% function test_centralSpinSystem()
clear;
clear centralSpinSystem;
clc;

Data.InputData = 'mTEMPO_Gly_70A.pdb';

System.radius = 12e-10; % m.
System.Electron.Coordinates = {28,29};

System.timepoints = 2^7;
System.dt = 0.5e-6; % s
System.carbon = false;

System.particleOptions = {...
  'H','TEM','flag', 'value',...
  'H','TEM','flag', 'value'};

Method.Criteria = {'dipole'};
Method.cutoff.dipole = 10^3; % Hz
Method.getNuclearContributions = true;
Method.getNuclearStatistics = true;

Method.Ori_cutoffs = true;

[System,Method, Data, ~] = setSystemDefaults(System,Method, Data);
Nuclei0 = parseNuclei(System,Method,Data,Data.InputData);
[~, System]= centralSpinSystem(System,Method,Data);
[Nuclei, System]= centralSpinSystem(System,Method,Data);

%% Check number
if Nuclei0.number~=Nuclei.number
  clc
  disp('    1H   13C   14N')
  disp(' ')

  H=1; C =2;  N=3;
  num0 = [0,0,0];
  
  for ii = 1:Nuclei0.number
    if strcmp(Nuclei0.Type{ii},'1H')
      num0(H) = num0(H) +1;
    elseif strcmp(Nuclei0.Type{ii},'14N')
      num0(N) = num0(N) +1;
    elseif strcmp(Nuclei0.Type{ii},'13C')
      num0(C) = num0(C) +1;
    end
  end
  disp(num0);
  
  num = [0,0,0];
  
  for ii = 1:Nuclei.number
    if strcmp(Nuclei.Type{ii},'1H')
      num(H) = num(H) +1;
    elseif strcmp(Nuclei.Type{ii},'14N')
      num(N) = num(N) +1;
    elseif strcmp(Nuclei.Type{ii},'13C')
      num(C) = num(C) +1;
    end
  end
  disp(num);
  
  error(['Error in test_centralSpinSystem(): ', ...
    'parseNuclei() generates ', num2str(Nuclei0.number), ' particles, but ',...
    'centralSpinSystem() generates ', num2str(Nuclei.number), ' particles.'])
else
  disp('number check: pass')
end
%% Check hydron types
if Nuclei0.number_1H_exchangeable == Nuclei.number_1H_exchangeable
  disp('number 1H_exchangeablecheck: pass');
else
  disp('number 1H_exchangeablecheck: fail');
end

if Nuclei0.number_1H_nonExchangeable == Nuclei.number_1H_nonExchangeable
  disp('number 1H_nonExchangeablecheck: pass');
else
  disp('number 1H_nonExchangeablecheck: fail');
end

if Nuclei0.number_2H_exchangeable == Nuclei.number_2H_exchangeable
  disp('number 2H_exchangeablecheck: pass');
else
  disp('number 2H_exchangeablecheck: fail');
end

if Nuclei0.number_2H_nonExchangeable == Nuclei.number_2H_nonExchangeable
  disp('number 2H_nonExchangeablecheck: pass');
else
  disp('number 2H_nonExchangeablecheck: fail');
end

%% pdb ID check
[~,idx0] = sort(Nuclei0.pdbID);
[~,idx] = sort(Nuclei.pdbID);
Delta_ID = Nuclei0.pdbID(idx0) - Nuclei.pdbID(idx);
if max(abs(Delta_ID(:)) ) >1e-12
  plot(Delta_ID)
  xlabel('index')
  ylabel('\DeltaID')
  disp('pdb ID: fail')
else
  disp('pdb ID: pass')
end

%% g factor check
Delta_g = Nuclei0.Nuclear_g(idx0) - Nuclei.Nuclear_g(idx);
if max(abs(Delta_g(:)) ) >1e-3
  plot(Delta_g)
  xlabel('index')
  ylabel('\Deltag')
  disp('g factor: fail')
else
  disp('g factor: pass')
end

%% StateMultiplicity check
Delta_StateMultiplicity = Nuclei0.StateMultiplicity(idx0) ...
  - Nuclei.StateMultiplicity(idx);
if max(abs(Delta_StateMultiplicity(:)) ) >1e-12
  plot(Delta_StateMultiplicity)
  xlabel('index')
  ylabel('\DeltaStateMultiplicity')
  disp('StateMultiplicity: fail')
else
  disp('StateMultiplicity: pass')
end
%% Spin check
Delta_Spin = Nuclei0.Spin(idx0) - Nuclei.Spin(idx);
if max(abs(Delta_Spin(:)) ) >1e-12
  plot(Delta_Spin)
  xlabel('index')
  ylabel('\DeltaSpin')
  disp('Spin: fail')
else
  disp('Spin: pass')
end
%% Index check
Same_Index = Nuclei0.Index == Nuclei.Index;
if any(~Same_Index)
  plot(Same_Index)
  xlabel('index')
  ylabel('\DeltaIndex')
  disp('Index: fail')
else
  disp('Index: pass')
end
%% Check Coordinates
DeltaCoor = Nuclei0.PDBCoordinates(idx0,:) - Nuclei.PDBCoordinates(idx,:);
err = vecnorm(DeltaCoor');
if max(err(:))>1e-12
  plot(err)
  xlabel('index')
  ylabel('|\DeltaR| (m)')
  disp('pdb coor: fail')
else
  disp('pdb coor: pass')
end
%% Check  Electron_pdbCoordinates

if vecnorm(Nuclei0. Electron_pdbCoordinates ...
    - Nuclei. Electron_pdbCoordinates) > 1e-12
  disp('electron coor: fail')
else
  disp('electron coor: pass')
end

%% Check Coordinates
DeltaCoor = Nuclei0.Coordinates(idx0,:) - Nuclei.Coordinates(idx,:);
err = vecnorm(DeltaCoor');
if max(err(:))>1e-12
  plot(err)
  xlabel('index')
  ylabel('|\DeltaR| (m)')
  disp('coor: fail')
else
  disp('coor: pass')
end

%% Check valid nuclei
matchValid =  Nuclei0.valid(idx0) == Nuclei.valid(idx);
if any(~matchValid)
  disp('valid nuclei: fail')
else
  disp('valid nuclei: pass')
end

%%
DeltaDistanceMatrix = Nuclei0.Statistics.DistanceMatrix(idx0,idx0) ...
  - Nuclei.Statistics.DistanceMatrix(idx,idx);

if max(abs(DeltaDistanceMatrix(:) ))>0
  disp('DistanceMatrix: fail')
else
  disp('DistanceMatrix: pass')
end

%%
DeltaAdjacency = Nuclei0.Adjacency(idx0,idx0,1) ...
  - Nuclei.Adjacency(idx,idx,1);

if max(abs(DeltaAdjacency(:) ))>0
  disp('Adjacency: fail')
else
  disp('Adjacency: pass')
end