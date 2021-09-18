% function test_centralSpinSystem_unitcell()
clear;
clear centralSpinSystem;
clc;

Data.InputData = 'CSD_MEMALA01.pdb';

System.radius = 15e-10; % m.
System.Electron.Coordinates = {6};

System.timepoints = 2^7;
System.dt = 0.5e-6; % s
System.carbon = false;

System.particleOptions = {...
  'H','TEM','flag', 'value',...
  'H','TEM','flag', 'value'};

Method.Criteria = {'dipole'};
Method.cutoff.dipole = 10^3; % Hz
Method.getNuclearContributions = true;

[System,Method, Data, ~] = setSystemDefaults(System,Method, Data);
Nuclei0 = parseNuclei(System,Method,Data,Data.InputData);
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

%% pdb ID check
[~,idx0] = sort(Nuclei0.ucpdbID);
[~,idx] = sort(Nuclei.pdbID);
Delta_ID = Nuclei0.ucpdbID(idx0) - Nuclei.pdbID(idx);
if max(abs(Delta_ID(:)) ) >1e-12
  plot(Delta_ID)
 % plot(1:Nuclei0.number, Nuclei0.ucpdbID(idx0), ...
 %   1:Nuclei0.number, Nuclei.pdbID(idx))
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
r0 = Nuclei0.PDBCoordinates;
% [~,idx_x0] = sort(r0(:,1));  r0 = r0(idx_x0,:);
% [~,idx_y0] = sort(r0(:,2));  r0 = r0(idx_y0,:);
% [~,idx_z0] = sort(r0(:,3));  r0 = r0(idx_z0,:);

r = Nuclei.PDBCoordinates;
% [~,idx_x] = sort(r(:,1));  r = r(idx_x,:);
% [~,idx_y] = sort(r(:,2));  r = r(idx_y,:);
% [~,idx_z] = sort(r(:,3));  r = r(idx_z,:);

DeltaCoor = r0(idx0,:) - r(idx,:);

% DeltaCoor = Nuclei0.PDBCoordinates(idx0,:) - Nuclei.PDBCoordinates(idx,:);
err = vecnorm(DeltaCoor');
if max(err(:))>1e-12
  plot(err)
  xlabel('index')
  ylabel('|\DeltaR| (m)')
  disp('pdb coor: fail')
  
  N = Nuclei0.pdbNumber;
  
  if max(err(1:N))>1e-12
     disp('pdb coor in first unit cell: fail');
  else
     disp('pdb coor in first unit cell: pass');
  end

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

