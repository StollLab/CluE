% function test_parsePDBfile()
% Data.InputData = 'TEMPO_Gly_70A.pdb';
Data.InputData = 'CSD_MEMALA01.pdb';

System.radius = 12e-10; % m.
System.Electron.Coordinates = {28,29};

System.timepoints = 2^7;
System.dt = 0.5e-6; % s
System.carbon = false;

Method.Criteria = {'dipole'};
Method.cutoff.dipole = 10^3; % Hz
Method.getNuclearContributions = true;

[System,Method, Data, ~] = setSystemDefaults(System,Method, Data);


pdbData = parsePDBfile(Data.InputData, System.angstrom);


% [pdbData0.pdbCoordinates,...
%   pdbData0.Type,pdbData0.UnitCell,pdbData0.Connected,...
%   pdbData0.Indices_nonSolvent,pdbData0.pdbID,pdbData0.MoleculeID,...
%   pdbData0.numberH,pdbData0.isSolvent,pdbData0.isWater,pdbData0.Exchangeable,...
%   pdbData0.VanDerWaalsRadii] ...
  pdbData0 = parsePDB(Data.InputData,System);

%% Check number
number0 = numel(pdbData0.Type);
if pdbData.number ~= number0
  error(['Error in tesr_parsePDBfile(): ',...
    'parsePDFfile() found ', num2str(pdbData.number),...
    ' nuclei, but parsePDB() found ', num2str(number0), ' nuclei.'])
else
  disp('number: pass');
end
%% Check Coordinates

coordinates = [pdbData.x,pdbData.y,pdbData.z];  
DeltaCoor = pdbData0.Coordinates - coordinates;
err = vecnorm(DeltaCoor');
if max(err(:))>1e-12
  plot(err)
  xlabel('index')
  ylabel('|\DeltaR| (m)')
  disp('pdb coor: fail')
else
  disp('pdb coor: pass')
end

%% Check ID
Delta_ID = pdbData0.pdbID - pdbData.serial';
if max(abs(Delta_ID(:)) ) >1e-12
  plot(Delta_ID)
  xlabel('index')
  ylabel('\DeltaID')
  disp('pdb ID: fail')
else
  disp('pdb ID: pass')
end