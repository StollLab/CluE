%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function tests = test_pairwiseTensors
tests = functiontests(localfunctions);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function setupOnce(testCase)
% oldpath = path;
% path('../',oldpath);
% testFunctionName = 'pairwiseTensors';

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function test_h18TEMPO(testCase)

testFunctionName = 'pairwiseTensors';
molName = 'h18TEMPO';
Data.writeSpinPDB = false;
testName = [testFunctionName,'_',molName];

Data.InputData = 'assets/TEMPO.pdb';

System.Electron.Coordinates = {28, 29};
System.X = {28, 29};
System.g = [2.0097, 2.0064,2.0025];
System.spinCenter = 'TEMPO';
System.magneticField  = 1.2; % T.
System.radius = 12e-10;
System.timepoints = 2^7;
System.dt = 0.05*1e-6; % s.

Method.Criteria = {'dipole'};
Method.cutoff.dipole = 0;

[System, Method, Data] = setSystemDefaults(System,Method,Data);
pdb = parsePDBfile(Data.InputData, System.angstrom);
[Nuclei, System] = centralSpinSystem(System,Method,Data,pdb);

Cluster = 1:Nuclei.number;
tensors = pairwisetensors_gpu(...
  Nuclei.Nuclear_g,...
  Nuclei.Coordinates,...
  Cluster,...
  Nuclei.Atensor,...
  System.magneticField,...
  System.ge,...
  System.g(3),...
  System.muB, ...
  System.muN, ...
  System.mu0, ...
  System.hbar,...
  System.theory, ...
  0, 0, 0, [], []);

  N = Nuclei.number+1;
  T_zz_neighbor = zeros(N*(N-1)/2,1);
  T_zz_vertex = zeros(N,1);

  vCounter = 0;
  nCounter = 0;
  for ii = 1:Nuclei.number+1
    vCounter = vCounter + 1;
    T_zz_vertex(vCounter) = tensors(3,3,ii,ii);
    for jj = ii+1:Nuclei.number+1
      nCounter = nCounter + 1;
      T_zz_neighbor(nCounter)  = tensors(3,3,ii,jj);
    end
  end

  T = array2table([T_zz_neighbor]);                               
  T.Properties.VariableNames(1) = {'T_zz'};      

  outFileName = ['output/output_',testName,'_neighbor' , '.csv'];
  writetable(T,outFileName);

  verified_Tn = readmatrix(['verified/verified_',testName,'_neighbor', '.csv']);
  error_n = max(abs(verified_Tn - T_zz_neighbor));

  tol = max(abs(verified_Tn))*1e-12;

  verifyEqual(testCase,error_n,0,'AbsTol',tol)

  T = array2table([T_zz_vertex]);                               
  T.Properties.VariableNames(1) = {'T_zz'};      
  outFileName = ['output/output_',testName,'_vertex' , '.csv'];
  writetable(T,outFileName);

  verified_Tv = readmatrix(['verified/verified_',testName,'_vertex', '.csv']);
  error_v = max(abs(verified_Tv - T_zz_vertex));
 
  tol = max(abs(verified_Tv))*1e-12;
  verifyEqual(testCase,error_v,0,'AbsTol',tol)
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function test_d18TEMPO(testCase)

testFunctionName = 'pairwiseTensors';
molName = 'd18TEMPO';
testName = [testFunctionName,'_',molName];

Data.InputData = 'assets/TEMPO.pdb';
Data.writeSpinPDB = false;
System.particleOptions = {'hydrogen','TEM', 'abundance', 0};
System.Electron.Coordinates = {28, 29};
System.X = {28, 29};
System.g = [2.0097, 2.0064,2.0025];
System.spinCenter = 'TEMPO';
System.magneticField  = 1.2; % T.
System.radius = 12e-10;
System.timepoints = 2^7;
System.dt = 0.05*1e-6; % s.

Method.Criteria = {'dipole'};
Method.cutoff.dipole = 0;

[System, Method, Data] = setSystemDefaults(System,Method,Data);
pdb = parsePDBfile(Data.InputData, System.angstrom);
[Nuclei, System] = centralSpinSystem(System,Method,Data,pdb);

Cluster = 1:Nuclei.number;
tensors = pairwisetensors_gpu(...
  Nuclei.Nuclear_g,...
  Nuclei.Coordinates,...
  Cluster,...
  Nuclei.Atensor,...
  System.magneticField,...
  System.ge,...
  System.g(3),...
  System.muB, ...
  System.muN, ...
  System.mu0, ...
  System.hbar,...
  System.theory, ...
  0, 0, 0, [], []);

  N = Nuclei.number+1;
  T_zz_neighbor = zeros(N*(N-1)/2,1);
  T_zz_vertex = zeros(N,1);

  vCounter = 0;
  nCounter = 0;
  for ii = 1:Nuclei.number+1
    vCounter = vCounter + 1;
    T_zz_vertex(vCounter) = tensors(3,3,ii,ii);
    for jj = ii+1:Nuclei.number+1
      nCounter = nCounter + 1;
      T_zz_neighbor(nCounter)  = tensors(3,3,ii,jj);
    end
  end

  T = array2table([T_zz_neighbor]);                               
  T.Properties.VariableNames(1) = {'T_zz'};      

  outFileName = ['output/output_',testName,'_neighbor' , '.csv'];
  writetable(T,outFileName);

  verified_Tn = readmatrix(['verified/verified_',testName,'_neighbor', '.csv']);
  error_n = max(abs(verified_Tn - T_zz_neighbor));
  tol = max(abs(verified_Tn))*1e-12;
  verifyEqual(testCase,error_n,0,'AbsTol',tol)

  T = array2table([T_zz_vertex]);                               
  T.Properties.VariableNames(1) = {'T_zz'};      
  outFileName = ['output/output_',testName,'_vertex' , '.csv'];
  writetable(T,outFileName);

  verified_Tv = readmatrix(['verified/verified_',testName,'_vertex', '.csv']);
  error_v = max(abs(verified_Tv - T_zz_vertex));
 
  tol = max(abs(verified_Tv))*1e-12;
  verifyEqual(testCase,error_v,0,'AbsTol',tol)
end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function teardownOnce(testCase) 
% path(oldpath)
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

