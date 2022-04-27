%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function tests = test_getAdjacencyMatrix
tests = functiontests(localfunctions);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function test_h18TEMPO(testCase)

testFunctionName = 'getAdjacencyMatrix';
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

Method.order = 2;

[System, Method, Data] = setSystemDefaults(System,Method,Data);
pdb = parsePDBfile(Data.InputData, System.angstrom);
[Nuclei, System] = centralSpinSystem(System,Method,Data,pdb);
Nuclei.Statistics = getPairwiseStatistics(System, Method, Nuclei);


Method.Criteria = {'dipole'};
Method.cutoff.dipole = 1e3*[1,1];
Method.Ori_cutoffs = true;
Adjacency = getAdjacencyMatrix(System, Nuclei,Method);
adjacency_dipole = Adjacency(:,:,1);

T = array2table(adjacency_dipole);
for ii =1:numel(Nuclei.pdbID)
  T.Properties.VariableNames(ii) = {['pdbID_', num2str(Nuclei.pdbID(ii))]};
end
writetable(T,'output/output_adjacency_dipole_1kHz.csv');

Method.Criteria = {'distance'};
Method.cutoff.rMax = 4e-10;
Method.cutoff.rMin = 0;
Method.Ori_cutoffs = true;
Adjacency = getAdjacencyMatrix(System, Nuclei,Method);
adjacency_distance = Adjacency(:,:,1);


T = array2table(adjacency_distance);
for ii =1:numel(Nuclei.pdbID)
  T.Properties.VariableNames(ii) = {['pdbID_', num2str(Nuclei.pdbID(ii))]};
end
writetable(T,'output/output_adjacency_distance_4A.csv');

verified_dipole = readmatrix('verified/verified_adjacency_dipole_1kHz.csv');
verified_distance = readmatrix('verified/verified_adjacency_distance_4A.csv');

error_dipole = abs(verified_dipole - adjacency_dipole);
error_distance = abs(verified_distance- adjacency_distance);

verifyEqual(testCase,sum(error_dipole(:)),0,'AbsTol',1e-6)
verifyEqual(testCase,sum(error_distance(:)),0,'AbsTol',1e-6)



end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>