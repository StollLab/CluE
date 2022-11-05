%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function tests = test_findClusters
tests = functiontests(localfunctions);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function test_cube(testCase)
%                                                                               
% cube connectivity                                                             
%                                                                               
%    5-------6                                                                  
%    | \   / |                                                                  
%    |  1-2  |                                                                  
%    |  | |  |                                                                  
%    |  4-3  |                                                                  
%    | /   \ |                                                                  
%    8-------7                                                                  
%                                                                               
% Correct cluster counts for a cube adjacency matrix.                           
% Cluster size:            1, 2, 3, 4, 5, 6, 7, 8                       
correctNumberClusters = [  8,12,24,38,48,28, 8, 1];

Nuclei.Adjacency = [...
  1,1,0,1,1,0,0,0; ...                                                        
  1,1,1,0,0,1,0,0; ...                                                        
  0,1,1,1,0,0,1,0; ...                                                        
  1,0,1,1,0,0,0,1; ...                                                       
  1,0,0,0,1,1,0,1; ...                                                       
  0,1,0,0,1,1,1,0; ...                                                       
  0,0,1,0,0,1,1,1; ...                                                        
  0,0,0,1,1,0,1,1];

Nuclei.number = 8;
Nuclei.AntiAdjacency = false(Nuclei.number,Nuclei.number);
Nuclei.valid = true(Nuclei.number,1);
Nuclei.Type = {'1H','1H','1H','1H','1H','1H','1H','1H'};


Method.graphCriterion = 'connected';
Method.emptyClusterSetsOkay = false;
Method.Criteria = {}; 
Method.includeAllSubclusters = false;

Clusters = findClusters_treeSearch(Nuclei,8,1,...         
  [], Method);

numberClusters = zeros(1,8);
for isize = 1:8
  numberClusters(isize) = size(Clusters{isize},1);
end

pass = all(numberClusters==correctNumberClusters);

verifyEqual(testCase,pass,true);
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>