function pass = test_findClusters()
pass = true;
oldpath = path;
path('../',oldpath);

%==========================================================================
%  2 x 2 x 2 cube 
%==========================================================================

side = 2;
threshold = 1;

Nuclei.ValidPair = getAdjacencyMatrix(side,threshold);
Nuclei.number = side^3;        
cluster2 = findClusters(Nuclei,2);
% 1
if norm(size(cluster2)-[12,2]) ~= 0
  pass = false;
  fprintf('\nfindClusters() --> failed at test 1.\n');
  return;
end
% 2
Nuclei.ValidPair = getAdjacencyMatrix(side,1.5);
cluster2 = findClusters(Nuclei,2);
if norm(size(cluster2)-[24,2]) ~= 0
  pass = false;
  fprintf('\nfindClusters() --> failed at test 2.\n');
  return;
end
% 3
Nuclei.ValidPair = getAdjacencyMatrix(side,2);
cluster2 = findClusters(Nuclei,2);
if norm(size(cluster2)-[28,2]) ~= 0
  pass = false;
  fprintf('\nfindClusters() --> failed at test 3.\n');
  return;
end

% 3.1
Nuclei.ValidPair = getAdjacencyMatrix(side,inf);
cluster2 = findClusters(Nuclei,2);
if norm(size(cluster2)-[NchooseK(Nuclei.number,2),2]) ~= 0
  pass = false;
  fprintf('\nfindClusters() --> failed at test 3.1.\n');
  return;
end

%==========================================================================
%  2 x 2 x 2 cube 
%  3-clusters
%==========================================================================

% 4
Nuclei.ValidPair = getAdjacencyMatrix(side,1);
cluster3 = findClusters(Nuclei,3);
if norm(size(cluster3)-[24,3]) ~= 0
  cluster3
  size(cluster3)
  pass = false;
  fprintf('\nfindClusters() --> failed at test 4.\n');
  return;
end

% 5
Nuclei.ValidPair = getAdjacencyMatrix(side,1.5);
cluster3 = findClusters(Nuclei,3);
if norm(size(cluster3)-[56,3]) ~= 0
  cluster3
  size(cluster3)
  pass = false;
  fprintf('\nfindClusters() --> failed at test 5.\n');
  return;
end
% 6
Nuclei.ValidPair = getAdjacencyMatrix(side,2);
cluster3 = findClusters(Nuclei,3);
if norm(size(cluster3)-[56,3]) ~= 0
  size(cluster3)
  pass = false;
  fprintf('\nfindClusters() --> failed at test 6.\n');
  return;
end

%==========================================================================
%  2 x 2 x 2 cube 
%  4-clusters
%==========================================================================

% 6.1
Nuclei.ValidPair = getAdjacencyMatrix(side,2);
cluster4 = findClusters(Nuclei,4);
if norm(size(cluster4)-[70,4]) ~= 0
  size(cluster4)
  pass = false;
  fprintf('\nfindClusters() --> failed at test 6.1.\n');
  return;
end
% 6.2
Nuclei.ValidPair = getAdjacencyMatrix(side,1);
cluster4 = findClusters(Nuclei,4);
if norm(size(cluster4)-[38,4]) ~= 0
  cluster4
  size(cluster4)
  pass = false;
  fprintf('\nfindClusters() --> failed at test 6.2.\n');
  return;
end
%==========================================================================
%  3 x 3 x 3 cube 
%  2-clusters
%==========================================================================

% 7
side = 3;
Nuclei.number = side^3;  
Nuclei.ValidPair = getAdjacencyMatrix(side,1);
cluster2 = findClusters(Nuclei,2);
if norm(size(cluster2)-[54,2]) ~= 0
  cluster2
  size(cluster2)
  pass = false;
  fprintf('\nfindClusters() --> failed at test 7.\n');
  return;
end

%==========================================================================
%  3 x 3 x 3 cube 
%  3-clusters
%==========================================================================

% 8
Nuclei.ValidPair = getAdjacencyMatrix(side,1);
cluster3 = findClusters(Nuclei,3);
if norm(size(cluster3)-[171,3]) ~= 0
  cluster3
  size(cluster3)
  pass = false;
  fprintf('\nfindClusters() --> failed at test 4.\n');
  return;
end
%==========================================================================
%  10 x 10 x 10 cube 
%  4-clusters
%==========================================================================
  
% 9
% Takes about a minute to run if test 9 is included.
if true

  Nuclei.number = side^3;
  side = 10;
  Nuclei.ValidPair = getAdjacencyMatrix(side,1);
  cluster4 = findClusters(Nuclei,4);
  q = size(cluster4,1)/150000; % |clusters4| should be approximately 1000*6*5*5 = 150000, and is 47940, from code.
  log_q = log(q)/log(10);
  if  ~(log_q>0) && ~(log_q<1)
    q*150000
    pass = false;
    fprintf('\nfindClusters() --> failed at test 4.\n');
    return;
  end
end
fprintf('\nfindClusters() --> pass.\n');
end

function Nuclei = test_findClusters_generateNuc(side)
numNuc=side^3;
x = side;
y = x;
z = x;
for ii = numNuc:-1:1
  Nuclei.Index(ii) = ii;
  z = mod(z-1, side);
  if z==(side-1)
    y = mod(y-1, side);
    if y==(side-1)
      x = mod(x-1, side);
    end
  end
  Nuclei.Coordinates(ii,:) = [x,y,z];
end
end


function adj = getAdjacencyMatrix(side,threshold)
numNuc=side^3;

adj = eye(numNuc);
n = 0;
Coordinates = zeros(numNuc,3);
for ix = 1:side
  for jy = 1:side
    for kz = 1:side
      n = n+1;
      Coordinates(n,:) = [ix, jy, kz];      
    end
  end
end

for ii = 1:numNuc
  for jj = 1:ii-1
  r = norm(Coordinates(ii,:)-Coordinates(jj,:));
  if r<=threshold
    adj(ii,jj) = 1;
    adj(jj,ii) = 1;
  end
  end
end
end

