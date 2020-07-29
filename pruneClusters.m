function Clusters = pruneClusters(Nuclei,inClusters,prune)
prune  = setPruneDefaults(prune);

Coordinates.r = vecnorm(Nuclei.Coordinates,2,2);
numberClusters = Nuclei.numberClusters;
DistanceMatrix = Nuclei.DistanceMatrix;
Method_order = Method.order;

ClusterGeo = getClusterGeoStats(inClusters,Coordinates, DistanceMatrix,numberClusters,Method_order); 

ModulationDepth = Nuclei.Statistics{1}.Modulation_Depth_p; % = matrix(N);
Hyperfine = abs(Nuclei.Statistics{1}.Hyperfine_perpendicular); % = matrix(N,1);
Nuclear_Dipole = abs(Nuclei.Statistics{1}.Nuclear_Dipole_perpendicular); % = matrix(N);
Frequency_Pair = Nuclei.Statistics{1}.Frequency_Pair_p; % = matrix(N);
DeltaHyperfine = abs(Nuclei.Statistics{1}.DeltaHyperfine_perpendicular); % = matrix(N);
Adjacency = Nuclei.Adjacency; % = matrix(N);

% ENUM
MIN = 1; MAX = 2; EDGE = 3; CRIT = 4;
ClusterH = cell(1,Method_order);

isize = 1;
ClusterH{isize}.ENUM = ['SELF']; 
ClusterH{isize}.Hyperfine  = Hyperfine;

for isize = 2:Method_order
  
  nC_ = numberClusters(isize);
  ClusterH{isize} = getClusterHStats(Hyperfine,DeltaHyperfine,Nuclear_Dipole, ModulationDepth,Frequency_Pair,Adjacency,Coordinates,inClusters,ClusterGeo,numberClusters,isize,TM_powder);



end

Clusters = inClusters;
for isize = 2:Method_order
  
C = inClusters{isize};

nC_ = numberClusters(isize);
keep = true(nC_,1);  

keep = keep & (ClusterGeo{isize}.proximity >= prune.rmin);
keep = keep & (ClusterGeo{isize}.distance <= prune.rmax);
keep = keep & (ClusterH{isize}.Nuclear_Dipole(:,prune.criterion) >= prune.DDmin);
keep = keep & (ClusterH{isize}.Nuclear_Dipole(:,prune.criterion) <= prune.DDmax);

Clusters{isize} = C(keep,:);

end

end


function prune  = setPruneDefaults(prune)

if ~isfield(prune,'rmin')
  prune.rmin = 0;
end  
if ~isfield(prune,'rmax')
  prune.rmin = inf;
end  

if ~isfield(prune,'DDmin')
  prune.DDmin = 0;
end  
if ~isfield(prune,'DDmax')
  prune.DDmax = inf;
end 

% ENUM
MIN = 1; MAX = 2; EDGE = 3; CRIT = 4;

if ~isfield(prune,'criterion')
  prune.criterion = CRIT;
end
end


