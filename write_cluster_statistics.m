function r_max = write_cluster_statistics(clusters,Nuclei,filename)

radii = vecnorm(Nuclei.Coordinates')';

counts = zeros(numel(radii),numel(clusters));
headers = cell(1,numel(clusters)+1);
headers{1} = 'radius';
for clustersize = 1:numel(clusters)
  headers{1 + clustersize} = ['counts_',int2str(clustersize),'_clusters'];
  for nuc_idx = clusters{clustersize}(:)'
    counts(nuc_idx,clustersize) = counts(nuc_idx,clustersize) + 1;
  end
end

r_max = -1;
if numel(clusters) >= 2
  sele = counts(:,2) > 0;
  r_max = max(radii(sele));
end

T = array2table([radii,counts]);
T.Properties.VariableNames(1:numel(headers)) = headers;
writetable(T,[filename , '.csv']);
end