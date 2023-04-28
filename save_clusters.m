function save_clusters(filename,clusters)

max_size = numel(clusters);

fileID = fopen(filename,'w');

for cluster_size = 1:max_size
  fprintf(fileID,'n_%i = %i;\n',cluster_size,size(clusters{cluster_size},1) );
end

fprintf(fileID,'\n#[clusters]\n');

for cluster_size = 1:max_size

  n_clusters = size(clusters{cluster_size},1);

  str = repmat('%i, ',1,cluster_size);
  str = ['[', str(1:end-2), '];\n'];
  
  for icluster = 1:n_clusters
    fprintf(fileID,str,clusters{cluster_size}(icluster,:));
  end
end

fclose(fileID);
end
