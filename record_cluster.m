function record_cluster(cluster, OutputData)

filename = [OutputData(1:end-4),'_cluster_records.txt'];
fileID = fopen(filename,'a');

cluster(cluster<=0) = [];
clusterSize = length(cluster);

switch clusterSize
  case 1
    fprintf(fileID,'%d\n', cluster);
  case 2
    fprintf(fileID,'%d %d\n', cluster);
  case 3
    fprintf(fileID,'%d %d %d\n', cluster);
  case 4
    fprintf(fileID,'%d %d %d %d\n', cluster);
end
fclose(fileID);
end