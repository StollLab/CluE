function Clusters = combineClusters(Clusters1,Clusters2)

clusterSize = min(numel(Clusters1),numel(Clusters2));

if numel(Clusters1) >= numel(Clusters2)
  Clusters = Clusters1;
else
  Clusters = Clusters2;
end
  
for isize = 1:clusterSize
  C_ = [ Clusters1{isize}; Clusters2{isize}];
  
  % Sort nuclei in each cluster by nuclei index
  C_ = sort(C_,2);
  
  % Sort clusters
  C_ = sortrows(C_);
  
  % Remove duplicates
  keep = [true; any(C_(1:end-1,:)~=C_(2:end,:),2)];
  
  Clusters{isize} = C_(keep,:);
end


end