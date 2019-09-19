function pass = isSubcluster(Subcluster,Cluster)

pass = true;
for ii=Subcluster
  if ~any(ii == Cluster)
    pass = false;
    return
  end
end
end