function analyze_methyls(Nuclei,file_name)

fileID = fopen(file_name,'w');
header = ['methyl_group,','center_x,center_y,center_z,',...
    'normal_x,normal_y,normal_z,',...
    'plane_x,plane_y,plane_z\n'];
fprintf(fileID,header);

methyl_ids = unique(Nuclei.MethylID);

for id = methyl_ids
  if id>=0
    break;
  end

  r_carbon = Nuclei.Coordinates(Nuclei.MethylID==id,:);


  hydrogens = Nuclei.MethylID==-id;
  cluster = Nuclei.Index(hydrogens);
  assert(numel(cluster)==3);
  methyl_name = ['methyl_', sprintf('%d_',cluster)];
  methyl_name = methyl_name(1:end-1);
  
  r_hydrogens = Nuclei.Coordinates(hydrogens,:);

  center = mean(r_hydrogens,1);
  center_str = sprintf('%d,',center);

  axis = center - r_carbon;
  axis = axis/norm(axis);

  plane_dir = r_hydrogens(1,:) - center;
  plane_dir = plane_dir/norm(plane_dir);
  plane_dir_str = sprintf('%d,',plane_dir);
 
  normal = cross(r_hydrogens(2,:) - center,plane_dir);
  normal = normal/norm(normal);
  if normal'*axis < 0
    normal = -normal;
  end
  normal_str = sprintf('%d,',normal);
  line = [methyl_name,',',center_str,normal_str,plane_dir_str];
  line = [line(1:end-1),'\n'];
  
  fprintf(fileID,line);

end

fclose(fileID);
end
