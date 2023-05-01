function analyze_methyls(Nuclei,file_name)

fileID = fopen(file_name,'w');
header = 'methyl_group,center_x,center_y,center_z,normal_x,normal_y,normal_z\n';
fprintf(fileID,header);

methyl_ids = unique(Nuclei.MethylID);

for id = methyl_ids
  if id>=0
    break;
  end

  r_carbon = Nuclei.Coordinates(Nuclei.MethylID==id,:);


  hydrogens = Nuclei.MethylID==-id;
  cluster = Nuclei.Index(hydrogens);
  methyl_name = ['methyl_', sprintf('%d_',cluster)];
  r_hydrogens = Nuclei.Coordinates(hydrogens,:);

  center = mean(r_hydrogens,1);
  center_str = sprintf('%d,',center);

  normal = center - r_carbon;
  normal = normal/norm(normal);
  normal_str = sprintf('%d,',normal);
 
  line = [methyl_name,',',center_str,normal_str];
  line = [line(1:end-1),'\n'];
  
  fprintf(fileID,line);

end

fclose(fileID);
end