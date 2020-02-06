for clusterSize = 1:6
  opNames = getOpNames(clusterSize);
  
  fileID = fopen(['OpList_clusterSize_',num2str(clusterSize),'.txt'],'w');
  N = numel(opNames);
  
  for ii = 1:N
    nstr = [opNames{ii}, ' = ', num2str(ii), ';'];
    if ii==1 || (ii>2 && mod(ii,3)==1)
      nstr = [nstr, '\n'];
    else
      nstr = [nstr, ' '];
    end
    fprintf(fileID,nstr);
  end
  fclose(fileID);
end