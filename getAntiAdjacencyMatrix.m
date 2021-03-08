function AntiAdjacency = getAntiAdjacencyMatrix(System, Nuclei,Method)

% Initialize a matrix of valid edges.
AntiAdjacency = false(Nuclei.number,Nuclei.number,Method.extraOrder);

if ~System.Methyl.methylMethylCoupling
  isMethyl = strcmp(Nuclei.Type,'CH3');
  mm = isMethyl + isMethyl';
  mm = mm - diag(diag(mm)) > 1;
  
  for isize = 1:Method.extraOrder
    AntiAdjacency(:,:,isize) = AntiAdjacency(:,:,isize) | mm;
  end
end

end