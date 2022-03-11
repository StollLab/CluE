

function make_ndx(pdbFileName,sysOnly)

pdb = parsePDBfile(pdbFileName);

groups = cell(1,1);
indices = cell(1,1);

groups{1}  = 'System';
indices{1} = [];


for ipdb =1:pdb.number
    grpIdx = 1;
    indices{grpIdx} = [indices{grpIdx}; pdb.serial(ipdb)];
    if ~sysOnly
        for igroup = 2:numel(groups)
            if strcmp(groups{igroup}, pdb.resName{ipdb})
                grpIdx = igroup;
                break;
            end
        end

        if grpIdx==1
            grpIdx = numel(groups)+1;
            groups{grpIdx}  = { pdb.resName{ipdb}};
            indices{grpIdx} = {[]};
        end
        indices{grpIdx} = [indices{grpIdx}; pdb.serial(ipdb)];
    end
end

outfilename = [pdbFileName(1:end-4),'.ndx'];
fileID = fopen(outfilename,'w');
for ii=1:numel(groups)
    indices{ii} = sort(unique(indices{ii}));

    groupStr = ['\n[ ',groups{ii}, ' ]\n'];

    disp(groupStr);
    fprintf(fileID,groupStr);


    for iline =1:10:numel(indices{ii})
        line = '';
        for jj = 1:10
            if iline+jj-1 > numel(indices{ii})
                break;
            end
            line = [line, num2str(indices{ii}(iline+jj-1)),' 	'];
        end
        line = [line,'\n'];
        disp(line);
        fprintf(fileID,line);
    end


end


end

%{

function make_ndx(particles,uniqueResidueList,pdbID,pdb_number)
  
groups = cell(numel(uniqueResidueList) + 1,1);
indices = cell(numel(uniqueResidueList) + 1,1);

groups{1}  = 'System';
indices{1} = [];

for itype = 1:numel(particle)
  indices{1} = [indices{1}; particles{itype}.members];
end

for ii =1:numel(uniqueResidueList)
  groups{ii+1} = uniqueResidueList{ii};
  indices{ii+1} = [];

  for itype = 1:numel(particle)
    if strcmp(uniqueResidueList{ii},particles{itype}.resName)
        
        newIndices = pdbID(particles{itype}.members);
        newIndices(newIndices>pdb_number) = []; 

        indices{ii+1} = [indices{ii+1}; newIndices];
    end
  end
end

fileID = fopen(outfilename,'w');
for ii=1:numel(groups)
  indices{ii} = sort(unique(indices{ii}));

  groupStr = ['\n[ ',groups{ii}, ' ]\n'];

  fprintf(fileID,groupStr);
  
  line = '';
  for iline =1:10:numel(indices{ii})
    for jj = 1:10
      if iline+jj-1 > numel(indices{ii})
        break;
      end
      line = [line, num2str(indices{ii}(iline+jj-1)),' 	'];
    end
    line = [line,'\n'];
    fprintf(fileID,line);
  end


end


end


%}