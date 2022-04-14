function writeSpinPDB(Nuclei,NuclerContribution,outfilename,options)
  %{ 
    Write the nuclear spin contributions to a PDB file.  
    Store the contributions in the occupancy catagory.
  
    http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
        -------------------------------------------------------------------------------------
     1 -  6        Record name   "ATOM  "
     7 - 11        Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator.
    18 - 20        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues.
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    77 - 78        LString(2)    element      Element symbol, right-justified.
    79 - 80        LString(2)    charge       Charge  on the atom.
  %}
   if ~isfield(options,'Honly')
     options.Honly = false;
   end

   gmx_indices = cell(3,1);
   gmx_groups = cell(1,1);
   gmx_groups{1} = 'System';
   gmx_groups{2} = 'ROH';
   gmx_groups{3} = 'RCH';

   fileID = fopen(outfilename,'w');
   printElectronLine = false;
   if printElectronLine
       line = getElectronLine(Nuclei);
       fprintf(fileID,[line, '\n']);
   end
   N = Nuclei.number;
   counter = 0;
   for n = 1:N

     if options.Honly && ~strcmp(Nuclei.Element{n},'H')
       continue;
     end

     line = 'ATOM';
     line = padstring(line,6);
    
     % serial number
     counter = counter + 1;
     str = num2str(counter);
     gmx_indices{1} = [gmx_indices{1}; counter];


     str = padstring(str,11-6,'left');
     line = [line, str];
     line = padstring(line,11);
     
     str = Nuclei.Element{n};
%      str = [str,dec2hex(n)];
     str = [' ', padstring(str,16-13,'right')];
     line = [line, str];
     line = padstring(line,16);
     
     %line = [line, ' '];
     line = padstring(line,17);
     if Nuclei.Exchangeable(n)
       line = [line, 'ROH A '];
       gmx_indices{2} = [gmx_indices{2}; counter];
     else
       line = [line, 'RCH A '];
       gmx_indices{3} = [gmx_indices{3}; counter];
     end
     
     %line = [line, ' '];
     line = padstring(line,22);
     
     str = '0';
     str = padstring(str,26-23,'left');
     line = [line, str];
     line = padstring(line,26);
     
     line = padstring(line,30);
     for ix = 1:3
       x = Nuclei.PDBCoordinates(n,ix)*1e10;
       str = num2str(x);
       sgnStr = ' ';
       if x < 0
           sgnStr = '';
       end
       str = [sgnStr, str];

       if length(str)>7
           periodIdx = find(str=='.');
           decStr = str(periodIdx+1:end);
           str = [sgnStr, num2str( round(x,7 - length(decStr)  ))];
       end
       str = padstring(str,8);

       if abs(x-str2num(str))>1e-3
         error('Could not get nuclear coordinates.');
       end
       
       line = [line, str];
     end
     
     str = num2str(round( abs(NuclerContribution(n)),2),'%f');
     str = str(1:4);
     
     if abs(NuclerContribution(n)-str2num(str))>1e-3
       error('Could not get nuclear coordinates.');
     end
     line = padstring(line,54);
     
     str = padstring(str,6,'left');
     line = [line, str(1:6)];
     line = padstring(line,60);
     
     line = [line, '  0.00'];
     line = padstring(line,76);
     str = Nuclei.Element{n};
     str = padstring(str,78-76,'left');
     line = [line, str];
     
     
     line = padstring(line,80);
     fprintf(fileID,[line, '\n']);
   end
 
   line = 'END';
   fprintf(fileID,line);
   write_gmx_ndx();

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    function write_gmx_ndx()
    whiteSpaces = '';
    for ii = 1:numel(num2str(max(gmx_indices{1})))+1
        whiteSpaces = [whiteSpaces,' '];
    end

    ndxfilename = [outfilename(1:end-4),'.ndx'];
    gmxfileID = fopen(ndxfilename,'w');
    for ii=1:numel(gmx_groups)
        gmx_indices{ii} = sort(unique(gmx_indices{ii}));

        if ii ==1
            groupStr = ['[ ',gmx_groups{ii}, ' ]\n'];
        else
            groupStr = ['\n\n[ ',gmx_groups{ii}, ' ]\n'];
        end
        %disp(groupStr);
        fprintf(gmxfileID,groupStr);


        for iline =1:15:numel(gmx_indices{ii})
            line = '';
            for jj = 1:15
                if iline+jj-1 > numel(gmx_indices{ii})
                    break;
                end
                numStr = num2str(gmx_indices{ii}(iline+jj-1));
                numStr = [whiteSpaces(1:end - numel(numStr) ),...
                    numStr];

                line = [line, numStr];
            end
            line = [line,'\n'];
            %disp(line);
            fprintf(gmxfileID,line);
        end


    end

    return;
    end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
return;
end

function out = padstring(str,len,side)

if nargin<3, side = 'right'; end

nPadChars = len-length(str);

chr = ' ';
padding = repmat(chr,1,nPadChars);

switch side
  case 'right'
    out = [str padding];
  case 'left'
    out = [padding str];
  otherwise
    error('Unrecognized side specification - must be ''left'' or ''right''.');
end

end


function line = getElectronLine(Nuclei)

line = 'ATOM';
line = padstring(line,6);

str = '0';
str = padstring(str,11-6,'left');
line = [line, str];
line = padstring(line,11);

str = 'e';

str = padstring(str,16-11,'left');
line = [line, str];
line = padstring(line,16);

%line = [line, ' '];
line = padstring(line,17);

line = [line, 'RES'];

%line = [line, ' '];
line = padstring(line,22);

str = '0';
str = padstring(str,26-22,'left');
line = [line, str];
line = padstring(line,26);

line = padstring(line,30);
for ix = 1:3
  x = Nuclei.Electron_pdbCoordinates(ix)*1e10;
  if x<0
    str = ['-', num2str(floor(-x))];
  else
    str = num2str(floor(x));
  end
  
  if mod(abs(x),1) >= 1e-4
    str2 = num2str(mod(abs(x),1));
  else
    str2 = '0.000';
  end
  
  while length(str2) < 8
    str2 = [str2,'0'];
  end
  
  if str2(2)~='.'
    str2 = ['0.',str2];
  end
  
  str = [str,str2(2:end)];
  if x >= 0
    str = [' ',str];
  end
  
  if length(str) > 7
    str = str(1:7);
  end
  
  str = padstring(str,8);
  if abs(x-str2num(str))>1e-3
    error('Could not get nuclear coordinates.');
  end
  
  line = [line, str];
end

str = '1.00';

line = padstring(line,54);

str = padstring(str,6,'left');
line = [line, str(1:6)];
line = padstring(line,60);

line = [line, '  0.00'];
line = padstring(line,76);
str = 'e';
str = padstring(str,78-76,'left');
line = [line, str];


line = padstring(line,80);
end



