function writeSpinPDB(Nuclei,NuclerContribution,outfilename)
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
   
   fileID = fopen(outfilename,'w');
   N = Nuclei.number;
   for n = 1:N
     line = 'ATOM';
     line = padstring(line,6);
     
     str = num2str(n);
     str = padstring(str,11-6,'left');
     line = [line, str];
     line = padstring(line,11);
     
     str = Nuclei.Element{n};
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
       x = Nuclei.PDBCoordinates(n,ix)*1e10;
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
     str = num2str(abs(NuclerContribution(n)),'%f');
     
     if abs(NuclerContribution(n)-str2num(str))>1e-3
       error('Could not get nuclear coordinates.');
     end
     str = padstring(str,6,'left');
     line = [line, str];
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