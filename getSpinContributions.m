%function getSpinContributions(System,Nuclei,Signal, AuxiliarySignal,Clusters,order,TM)
function NuclearContribution = getSpinContributions(System, Nuclei, Signals, AuxiliarySignal, Clusters,Method,OutputData,gridWeight,gammaGridSize, TM)
  
  
  order = Method.order;
  % Loop through all the spins and calculate the contribution to the total TM. 
  
  % Determine how to find v(TM).
  TM_weights = getTMweights(2*System.Time,TM);
  
  % Find all clusters containing each nucleus.
  N = Nuclei.number;
  SuperClusterIndices = findSuperClusters(N,Clusters,order);
  
  % Find nucler spin contributions.  
  nSignal = length(Signals);
  nuclearContribution = zeros(N,nSignal);

  % loop through all orientations' signals
  signal_weights = zeros(nSignal,1);
  powder_contribution = zeros(N,1);
  for irot = 1:nSignal
    Signal = Signals{irot};
    signal_weights(irot) = abs(Signal(1));
    Signal = abs(Signal)/signal_weights(irot);

    ln_v_TM = log(sum(abs(Signal).*TM_weights));
    
    % loop through all nuclei.
    for inucleus =1:N
      if strcmp(Method.method,'rCCE')        
        xn_ = log(sum(TM_weights.*AuxiliarySignal{irot}{inucleus} ));
        nuclearContribution(inucleus,irot) = xn_/ln_v_TM;
      else
        nuclearContribution(inucleus,irot) = contributionMetric(inucleus,Signal,AuxiliarySignal,SuperClusterIndices,order,TM_weights);
      end
      
      powder_contribution(inucleus) = powder_contribution(inucleus) + signal_weights(irot)*nuclearContribution(inucleus,irot);
      spin_distance(inucleus) = norm(Nuclei.Coordinates(inucleus,:));
    end
    
  end
  
  % Save as PDB file.
  outfilename = [OutputData(1:end-4), '_NuclearSpinContribution.pdb'];
  
  total_contribution = sum(powder_contribution);
  fprintf('Total contribution = %d.\n',total_contribution);
  max_contribution = max(powder_contribution); 
  fprintf('Max contribution = %d.\n',max_contribution);
  min_contribution= min(powder_contribution);
  fprintf('Min contribution = %d.\n',min_contribution);
  if abs(1-total_contribution)>1e-6
    fprintf('Nuclear contribution do not sum to unity.  Normalizing.');
    fprintf('Normalizing contributions.\n');
    powder_contribution = powder_contribution/total_contribution;
    total_contribution = sum(powder_contribution);
    fprintf('Total contribution = %d.\n',total_contribution);
    max_contribution = max(powder_contribution);
    fprintf('Max contribution = %d.\n',max_contribution);
    min_contribution= min(powder_contribution);
    fprintf('Min contribution = %d.\n',min_contribution);
  end
  NuclearContribution.nuclearContribution = nuclearContribution;
  NuclearContribution.powder_contribution = powder_contribution;
  NuclearContribution.spin_distance=spin_distance;
  
  writeNuclerContributionPDB(Nuclei,powder_contribution,outfilename);
  save([OutputData(1:end-4), '_NuclearSpinContribution.mat'],'');
end

function xn = contributionMetric(n,Signal,AuxiliarySignal,SuperClusterIndices,order,TM_weights)
  % Calculate the nucler spin contribution xn, for nucleus n, according to
  % xn := 1/ln(v(TM)) sum_{C} ln(vn(TM))/|C|.
  ln_v_TM = log(sum(abs(Signal).*TM_weights));
  xn = 0;
  % Loop over all investigated cluster sizes > 1. 
  for clusterSize = 2:order
    
    % Loop over all clusters of size clusterSize containing nucleus n.
    for iCluster = SuperClusterIndices{n}{clusterSize}
      
      Nsig = size(AuxiliarySignal,2);
      
      ln_v_nC = zeros(size( AuxiliarySignal{1}{clusterSize,iCluster} ));
      
      for isignal = 1:Nsig
        %         index_gamma = mod(isignal-1,gammaGridSize) + 1;
        %         igrid = 1 + (isignal - index_gamma)/gammaGridSize;
        %         vnC = vnC + gridWeight(igrid,index_gamma)*abs(AuxiliarySignal{isignal}{clusterSize,iCluster}); % simpler variable name.
        %         normFactor = normFactor + gridWeight(igrid,index_gamma);
        ln_v_nC = ln_v_nC + log(abs(AuxiliarySignal{isignal}{clusterSize,iCluster})); % simpler variable name.
        
      end
      ln_vn_at_TM = sum(TM_weights.*ln_v_nC); % interpolating vnC at the time TM.
      if isnan(ln_vn_at_TM)
        error('ln_vn_at_TM is NaN.');
      end
      % Apply metric definition.
      xn_ = ln_vn_at_TM/clusterSize/ln_v_TM;
      xn = xn + xn_;
      
    end
  end
 
  
 if xn < 0
   fprintf('Nuclear spin contribution is  %d < 0;\n',xn);
   fprintf('r0 may be to small.\n');
 end
   

end

function TM_weights = getTMweights(Time,TM)
  % Currently only a two-point interpolation is employed, but the 
  % implementation should be flexible enough to allow for a higher-order 
  % interpolation without changing the code outside this function. 

  TM_weights = zeros(size(Time));
  sup_index = sum(Time<=TM); % highest index with t <= TM.
  if sup_index==length(Time)
    warning("The TM >= max time.");
    TM_weights(sup_index) = 1;
    return;
  end
  Dt = Time(sup_index+1)-Time(sup_index); % smallest time interval.
  dt = TM -Time(sup_index); % time between TM and the previous time point.
  % Define weight to place v(TM) as the point t = TM, on the line conecting
  % the signal values to either side of t = TM.
  TM_weights(sup_index) = 1-dt/Dt;
  TM_weights(sup_index+1) = dt/Dt;
end


function Indices = findSuperClusters(Nuclei_number,Clusters,clusterSize)
  % Indices{n}{size} = list of all jCluster such that
  % Clusters{size}(jCluster,:) contains nucleus n;
  
  % Clusters{clusterSize} =  number of clusterSize-clusters by clusterSize 
  % matrix with each row being a list of nuclear indices.
  Indices{Nuclei_number}{clusterSize}=[];  
  for n=1:Nuclei_number
    for isize = 2:clusterSize
      containingClusters = sum(Clusters{isize}==n,2);
      Indices{n}{isize} = find(containingClusters)';
    end
  end
end

function writeNuclerContributionPDB(Nuclei,NuclerContribution,outfilename)
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
