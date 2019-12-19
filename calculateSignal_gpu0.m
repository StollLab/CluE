% Clusters = Clusters(cluster index , 1:size ,order)
% Clusters(cluster index , size > order ,order) = 0.

function [Signal, AuxiliarySignal_1,AuxiliarySignal_2,AuxiliarySignal_3,AuxiliarySignal_4,Signals] ... 
       = calculateSignal_gpu0(System, Method, Nuclei,Clusters, timepoints,dt,t0)

maxSize = 6;     
     
% ENUM
FID = 1; HAHN = 2; CPMG = 3; CPMG_CONST = 4; CPMG_2D = 5;

Method_order = Method.order;  
System_full_Sz_Hyperfine = System.full_Sz_Hyperfine;
maxClusterSize = min(maxSize,Method_order);
% Convert variable to gpu compatible forms
dimensionality = 1;
if strcmp(System.experiment,'FID')
  EXPERIMENT = FID;
  total_time = System.Time(end);
elseif strcmp(System.experiment,'Hahn')
  EXPERIMENT = HAHN;
  total_time = 2*System.Time(end);
elseif strcmp(System.experiment,'CPMG')
  EXPERIMENT = CPMG;
  total_time = 4*System.Time(end);
elseif strcmp(System.experiment,'CPMG-const')
  EXPERIMENT = CPMG_CONST;
  total_time = 4*System.Time(end);
elseif strcmp(System.experiment,'CPMG-2D')
  EXPERIMENT = CPMG_2D;
  total_time = 4*System.Time(end);
  dimensionality = 2;
else
  error('The experiment ''%s'' is not supported.',System.experiment);
end
Nuclei_Coordinates = Nuclei.Coordinates;
Nuclei_Abundance = Nuclei.Abundance;
Nuclei_Spin = Nuclei.Spin;
Nuclei_g = Nuclei.Nuclear_g;
NumberStates = Nuclei.NumberStates;
state_multiplicity = Nuclei.StateMultiplicity;
ZeemanStates = Nuclei.ZeemanStates;
max_basis = max(NumberStates);
States = zeros(max_basis,Nuclei.number);
 
Qtensors = Nuclei.Qtensor;

% Unpackspin operators.
Op = Nuclei.SpinOperators;
SpinXiXjOps = Nuclei.SpinXiXjOperators;

% Set 1-cluster operators.
if maxClusterSize > 0
  Spin2Op1 = Op{2}{1};
  Spin3Op1 = Op{3}{1};
  Spin4Op1 = Op{4}{1};
  SpinXiXjOp_1 = SpinXiXjOps{1};
else
  Spin2Op1 = [];
  Spin3Op1 = [];
  Spin4Op1 = [];
  SpinXiXjOp_1 = [];
end

% Set 2-cluster operators.
if maxClusterSize > 1
  Spin2Op2 = Op{2}{2};
  Spin3Op2 = Op{3}{2};
  Spin4Op2 = Op{4}{2};
  SpinXiXjOp_2 = SpinXiXjOps{2};
else
  Spin2Op2 = [];
  Spin3Op2 = [];
  Spin4Op2 = [];
  SpinXiXjOp_2 =[];
end

% Set 3-cluster operators.
if maxClusterSize > 2
  Spin2Op3 = Op{2}{3};
  Spin3Op3 = Op{3}{3};
  Spin4Op3 = Op{4}{3};
  SpinXiXjOp_3 = SpinXiXjOps{3};
else
  Spin2Op3 = [];
  Spin3Op3 = [];
  Spin4Op3 = [];
  SpinXiXjOp_3 =[];
end

% Set 4-cluster operators.
if maxClusterSize > 3
  Spin3Op4 = Op{3}{4};
  Spin2Op4 = Op{2}{4};
  Spin4Op4 = Op{4}{4};
  SpinXiXjOp_4 = SpinXiXjOps{4};
else
  Spin3Op4 = [];
  Spin2Op4 = [];
  Spin4Op4 = [];
  SpinXiXjOp_4 =[];
end

% Set 5-cluster operators.
if maxClusterSize > 4
  Spin2Op5 = Op{2}{5};
  Spin3Op5 = Op{3}{5};
  Spin4Op5 = Op{4}{5};
  SpinXiXjOp_5 = SpinXiXjOps{5};
else
  Spin2Op5 = [];
  Spin3Op5 = [];
  Spin4Op5 = [];
  SpinXiXjOp_5 =[];
end

% Set 6-cluster operators.
if maxClusterSize > 5
  Spin2Op6 = Op{2}{6};
  Spin3Op6 = Op{3}{6};
  Spin4Op6 = Op{4}{6};
  SpinXiXjOp_6 = SpinXiXjOps{6};
else
  Spin2Op6 = [];
  Spin3Op6 = [];
  Spin4Op6 = [];
  SpinXiXjOp_6 =[];
end

numberClusters = Nuclei.numberClusters(1:maxClusterSize);
maxNumberClusters = max(numberClusters(1:maxClusterSize));


% CluserArray(iCluster,:,clusterSize) = nuclear indices.
ClusterArray = zeros(maxNumberClusters,maxClusterSize,maxClusterSize);
 

for isize = 1:Method_order
  
  for ii = 1:Nuclei.numberClusters(isize)
    ClusterArray(ii,1:isize,isize) = Clusters{isize}(ii,:);
  end
  
      switch isize 
        
        % SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
        % the jth cluster of size subCluster_size that is a subcluster of
        % the ith ccluster of size clusterSize.
        
        case 1
          Coherences_1 = ones(numberClusters(isize),timepoints^dimensionality);  
        case 2
          SubclusterIndices_2 = zeros(nchoosek(isize,1), isize , Nuclei.numberClusters(isize)); 
          Coherences_2 = ones(numberClusters(isize),timepoints^dimensionality);
          
        case 3
          SubclusterIndices_3 = zeros(nchoosek(isize,1), isize , Nuclei.numberClusters(isize));
          Coherences_3 = ones(numberClusters(isize),timepoints^dimensionality);
        case 4
          SubclusterIndices_4 = zeros(nchoosek(isize,2), isize , Nuclei.numberClusters(isize));
          Coherences_4 = ones(numberClusters(isize),timepoints^dimensionality);
        case 5
          SubclusterIndices_5 = zeros(nchoosek(isize,3), isize , Nuclei.numberClusters(isize));
          Coherences_5 = ones(numberClusters(isize),timepoints^dimensionality);
        case 6
          SubclusterIndices_6 = zeros(nchoosek(isize,3), isize , Nuclei.numberClusters(isize));
          Coherences_6 = ones(numberClusters(isize),timepoints^dimensionality);
      end
      
end

% Define placeholder variables.
for isize = Method_order+1:maxSize
  
      switch isize 
        
        case 1
          Coherences_1 = 1;  
        case 2
          SubclusterIndices_2 = []; 
          Coherences_2 = 1;
          
        case 3
          SubclusterIndices_3 = [];
          Coherences_3 = 1;
        case 4
          SubclusterIndices_4 = 1;
          Coherences_4 = 1;
          
        case 5
          SubclusterIndices_5 = [];
          Coherences_5 = 1;
        case 6
          SubclusterIndices_6 = 1;
          Coherences_6 = 1;
      end
      
end
ge=System.gMatrix(3,3);
magneticField = System.magneticField;
muB = System.muB;
muN = System.muN;
mu0 = System.mu0;
hbar = System.hbar;

Nuclei_ValidPair = Nuclei.ValidPair;
% ENUM
CONNECTED = 0;  COMPLETE = 1;

graphCriterion = CONNECTED;
if strcmp(Nuclei.graphCriterion,'complete')
  graphCriterion = COMPLETE;
end

useHamiltonian = System.useHamiltonian;

% Methyl Groups
IsMethyl = strcmp(Nuclei.Type,'CH3');
if ~isempty(Nuclei.Methyl_Data)
  Methyl_P3 = Nuclei.Methyl_Data.Projection_3;
  Methyl_P4 = Nuclei.Methyl_Data.Projection_4;
  Methyl_P5 = Nuclei.Methyl_Data.Projection_5;
  Methyl_P6 = Nuclei.Methyl_Data.Projection_6;
end
Methyl_gA = [0,0];
%--------------------------------------------------------------------------
% gpu code
%--------------------------------------------------------------------------


for clusterSize = 1:Method_order
   
  % Find coherences
  numClusters = numberClusters(clusterSize);
  
  for iCluster = numClusters:-1:1
    Cluster = ClusterArray(iCluster,1:clusterSize,clusterSize); 
   
    if Cluster(1,1) == 0 || ~validateCluster(Cluster,Nuclei_ValidPair,graphCriterion)
      continue;
    end
    
    thisClusterSize = clusterSize + 2*sum(IsMethyl(Cluster));
    thisCluster = zeros(1,thisClusterSize);
    thisIndex = 0;
    MethylID = zeros(1,thisClusterSize);
    methyl_number = 0;
    for ii = 1:clusterSize

      thisIndex = thisIndex + 1;
      
      if IsMethyl(Cluster(ii))
        methyl_number = methyl_number + 1;
        Methyl_gA(methyl_number) = Nuclei_Abundance(Cluster(ii));
        MethylID(thisIndex:thisIndex+2) = methyl_number;
        thisCluster(thisIndex) = Cluster(ii) + 1;
        thisIndex = thisIndex + 1;
        thisCluster(thisIndex) = Cluster(ii) + 2;
        thisIndex = thisIndex + 1;
        thisCluster(thisIndex) = Cluster(ii) + 3;
      else
        thisCluster(thisIndex) = Cluster(ii);
      end
      
    end
    
    switch clusterSize 
      % SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
      % the jth cluster of size subCluster_size that is a subcluster of
      % the ith ccluster of size clusterSize.
      case 2
        SubclusterIndices_2(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);   
      case 3
        SubclusterIndices_3(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 4
        SubclusterIndices_4(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 5
        SubclusterIndices_5(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      case 6
        SubclusterIndices_6(:,:,iCluster) = findSubclusters_gpu(ClusterArray,clusterSize,iCluster,clusterSize);
      otherwise
        continue;
        
    end
    
    % Select the appropriate spin operator.
     switch thisClusterSize
      case 1
        switch Nuclei_Spin(Cluster(1))
          case 1/2
            SpinOp = Spin2Op1;
            SpinXiXjOp = [];
          case 1
            SpinOp = Spin3Op1;
            SpinXiXjOp = SpinXiXjOp_1;
          case 3/2
            SpinOp = Spin4Op1;
            SpinXiXjOp = [];
        end
        
       case 2
         switch Nuclei_Spin(Cluster(1))
           case 1/2
             SpinOp = Spin2Op2;
             SpinXiXjOp = [];
           case 1
             SpinOp = Spin3Op2;
             SpinXiXjOp = SpinXiXjOp_2;
           case 3/2
             SpinOp = Spin4Op2;
             SpinXiXjOp = [];
         end
         
       case 3
         switch Nuclei_Spin(Cluster(1))
           case 1/2
             SpinOp = Spin2Op3;
             SpinXiXjOp = [];
           case 1
             SpinOp = Spin3Op3;
             SpinXiXjOp = SpinXiXjOp_3;
           case 3/2
             SpinOp = Spin4Op3;
             SpinXiXjOp = [];
         end
         
       case 4
         switch Nuclei_Spin(Cluster(1))
           case 1/2
             SpinOp = Spin2Op4;
             SpinXiXjOp = [];
           case 1
             SpinOp = Spin3Op4;
             SpinXiXjOp = SpinXiXjOp_4;
           case 3/2
             SpinOp = Spin4Op4;
             SpinXiXjOp = [];
         end
         
       case 5
         switch Nuclei_Spin(Cluster(1))
           case 1/2
             SpinOp = Spin2Op5;
             SpinXiXjOp = [];
           case 1
             SpinOp = Spin3Op5;
             SpinXiXjOp = SpinXiXjOp_5;
           case 3/2
             SpinOp = Spin4Op5;
             SpinXiXjOp = [];
         end
         
       case 6
         switch Nuclei_Spin(Cluster(1))
           case 1/2
             SpinOp = Spin2Op6;
             SpinXiXjOp = [];
           case 1
             SpinOp = Spin3Op6;
             SpinXiXjOp = SpinXiXjOp_6;
           case 3/2
             SpinOp = Spin4Op6;
             SpinXiXjOp = [];
         end
         
       otherwise
        continue;
         
     end
     
    [tensors,zeroIndex] = pairwisetensors_gpu(Nuclei_g, Nuclei_Coordinates,thisCluster,magneticField, ge, muB, muN, mu0, hbar,useHamiltonian,MethylID);
    
    [Halpha,Hbeta] = ...
      assembleHamiltonian_gpu(tensors,SpinOp,SpinXiXjOp,thisCluster,...
      System_full_Sz_Hyperfine,zeroIndex,thisClusterSize,...
      methyl_number,Qtensors,state_multiplicity);
  
    % get density matrix
    DensityMatrix0 = getDensityMatrix(ZeemanStates,NumberStates,Cluster);
    
    switch methyl_number
      case 0
        PA = eye(size(Halpha));
      case 1
        gA = Methyl_gA(1);
        gE = 1-gA;
        methyl_index = find(MethylID==1,1);
        switch thisClusterSize
          case 3
            PA = Methyl_P3(:,:,1);
            PEap = Methyl_P3(:,:,2);
            PEam = Methyl_P3(:,:,3);
            PEbp = Methyl_P3(:,:,4);
            PEbm = Methyl_P3(:,:,5);
           
          case 4      
            switch methyl_index
              case 1
                jmethyl = 6;
              case 2
                jmethyl = 1;
            end
                PA = Methyl_P4(:,:,jmethyl);
                PEap = Methyl_P4(:,:,jmethyl+1);
                PEam = Methyl_P4(:,:,jmethyl+2);
                PEbp = Methyl_P4(:,:,jmethyl+3);
                PEbm = Methyl_P4(:,:,jmethyl+4);
            
          case 5
            switch methyl_index
              case 1
                jmethyl = 11;
              case 2
                jmethyl = 6;
             case 3
                jmethyl = 1;
            end
                PA   = Methyl_P5(:,:,jmethyl);
                PEap = Methyl_P5(:,:,jmethyl+1);
                PEam = Methyl_P5(:,:,jmethyl+2);
                PEbp = Methyl_P5(:,:,jmethyl+3);
                PEbm = Methyl_P5(:,:,jmethyl+4);
            
          case 6
            switch methyl_index
              case 1
                jmethyl = 16;
              case 2
                jmethyl = 11;
             case 3
                jmethyl = 6;
             case 4
                jmethyl = 1;
            end
                PA   = Methyl_P6(:,:,jmethyl);
                PEap = Methyl_P6(:,:,jmethyl+1);
                PEam = Methyl_P6(:,:,jmethyl+2);
                PEbp = Methyl_P6(:,:,jmethyl+3);
                PEbm = Methyl_P6(:,:,jmethyl+4);
            
        end
      case 2
        gA = Methyl_gA(1)*Methyl_gA(2);
        gE = (1-Methyl_gA(1))*(1-Methyl_gA(2));
        
        jmethyl = 21;
        PA1   = Methyl_P6(:,:,jmethyl);
        PEap1 = Methyl_P6(:,:,jmethyl+1);
        PEam1 = Methyl_P6(:,:,jmethyl+2);
        PEbp1 = Methyl_P6(:,:,jmethyl+3);
        PEbm1 = Methyl_P6(:,:,jmethyl+4);
        
        PA2   = Methyl_P6(:,:,jmethyl+5);
        PEap2 = Methyl_P6(:,:,jmethyl+6);
        PEam2 = Methyl_P6(:,:,jmethyl+7);
        PEbp2 = Methyl_P6(:,:,jmethyl+8);
        PEbm2 = Methyl_P6(:,:,jmethyl+9);
        
        PA = PA1*PA2; 
        PEap = PEap1*PEap2;
        PEam = PEam1*PEam2;
        PEbp = PEbp1*PEbp2;
        PEbm = PEbm1*PEbm2;
        
    end
    
      Hb=PA*Hbeta*PA;
      Ha=PA*Halpha*PA;
      DensityMatrix=PA*DensityMatrix0*PA;

    % get cluster coherence
    
    switch clusterSize
      case 1
        Coherences_1(iCluster,:) = propagate(total_time, DensityMatrix,timepoints,dt,t0,Hb,Ha,EXPERIMENT); 
      case 2
        Coherences_2(iCluster,:) = propagate(total_time, DensityMatrix,timepoints,dt,t0,Hb,Ha,EXPERIMENT); 
      case 3
        Coherences_3(iCluster,:) = propagate(total_time, DensityMatrix,timepoints,dt,t0,Hb,Ha,EXPERIMENT); 
      case 4
        Coherences_4(iCluster,:) = propagate(total_time, DensityMatrix,timepoints,dt,t0,Hb,Ha,EXPERIMENT); 
      case 5
        Coherences_5(iCluster,:) = propagate(total_time, DensityMatrix,timepoints,dt,t0,Hb,Ha,EXPERIMENT); 
      case 6
        Coherences_6(iCluster,:) = propagate(total_time, DensityMatrix,timepoints,dt,t0,Hb,Ha,EXPERIMENT); 

    end 
    
    if methyl_number==0
      continue;
    end
    
    Hb_E =   PEap*Hbeta*PEap + PEam*Hbeta*PEam ...
           + PEbp*Hbeta*PEbp + PEbm*Hbeta*PEbm ...
           + PEap*Hbeta*PEbm + PEam*Hbeta*PEbp ...
           + PEbp*Hbeta*PEam + PEbm*Hbeta*PEap;
    
    Ha_E =   PEap*Halpha*PEap + PEam*Halpha*PEam ...
           + PEbp*Halpha*PEbp + PEbm*Halpha*PEbm ...
           + PEap*Halpha*PEbm + PEam*Halpha*PEbp ...
           + PEbp*Halpha*PEam + PEbm*Halpha*PEap;
         
    DensityMatrix_E =   PEap*DensityMatrix0*PEap + PEam*DensityMatrix0*PEam ...
           + PEbp*DensityMatrix0*PEbp + PEbm*DensityMatrix0*PEbm ...
           + PEap*DensityMatrix0*PEbm + PEam*DensityMatrix0*PEbp ...
           + PEbp*DensityMatrix0*PEam + PEbm*DensityMatrix0*PEap;     
         
    %{
    % Get dimensionality of the Hilbert space.
    dimension = size(DensityMatrix,1);
    StateNumbers = 1:dimension;
    AStates = [];
    
    for imethyl = 1:methyl_number
      state_cycles = 2.^find(MethylID==imethyl);
      StateList = zeros(3,dimension);
      
      for iH = 1:3
        n_=state_cycles(iH);
        for ii = 1:dimension/n_
%           StateList(iH, 1+ii*n_ :(ii+1)*n_/2)  = 0;
          StateList(iH, (ii-1/2)*n_ + 1:ii*n_)  = 1;
        end
        
      end
      AStates = [AStates ,find((StateList(1,:)==StateList(2,:)) .* (StateList(1,:)==StateList(3,:)))];
    end
    AStates = unique(AStates);
    keepStates= setdiff(StateNumbers,AStates);
    DensityMatrix_E = DensityMatrix(keepStates,keepStates);
    DensityMatrix_E =DensityMatrix_E/trace(DensityMatrix_E);
    Hb_E = Hb(keepStates,keepStates);
    Ha_E = Ha(keepStates,keepStates);
    %}
    
    Coherences_E = propagate(total_time, DensityMatrix_E,timepoints,dt,t0,Hb_E,Ha_E,EXPERIMENT);
    
    switch clusterSize
      case 1
        Coherences_1(iCluster,:) = gA*Coherences_1(iCluster,:) + gE*Coherences_E;
      case 2
        Coherences_2(iCluster,:) = gA*Coherences_2(iCluster,:) + gE*Coherences_E;
      case 3
        Coherences_3(iCluster,:) = gA*Coherences_3(iCluster,:) + gE*Coherences_E;
      case 4
        Coherences_4(iCluster,:) = gA*Coherences_4(iCluster,:) + gE*Coherences_E;
      case 5
        Coherences_5(iCluster,:) = gA*Coherences_5(iCluster,:) + gE*Coherences_E;
      case 6
        Coherences_6(iCluster,:) = gA*Coherences_6(iCluster,:) + gE*Coherences_E;
    end
    
     if methyl_number < 2
      continue;
     end
     
     Hb_AE =   PEap2*Hbeta*PEap2 + PEam*Hbeta*PEam2 ...
       + PEbp2*Hbeta*PEbp2 + PEbm2*Hbeta*PEbm2 ...
       + PEap2*Hbeta*PEbm2 + PEam2*Hbeta*PEbp2 ...
       + PEbp2*Hbeta*PEam2 + PEbm2*Hbeta*PEap2;
     Hb_AE = PA1*Hb_AE*PA1;
     
     Ha_AE =   PEap2*Halpha*PEap + PEam2*Halpha*PEam2 ...
       + PEbp2*Halpha*PEbp + PEbm2*Halpha*PEbm2 ...
       + PEap2*Halpha*PEbm + PEam2*Halpha*PEbp2 ...
       + PEbp2*Halpha*PEam + PEbm2*Halpha*PEap2;
     Ha_AE = PA1*Ha_AE*PA1;
     
     DensityMatrix_AE =   PEap2*DensityMatrix0*PEap2 + PEam2*DensityMatrix0*PEam2 ...
       + PEbp2*DensityMatrix0*PEbp2 + PEbm2*DensityMatrix0*PEbm2 ...
       + PEap2*DensityMatrix0*PEbm2 + PEam2*DensityMatrix0*PEbp2 ...
       + PEbp2*DensityMatrix0*PEam2 + PEbm2*DensityMatrix0*PEap2;
     DensityMatrix_AE = PA1*DensityMatrix_AE*PA1;
     
     Coherences_AE = propagate(total_time, DensityMatrix_AE,timepoints,dt,t0,Hb_AE,Ha_AE,EXPERIMENT);
    
     
     Hb_EA =   PEap1*Hbeta*PEap1 + PEam*Hbeta*PEam1 ...
       + PEbp1*Hbeta*PEbp1 + PEbm1*Hbeta*PEbm1 ...
       + PEap1*Hbeta*PEbm1 + PEam1*Hbeta*PEbp1 ...
       + PEbp1*Hbeta*PEam1 + PEbm1*Hbeta*PEap1;
     Hb_EA = PA2*Hb_EA*PA2;
     
     Ha_EA =   PEap1*Halpha*PEap + PEam1*Halpha*PEam1 ...
       + PEbp1*Halpha*PEbp + PEbm1*Halpha*PEbm1 ...
       + PEap1*Halpha*PEbm + PEam1*Halpha*PEbp1 ...
       + PEbp1*Halpha*PEam + PEbm1*Halpha*PEap1;
     Ha_EA = PA2*Ha_EA*PA2;
     
     DensityMatrix_EA =   PEap1*DensityMatrix0*PEap1 + PEam1*DensityMatrix0*PEam1 ...
       + PEbp1*DensityMatrix0*PEbp1 + PEbm1*DensityMatrix0*PEbm1 ...
       + PEap1*DensityMatrix0*PEbm1 + PEam1*DensityMatrix0*PEbp1 ...
       + PEbp1*DensityMatrix0*PEam1 + PEbm1*DensityMatrix0*PEap1;
     DensityMatrix_EA = PA2*DensityMatrix_EA*PA2;
     
     Coherences_EA = propagate(total_time, DensityMatrix_EA,timepoints,dt,t0,Hb_EA,Ha_EA,EXPERIMENT);

     switch clusterSize
       case 2
         Coherences_2(iCluster,:) = Coherences_2(iCluster,:)  ...
           + Methyl_gA(1)*(1-Methyl_gA(2))*Coherences_AE ...
           +(1 - Methyl_gA(1))*Methyl_gA(2)*Coherences_EA;
       case 3
         Coherences_3(iCluster,:) = Coherences_3(iCluster,:) ...
           + Methyl_gA(1)*(1-Methyl_gA(2))*Coherences_AE ...
           +(1 - Methyl_gA(1))*Methyl_gA(2)*Coherences_EA;
       case 4
         Coherences_4(iCluster,:) = Coherences_4(iCluster,:) ...
           + Methyl_gA(1)*(1-Methyl_gA(2))*Coherences_AE ...
           +(1 - Methyl_gA(1))*Methyl_gA(2)*Coherences_EA;
       case 5
         Coherences_5(iCluster,:) = Coherences_5(iCluster,:) ...
           + Methyl_gA(1)*(1-Methyl_gA(2))*Coherences_AE ...
           +(1 - Methyl_gA(1))*Methyl_gA(2)*Coherences_EA;
       case 6
         Coherences_6(iCluster,:) = Coherences_6(iCluster,:)  ...
           + Methyl_gA(1)*(1-Methyl_gA(2))*Coherences_AE ...
           +(1 - Methyl_gA(1))*Methyl_gA(2)*Coherences_EA;
     end
  end
  
  
end


% Calculate signal
%-------------------------------------------------------------------------------
[Signals, AuxiliarySignal_1,AuxiliarySignal_2,AuxiliarySignal_3,AuxiliarySignal_4] ...
  = doClusterCorrelationExpansion_gpu(Coherences_1,Coherences_2,Coherences_3,Coherences_4,Coherences_5,Coherences_6,ClusterArray, ...
  SubclusterIndices_2,SubclusterIndices_3,SubclusterIndices_4,SubclusterIndices_5,SubclusterIndices_6,...
  timepoints,dimensionality, Method_order,numberClusters, Nuclei_Abundance);
% 
% if EXPERIMENT == CPMG_2D
%     Signal = reshape(Signals(Method_order,:).',timepoints,timepoints)';
% else
%     Signal = Signals(Method_order,:);
% end

    Signal = Signals(Method_order,:);
end

% ========================================================================
% Density Matrix Function
% ========================================================================

function DensityMatrix = getDensityMatrix(States,NumberStates, Cluster)

dimension = prod(NumberStates(Cluster));

DensityMatrix = eye(dimension);

switchIndex = dimension;
for inucleus = Cluster
  switchIndex = switchIndex/NumberStates(inucleus);
  counter = 0;
  stateNumber = 1;
  for jj = 1:dimension
    if counter == switchIndex
      stateNumber = mod(stateNumber,NumberStates(inucleus)) + 1;
    end
    counter = mod(counter,switchIndex) + 1;
    DensityMatrix(jj,jj) = DensityMatrix(jj,jj)*States(inucleus,stateNumber)^2;
  end
end
%DensityMatrix = DensityMatrix.^2;
if abs(trace(DensityMatrix)-1)>1e-9
  error('State is not normalized.');
end

end

% ========================================================================
% Propagate Function
% ========================================================================
function Signal = propagate(total_time,DensityMatrix,timepoints,dt,t0,Hamiltonian_beta,Hamiltonian_alpha,EXPERIMENT)

% ENUM
FID = 1; HAHN = 2; CPMG = 3; CPMG_CONST = 4; CPMG_2D = 5;

DensityMatrix = eye(size(DensityMatrix)).*DensityMatrix/trace(DensityMatrix);
vecDensityMatrixT = reshape(DensityMatrix.',1,[]);

Hamiltonian_beta =(Hamiltonian_beta+Hamiltonian_beta')/2; 
Hamiltonian_alpha =(Hamiltonian_alpha+Hamiltonian_alpha')/2;

dU_beta = propagator_eig(Hamiltonian_beta,dt);
dU_alpha = propagator_eig(Hamiltonian_alpha,dt);

nStates = length(Hamiltonian_beta);
if t0 > 0
  U_beta = propagator_eig(Hamiltonian_beta,t0);
  U_alpha = propagator_eig(Hamiltonian_alpha,t0);
else
  U_beta = eye(nStates);
  U_alpha = eye(nStates);
end

if EXPERIMENT==CPMG
  if t0 > 0
    U_beta_2 = propagator_eig(Hamiltonian_beta,t0);
    U_alpha_2 = propagator_eig(Hamiltonian_alpha,t0);
  else
    U_beta_2 = eye(nStates);
    U_alpha_2 = eye(nStates);
  end
end

if EXPERIMENT == CPMG_CONST
  U_beta_2 = propagator_eig(Hamiltonian_beta,total_time);
  U_alpha_2 = propagator_eig(Hamiltonian_alpha,total_time);
end


v= ones(1,timepoints);
for iTime = 1:timepoints
  
  switch EXPERIMENT
    case FID
      U_ = U_beta'*U_alpha;
      v(iTime) = vecDensityMatrixT*U_(:);
      
    case HAHN
      U_ = U_beta'*U_alpha'*U_beta*U_alpha;
      v(iTime) = vecDensityMatrixT*U_(:);
      
    case CPMG
      U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
      
      v(iTime) = vecDensityMatrixT*U_(:);
      
      U_beta_2 = dU_beta*U_beta_2;
      U_alpha_2 = dU_alpha*U_alpha_2;
      
    case CPMG_CONST
      
      U_beta = propagator_eig(Hamiltonian_beta,(iTime-1)*dt);
      U_alpha = propagator_eig(Hamiltonian_alpha,(iTime-1)*dt);
      
      U_beta_2 = propagator_eig(Hamiltonian_beta, total_time/4-(iTime-1)*dt);
      U_alpha_2 = propagator_eig(Hamiltonian_alpha, total_time/4-(iTime-1)*dt);
      U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
      
      v(iTime) = vecDensityMatrixT*U_(:);
      
      
    case CPMG_2D
      
      U_beta_2 = eye(nStates);
      U_alpha_2 = eye(nStates);
      
      for jTime = 1:timepoints
        
        U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
        
        v(iTime,jTime) = vecDensityMatrixT*U_(:);
        
        U_beta_2 = dU_beta*U_beta_2;
        U_alpha_2 = dU_alpha*U_alpha_2;
      end
      
  end
  
  
  U_beta = dU_beta*U_beta;
  U_alpha = dU_alpha*U_alpha;
  
  
end
 
if EXPERIMENT == CPMG_2D
    Signal = reshape(v.',1,[]);
else
    Signal = v;
end

if any(abs(v) - 1 > 1e-9)
  error('Coherence error: coherence cannot be larger than 1.');
end 
end


% ========================================================================
% Subcluster Function
% ========================================================================

function Indices = findSubclusters(Clusters,clusterSize,iCluster,reference_clusterSize)
% from Indices{subCluster_size} = list of all jCluster such that Clusters{subCluster_size}(jCluster,:) is a subcluster of Clusters{clusterSize}(iCluster,:)
% to
% Given the ith cluster of size clusterSize, 
% Indices(jCluster,subCluster_size) = jth cluster of sizesubCluster_size
% that is a subcluster of the ith cluster of size clusterSize
%
% SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
% the jth cluster of size subCluster_size that is a subcluster of
% the ith ccluster of size clusterSize.
% In a valid cluster, all indices must be natural numbers.


Indices = zeros(nchoosek(reference_clusterSize, ceil(reference_clusterSize/2)),reference_clusterSize);
jCluster = 0;

% Each cluster is a subset of itself.
Indices(1,clusterSize)=iCluster;

if clusterSize==1
  % All non-empty subsets have been found.
  return;
end

% Loop over all cluster sizes up to clusterSize.
for index = 1:clusterSize
  % CluserArray(iCluster,:,clusterSize) = nuclear indices.
  Cluster = Clusters(iCluster, 1:clusterSize ,clusterSize);
%   Cluster(Cluster==0) = [];
  
  % Remove one element labeled by index from Clusters{clusterSize}(iCluster:) to get a sub-cluster with one less element.  
  SubCluster = [Cluster(1:index-1), Cluster(index+1:end)];
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:clusterSize -1, clusterSize -1)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2) );
  
  % Skip if there is no subclustere.
  if isempty(subclusterIndex)
    continue;
  end
  
  jCluster = jCluster + 1;
  Indices(jCluster,clusterSize -1) = subclusterIndex;
  
  if clusterSize > 2
    % Given the ith cluster of size clusterSize,
    % Indices(jCluster,subCluster_size) = jth cluster of sizesubCluster_size
    % that is a subcluster of the ith cluster of size clusterSize
    
    % Recursively find smaller subclusters
    Subindices = findSubclusters(Clusters,clusterSize -1,subclusterIndex,reference_clusterSize);
    
    % Assign subcluster indices. 
    Indices(:,1:clusterSize-1) = Subindices(:,1:clusterSize-1);
    
  end
  
end

end


% ========================================================================
% Subcluster Function gpu
% ========================================================================

function Indices = findSubclusters_gpu(Clusters,clusterSize,iCluster,reference_clusterSize)
% from Indices{subCluster_size} = list of all jCluster such that Clusters{subCluster_size}(jCluster,:) is a subcluster of Clusters{clusterSize}(iCluster,:)
% to
% Given the ith cluster of size clusterSize, 
% Indices(jCluster,subCluster_size) = jth cluster of sizesubCluster_size
% that is a subcluster of the ith cluster of size clusterSize
%
% SubclusterIndices_clusterSize(jCluster,subCluster_size, iCluster) =
% the jth cluster of size subCluster_size that is a subcluster of
% the ith ccluster of size clusterSize.
% In a valid cluster, all indices must be natural numbers.


Indices = zeros(nchoosek(reference_clusterSize, ceil(reference_clusterSize/2)),reference_clusterSize);
jCluster = 0;

% Each cluster is a subset of itself.
Indices(1,clusterSize)=iCluster;
    
numberSubClusters = NchooseK(clusterSize,1:clusterSize);



switch clusterSize
  case 1
    % All non-empty subsets have been found.
    return;
    
  case 3
    possibleSubClusters_2 = zeros(numberSubClusters(2),clusterSize);
%     possibleSubClusters_2 = [0,1,1; 1, 0,1; 1,1,0];
    
  case 4
    possibleSubClusters_2 = zeros(numberSubClusters(2),clusterSize);
    possibleSubClusters_3 = zeros(numberSubClusters(3),clusterSize);
%     
%     possibleSubClusters_2 = [0,0,1,1; 0,1,0,1; 1,0,0,1; 0,1,1,0; 1,0,1,0; 1,1,0,0];
%     possibleSubClusters_3 = [0,1,1,1; 1,0,1,1; 1,1,0,1; 1,1,1,0];
   
  case 5
    possibleSubClusters_2 = zeros(numberSubClusters(2),clusterSize);
    possibleSubClusters_3 = zeros(numberSubClusters(3),clusterSize); 
    possibleSubClusters_4 = zeros(numberSubClusters(4),clusterSize); 
    
  case 6
    possibleSubClusters_2 = zeros(numberSubClusters(2),clusterSize);
    possibleSubClusters_3 = zeros(numberSubClusters(3),clusterSize); 
    possibleSubClusters_4 = zeros(numberSubClusters(4),clusterSize);  
    possibleSubClusters_5 = zeros(numberSubClusters(5),clusterSize); 
end  

subcluster_count = zeros(1,clusterSize);
for isc = 1:2^(clusterSize)-1
  
  subcluster_str = dec2bin(isc);
  subcluster_str = pad(subcluster_str,clusterSize,'left','0');
  sc_size = sum(subcluster_str=='1');
  
  switch sc_size
    
    case 2
      subcluster_count(sc_size) = subcluster_count(sc_size) + 1; 
      possibleSubClusters_2(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 3
      subcluster_count(sc_size) = subcluster_count(sc_size) + 1; 
      possibleSubClusters_3(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 4
      subcluster_count(sc_size) = subcluster_count(sc_size) + 1; 
      possibleSubClusters_4(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
    case 5
      subcluster_count(sc_size) = subcluster_count(sc_size) + 1; 
      possibleSubClusters_5(subcluster_count(sc_size),1:clusterSize) = subcluster_str=='1';
  end
  
end

Cluster = Clusters(iCluster, 1:clusterSize ,clusterSize);
Indices(1:clusterSize,1) = Cluster;

% for jCluster = 1:numberSubClusters(1)
%   
%   SubCluster = Cluster(possibleSubClusters_1==1);
%   
%   % Search for a valid cluster that equals the subcluster.
%   Search = Clusters(:, 1:clusterSize -1, clusterSize -1)==SubCluster;
%  
%   % Locate where the subcluster is.
%   subclusterIndex = find( all(Search,2));
%   
%   Indices(jCluster,1) = subclusterIndex;
% end

if clusterSize == 2
  return;
end
  
for jCluster = 1:numberSubClusters(2)
  
  SubCluster = Cluster(possibleSubClusters_2(jCluster,:)==1);
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:2, 2)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2));
 
  if isempty(subclusterIndex)
    subclusterIndex = 0;
  end
  
  Indices(jCluster,2) = subclusterIndex;
end

if clusterSize == 3
  return;
end
  
for jCluster = 1:numberSubClusters(3)
  
  SubCluster = Cluster(possibleSubClusters_3(jCluster,:)==1);
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:3, 3)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2));
 
  if isempty(subclusterIndex)
    subclusterIndex = 0;
  end
  
  Indices(jCluster,3) = subclusterIndex;
end

%--------------------------------------------------------------------------

isize = 4;  
if clusterSize == isize
  return;
end

for jCluster = 1:numberSubClusters(isize)
  
  SubCluster = Cluster(possibleSubClusters_4(jCluster,:)==1);
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:isize, isize)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2));
 
  if isempty(subclusterIndex)
    subclusterIndex = 0;
  end
  
  Indices(jCluster,isize) = subclusterIndex;
end
%--------------------------------------------------------------------------

isize = 5;  
if clusterSize == isize
  return;
end

for jCluster = 1:numberSubClusters(isize)
  
  SubCluster = Cluster(possibleSubClusters_5(jCluster,:)==1);
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters(:, 1:isize, isize)==SubCluster;
 
  % Locate where the subcluster is.
  subclusterIndex = find( all(Search,2));
 
  if isempty(subclusterIndex)
    subclusterIndex = 0;
  end
  
  Indices(jCluster,isize) = subclusterIndex;
end

end
% ========================================================================
% Calculate propagator using diagonalization
% ========================================================================

function U = propagator_eig(Ham,t)
%Ham = (Ham+Ham')/2; % "hermitianize" Hamiltonian
[EigenVectors, EigenValues] = eig(Ham);
Udiag = exp(-2i*pi*diag(EigenValues)*t);
U = EigenVectors*diag(Udiag)*EigenVectors';
end


% ========================================================================
% New Function TEMPORARY
% ========================================================================

function Indices = findSubclusters0(Clusters,clusterSize,iCluster)
% Indices{size} = list of all jCluster such that Clusters{size}(jCluster,:) is a subcluster of Clusters{clusterSize}(iCluster,:)

% In a valid cluster, all indices must be natural numbers.
if Clusters{clusterSize}(iCluster,1) == 0
  Indices = [];
  return;
end

% Each cluster is a subset of itself.
Indices{clusterSize}=iCluster;

if clusterSize==1
  % All non-empty subsets have been found.
  return;
end

% Loop over all cluster sizes up to clusterSize.
for index = 1:clusterSize
  
  % Remove one element labeled by index from Clusters{clusterSize}(iCluster:) to get a sub-cluster with one less element.  
  SubCluster = [Clusters{clusterSize}(iCluster,1:index-1),Clusters{clusterSize}(iCluster,index+1:end)];
  
  % Search for a valid cluster that equals the subcluster.
  Search = Clusters{clusterSize -1}==SubCluster;
  for ii = 2:size(Search,2)
    Search(:,1) = Search(:,1).*Search(:,ii);
  end
  subclusterIndex = find(Search(:,1)==1);
  
  if isempty(subclusterIndex)
    continue;
  end
  Indices{clusterSize -1} = [Indices{clusterSize -1} , subclusterIndex];
  
  if clusterSize > 2
    Subindices = findSubclusters0(Clusters,clusterSize -1,subclusterIndex);
    for isize = (clusterSize-2):-1:1
      if ~isempty(Subindices{isize})
        Indices{isize} = [Indices{isize},Subindices{isize}];
      end
    end
  end
  
end

for isize = 1:clusterSize
  Indices{isize} = unique(Indices{isize});
end
end


% ========================================================================
% New Function
% ========================================================================
function Coherences = propagate0(System, Method, DensityMatrix,timepoints,dt,Hamiltonian_beta,Hamiltonian_alpha,isotopeProbability, linearTimeAxis,verbose)
vecDensityMatrixT = reshape(DensityMatrix.',1,[]);

if strcmp('time-domain',Method.propagationDomain)
  
if linearTimeAxis
  dU_beta = propagator_eig(Hamiltonian_beta,dt);
  dU_alpha = propagator_eig(Hamiltonian_alpha,dt);
    
  nStates = length(Hamiltonian_beta);
  U_beta = eye(nStates);
  U_alpha = eye(nStates);
  
  if strcmp(System.experiment,'CPMG')
    U_beta_2 = eye(nStates);
    U_alpha_2 = eye(nStates);
  end
  
  if strcmp(System.experiment,'CPMG-const')
    U_beta_2 = propagator_eig(Hamiltonian_beta,System.total_time);
    U_alpha_2 = propagator_eig(Hamiltonian_alpha,System.total_time);
  end
end

Signal_= ones(1,timepoints);
for iTime = 1:timepoints
  if ~linearTimeAxis
    t = System.Time(iTime);
    U_beta = propagator_eig(Hamiltonian_beta,t);
    U_alpha = propagator_eig(Hamiltonian_alpha,t);
  end
  
  if strcmp(System.experiment,'FID')
    U_ = U_beta'*U_alpha;
    Signal_(iTime) = vecDensityMatrixT*U_(:);
    
  elseif strcmp(System.experiment,'Hahn')
    U_ = U_beta'*U_alpha'*U_beta*U_alpha;
    Signal_(iTime) = vecDensityMatrixT*U_(:);
  elseif strcmp(System.experiment,'CPMG')
    
    U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
    
    Signal_(iTime) = vecDensityMatrixT*U_(:);
    
    U_beta_2 = dU_beta*U_beta_2;
    U_alpha_2 = dU_alpha*U_alpha_2;
    
  elseif strcmp(System.experiment,'CPMG-const')
    
    U_beta = propagator_eig(Hamiltonian_beta,(iTime-1)*dt);
    U_alpha = propagator_eig(Hamiltonian_alpha,(iTime-1)*dt);

    U_beta_2 = propagator_eig(Hamiltonian_beta,System.total_time/4-(iTime-1)*dt);
    U_alpha_2 = propagator_eig(Hamiltonian_alpha,System.total_time/4-(iTime-1)*dt);
    U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;
      
    Signal_(iTime) = vecDensityMatrixT*U_(:);
    
%     U_beta_2 = dU_beta'*U_beta_2;
%     U_alpha_2 = dU_alpha'*U_alpha_2;
    
  elseif strcmp(System.experiment,'CPMG-2D')
    
    U_beta_2 = eye(nStates);
    U_alpha_2 = eye(nStates);
    
    for jTime = 1:timepoints
    
      U_ = U_alpha_2'*U_beta_2'  *  (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;  
      
      Signal_(iTime,jTime) = vecDensityMatrixT*U_(:);
      
      U_beta_2 = dU_beta*U_beta_2;
      U_alpha_2 = dU_alpha*U_alpha_2;
    end
   
  else
    error('The experiment ''%s'' is not supported.',System.experiment);
  end
  
%   Signal_(iTime) = vecDensityMatrixT*U_(:);
  
  if linearTimeAxis
    U_beta = dU_beta*U_beta;
    U_alpha = dU_alpha*U_alpha;
  end
  
  
end
 
elseif strcmp('frequency-domain',Method.propagationDomain)
 
% Diagonalize sub-Hamiltonians for alpha and beta electron manifolds
[Vec_beta, E_beta] = eig(2*pi*Hamiltonian_beta);
[Vec_alpha, E_alpha] = eig(2*pi*Hamiltonian_alpha);
E_beta = diag(E_beta);
E_alpha = diag(E_alpha);
M = Vec_alpha'*Vec_beta; % Mims matrix = <a|b> overlap matrix
Madj = M';

prefactor = 1; % include orienation weights etc here

G = prefactor*M; % prepared density matrix before first varied evolution time
D = M; % detection matrix after last varied evolution time
T1left = Madj; % left transfer matrix due to pi pulse
T1right = Madj; % right transfer matrix due to pi pulse

x=2^20;
nPoints = x;

IncSchemeID = 2; % for 2p ESEEM
buffRe = zeros(1,nPoints); % real part of spectral histogram
buffIm = zeros(1,nPoints); % imag part of spectral histogram
dt = System.dt/2/pi; % time step, in microseconds

sf_peaks(IncSchemeID,buffRe,buffIm,dt,[1 2],[2 1],E_alpha,E_beta,G,D,T1left,T1right);

Signal_ =ifft(buffRe+1i*buffIm);

Signal_ = Signal_(1:timepoints);
Signal_=Signal_./Signal_(1);

elseif strcmp('frequency-domain --matlab',Method.propagationDomain)
  
  nStates = length(Hamiltonian_beta);
  
  [Vec_beta, E_beta] = eig(2*pi*Hamiltonian_beta);
  [Vec_alpha, E_alpha] = eig(2*pi*Hamiltonian_alpha);
  E_beta = diag(E_beta);
  E_alpha = diag(E_alpha);
  dE_beta = E_beta - E_beta.';
  dE_alpha = E_alpha - E_alpha.';
  M = Vec_alpha'*Vec_beta; % Mims matrix
  Madj = M';
  
  % frequency domain variables
  n_omega = 2^20 - 1;

  domega = 2*pi/dt/(n_omega-1);
  
  
  omega_amp = zeros(1,n_omega);
  
  for ii = 1:nStates
    for jj = 1:nStates
      for kk = 1:nStates
        for ll = 1:nStates

          amplitude = Madj(ii,jj) * M(jj,kk) * Madj(kk,ll) * M(ll,ii)/2;
          omega_ = -( dE_beta(ii,kk) + dE_alpha(jj,ll));
          idx = mod(round(omega_./domega),n_omega) + 1;
          
          omega_amp(idx) = omega_amp(idx) + amplitude;
          
        end
      end
    end
  end

  Signal_ = ifft(omega_amp,n_omega);
  Signal_ = Signal_(1:timepoints);
  Signal_ =Signal_ /Signal_(1);
end
if Method.gpu && gpuDeviceCount > 0
  Signal_ = gather(Signal_);
end

if any(abs(Signal_) - 1 > 1e-9)
  error('Coherence error: coherence cannot be larger than 1.');
end
% Coherences = 1 + isotopeProbability*(Signal_ - 1);
Coherences = Signal_;
end
