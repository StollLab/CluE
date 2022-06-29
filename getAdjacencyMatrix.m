function Adjacency = getAdjacencyMatrix(System, Nuclei,Method)

pwstat = Nuclei.Statistics;

% Orientation dependent cutoffs
Ori_cutoffs = Method.Ori_cutoffs;

% Get the highest spin value.
Nuclei.maxSpin = max(Nuclei.Spin);

% Initialize a matrix of valid edges.
Adjacency = true(Nuclei.number,Nuclei.number,Method.extraOrder);
for isize = 1:Method.extraOrder
  Adjacency(:,:,isize) = Nuclei.valid'*Nuclei.valid > 0;
  if ~Method.allowHDcoupling
    Adjacency(:,:,isize) = Adjacency(:,:,isize).*pwstat.Same_g;
  end
  % Loop over all criteria used to determine an edge.
  num_criteria = numel(Method.Criteria);
  for ii = 1:num_criteria
    switch Method.Criteria{ii}
      case 'distance'
        Max_R_ = pwstat.DistanceMatrix <= Method.neighborCutoff.rMax;
        Min_R_ = pwstat.DistanceMatrix > Method.neighborCutoff.rMin;
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*Max_R_.*Min_R_;
        
      case 'modulation'
        if Ori_cutoffs
          Min_Mod_ = pwstat.Modulation_Depth >= Method.neighborCutoff.modulation(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_Mod_;
        else
          Min_Mod_ = pwstat.Modulation_Depth_p ...
            >= Method.neighborCutoff.modulation(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_Mod_;
        end
        
      case 'dipole'
          if Ori_cutoffs
          Min_dipole_ = abs(pwstat.Nuclear_Dipole) ...
            >= Method.neighborCutoff.dipole(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_dipole_;
          else
            Min_dipole_ = abs(pwstat.Nuclear_Dipole_perpendicular) ...
              >= Method.neighborCutoff.dipole(isize);
            Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_dipole_;
          end
          
      case 'dipoleHalf'
        Min_dipole_ = Nuclei.Spin == 1/2;
        Min_dipole_ = Min_dipole_ & Min_dipole_';
        if Ori_cutoffs
          
          Min_dipole_ = ~Min_dipole_ ...
            | (Min_dipole_ & abs(pwstat.Nuclear_Dipole) ...
            >= Method.neighborCutoff.dipoleHalf(isize));
          
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_dipole_;
        else
          Min_dipole_ = ~Min_dipole_ ...
            | (Min_dipole_ & abs(pwstat.Nuclear_Dipole_perpendicular) ...
            >= Method.neighborCutoff.dipoleHalf(isize));
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_dipole_;
        end
          
      case 'dipoleOne'
        Min_dipole_ = Nuclei.Spin == 1;
        Min_dipole_ = Min_dipole_ & Min_dipole_';
        if Ori_cutoffs
          Min_dipole_ = ~Min_dipole_ | (Min_dipole_ ...
            & abs(pwstat.Nuclear_Dipole) >= Method.neighborCutoff.dipoleOne(isize));
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_dipole_;
        else
          Min_dipole_ = ~Min_dipole_ | (Min_dipole_ ...
            & abs(pwstat.Nuclear_Dipole_perpendicular) ...
            >= Method.neighborCutoff.dipoleOne(isize));
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_dipole_;
        end
        
      case 'bAmax'
        Min_bAmax_ = abs(pwstat.bAmax) >= Method.neighborCutoff.bAmax(isize);
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_bAmax_;

      case 'minAmax'
        Min_Amax_ = abs(pwstat.Amax) >= Method.neighborCutoff.minAmax(isize);
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_Amax_;  

      case 'maxAmax'
        Max_Amax_ = abs(pwstat.Amax) <= Method.neighborCutoff.maxAmax(isize);
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*Max_Amax_;  
        
      case 'minimum-frequency'
        if Ori_cutoffs
          Min_Freq_ = abs(pwstat.Frequency_Pair) ...
            > Method.neighborCutoff.minimum_frequency(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_Freq_;
        else
          Min_Freq_ = abs(pwstat.Frequency_Pair_p) ...
            > Method.neighborCutoff.minimum_frequency(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_Freq_;
        end
        
      case 'delta hyperfine'
        if Ori_cutoffs
          DeltaHyperfine = pwstat.Hyperfine - pwstat.Hyperfine';
          Min_DeltaA_ = abs(DeltaHyperfine) ...
            > Method.neighborCutoff.DeltaHyperfine(isize);
          %Max_DeltaA_ = abs(DeltaHyperfine) ...
          %  < Method.neighborCutoff.hyperfine_sup(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_DeltaA_;%.*Max_DeltaA_;
          clear('DeltaHyperfine');
        else
          DeltaHyperfine_perpendicular = pwstat.Hyperfine_perpendicular ...
            - pwstat.Hyperfine_perpendicular';
          Min_DeltaA_ = abs(DeltaHyperfine_perpendicular) ...
            > Method.neighborCutoff.hyperfine_inf(isize);
          Max_DeltaA_ = abs(DeltaHyperfine_perpendicular) ...
            < Method.neighborCutoff.hyperfine_sup(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_DeltaA_.*Max_DeltaA_;

          clear('DeltaHyperfine_perpendicular');
        end
        
      case 'Gaussian RMSD'
        if Ori_cutoffs
          Min_gRMSD = pwstat.GaussianRMSD > Method.neighborCutoff.gaussianRMSD;
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_gRMSD;
        else
          Min_gRMSD = pwstat.GaussianRMSD_p > Method.neighborCutoff.gaussianRMSD;
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_gRMSD;
        end
      
      case 'radius_nonSpinHalf'
        maxDistance = vecnorm(Nuclei.Coordinates,2,2);
        maxDistance = max(maxDistance,maxDistance');
        MaxR_ = Nuclei.Spin == 1/2;
        MaxR_ = MaxR_ & MaxR_';
        MaxR_ = MaxR_ | maxDistance <= Method.vertexCutoff.radius_nonSpinHalf;
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*MaxR_;
        clear('maxDistance');

        
      case 'same g'
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*pwstat.Same_g;
        
      case 'spin 1/2 only'
        Sele_ = Nuclei.Spin == 1/2;
        Sele_ = Sele_ & Sele_';
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*Sele_;
    end
  end
  
  if System.Methyl.method == 1
    isMethylHydron = Nuclei.MethylID > 0;
    Adjacency(isMethylHydron,:,isize) = false;
    Adjacency(:,isMethylHydron,isize) = false;
  elseif System.Methyl.method == 2
    % Protons within a methyl group should always be connected to each
    % other.
    isMethylHydron = Nuclei.MethylID > 0;
    isMethylHydron = (isMethylHydron + isMethylHydron')>0;
    isSameMethyl = (Nuclei.MethylID==Nuclei.MethylID');

    Sele_ = isMethylHydron.*isSameMethyl;
    Adjacency(:,:,isize) = (Adjacency(:,:,isize) + Sele_) > 0;

    
    isMethylCarbon = strcmp(Nuclei.Type,'CH3');
    Adjacency(isMethylCarbon,:,isize) = false;
    Adjacency(:,isMethylCarbon,isize) = false;

  end
  
end
end
