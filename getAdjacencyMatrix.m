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
        Max_R_ = pwstat.DistanceMatrix <= Method.r0;
        Min_R_ = pwstat.DistanceMatrix > Method.r_min;
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*Max_R_.*Min_R_;
        
      case 'modulation'
        if Ori_cutoffs
          Min_Mod_ = pwstat.ModulationDepth >= Method.cutoff.modulation(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_Mod_;
        else
          Min_Mod_ = pwstat.ModulationDepth_p >= Method.cutoff.modulation(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_Mod_;
        end
        
      case 'dipole'
          if Ori_cutoffs
          Min_dipole_ = abs(pwstat.Nuclear_Dipole) >= Method.cutoff.dipole(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_dipole_;
          else
            Min_dipole_ = abs(pwstat.Nuclear_Dipole_perpendicular) >= Method.cutoff.dipole(isize);
            Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_dipole_;
          end
          
      case 'dipoleHalf'
        Min_dipole_ = Nuclei.Spin == 1/2;
        Min_dipole_ = Min_dipole_ & Min_dipole_';
        if Ori_cutoffs
          Min_dipole_ = ~Min_dipole_ | (Min_dipole_ & abs(pwstat.Nuclear_Dipole) >= Method.cutoff.dipoleHalf(isize));
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_dipole_;
        else
          Min_dipole_ = ~Min_dipole_ | (Min_dipole_ & abs(pwstat.Nuclear_Dipole_perpendicular) >= Method.cutoff.dipoleHalf(isize));
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_dipole_;
        end
          
      case 'dipoleOne'
        Min_dipole_ = Nuclei.Spin == 1;
        Min_dipole_ = Min_dipole_ & Min_dipole_';
        if Ori_cutoffs
          Min_dipole_ = ~Min_dipole_ | (Min_dipole_ & abs(pwstat.Nuclear_Dipole) >= Method.cutoff.dipoleOne(isize));
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_dipole_;
        else
          Min_dipole_ = ~Min_dipole_ | (Min_dipole_ & abs(pwstat.Nuclear_Dipole_perpendicular) >= Method.cutoff.dipoleOne(isize));
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_dipole_;
        end
        
      case 'bAmax'
        Min_bAmax_ = abs(pwstat.bAmax) >= Method.cutoff.bAmax(isize);
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_bAmax_;

      case 'minAmax'
        Min_Amax_ = abs(pwstat.Amax) >= Method.cutoff.minAmax(isize);
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_Amax_;  

      case 'maxAmax'
        Max_Amax_ = abs(pwstat.Amax) <= Method.cutoff.maxAmax(isize);
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*Max_Amax_;  
        
      case 'minimum-frequency'
        if Ori_cutoffs
          Min_Freq_ = abs(pwstat.Frequency_Pair) >  Method.cutoff.minimum_frequency(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_Freq_;
        else
          Min_Freq_ = abs(pwstat.Frequency_Pair_p) >  Method.cutoff.minimum_frequency(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_Freq_;
        end
        
      case 'delta hyperfine'
        if Ori_cutoffs
          Min_DeltaA_ = abs(pwstat.DeltaHyperfine) > Method.cutoff.hyperfine_inf(isize);
          Max_DeltaA_ = abs(pwstat.DeltaHyperfine) < Method.cutoff.hyperfine_sup(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_DeltaA_.*Max_DeltaA_;
        else
          Min_DeltaA_ = abs(pwstat.DeltaHyperfine_perpendicular) > Method.cutoff.hyperfine_inf(isize);
          Max_DeltaA_ = abs(pwstat.DeltaHyperfine_perpendicular) < Method.cutoff.hyperfine_sup(isize);
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_DeltaA_.*Max_DeltaA_;
        end
        
      case 'Gaussian RMSD'
        if Ori_cutoffs
          Min_gRMSD = pwstat.GaussianRMSD > Method.cutoff.gaussianRMSD;
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_gRMSD;
        else
          Min_gRMSD = pwstat.GaussianRMSD_p > Method.cutoff.gaussianRMSD;
          Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_gRMSD;
        end
      
      case 'radius_nonSpinHalf'
        MaxR_ = Nuclei.Spin == 1/2;
        MaxR_ = MaxR_ & MaxR_';
        MaxR_ = MaxR_ | pwstat.maxDistance <= Method.cutoff.radius_nonSpinHalf;
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*MaxR_;
        
%       case 'methyl only' 
%         if Method.cutoff.methylOnly(isize)
%           if System.Methyl.method == 2
%             isMethyl = strcmp(Nuclei.Type,'CH3_1H');
%           else
%             isMethyl = strcmp(Nuclei.Type,'CH3');
%           end
%           isCH3 = isMethyl + isMethyl' == 2;
%           Adjacency(:,:,isize) = Adjacency(:,:,isize).*isCH3;
%         end
        
      case 'same g'
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*pwstat.Same_g;
        
      case 'spin 1/2 only'
        Sele_ = Nuclei.Spin == 1/2;
        Sele_ = Sele_ & Sele_';
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*Sele_;
    end
  end
  
  if System.Methyl.method == 2
    % Protons within a methyl group should always be connected to each
    % other.
    Sele_ = Nuclei.MethylID > 0;
    Sele_ = Sele_ & Sele_';
    Sele_ = Sele_.*(Nuclei.MethylID==Nuclei.MethylID');
    Adjacency(:,:,isize) = Adjacency(:,:,isize) | Sele_;

  end
  
end
end