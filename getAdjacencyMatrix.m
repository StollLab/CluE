function Adjacency = getAdjacencyMatrix(Nuclei,Method)

pwstat = Nuclei.Statistics;

% Orientation dependent cutoffs
Ori_cutoffs = Method.Ori_cutoffs;

% Get the highest spin value.
Nuclei.maxSpin = max(Nuclei.Spin);

% Initialize a matrix of valid edges.
Adjacency = ones(Nuclei.number,Nuclei.number,Method.order);
for isize = 1:Method.order
  Adjacency(:,:,isize) = Nuclei.valid'*Nuclei.valid > 0;
  Adjacency(:,:,isize) = Adjacency(:,:,isize).*pwstat.Same_g;
 
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
          
      case 'bAmax'
        Min_bAmax_ = abs(pwstat.bAmax) >= Method.cutoff.bAmax(isize);
        Adjacency(:,:,isize) = Adjacency(:,:,isize).*Min_bAmax_;
        
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
        
        
    end
  end
end
end