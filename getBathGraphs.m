function graphs = getBathGraphs(Nuclei)
N = Nuclei.number;

graphs.DistanceMatrix = zeros(N);
for ispin = 1:N
  for jspin = ispin+1:N
    graphs.DistanceMatrix(ispin,jspin) = norm(Nuclei.Coordinates(ispin,:) - Nuclei.Coordinates(jspin,:));
    graphs.DistanceMatrix(jspin,ispin) = graphs.DistanceMatrix(ispin,jspin);
  end
end
PairwiseStatistics = getPairwiseStatistics(System, Nuclei);


graphs.ValidPair = ones(N,N,Method.order);
for isize = 1:Method.order
  graphs.ValidPair(:,:,isize) = graphs.valid'*graphs.valid;
  graphs.ValidPair(:,:,isize) = graphs.ValidPair(:,:,isize).*PairwiseStatistics.Same_g;
  % exclude_spins = find( isnan(Nuclei.Spin));
  % Nuclei.ValidPair(exclude_spins,:) = 0;
  % Nuclei.ValidPair(:,exclude_spins) = 0;
  
  num_criteria = numel(Method.Criteria);
  for ii = 1:num_criteria
    switch Method.Criteria{ii}
      case 'neighbor'
        Max_Neighbor_ = PairwiseStatistics.DistanceMatrix <= Method.r0;
        Min_Neighbor_ = PairwiseStatistics.DistanceMatrix > Method.r_min;
        graphs.ValidPair(:,:,isize) = graphs.ValidPair(:,:,isize).*Max_Neighbor_.*Min_Neighbor_;
        
      case 'modulation'
        Min_Mod_ = PairwiseStatistics.ModulationDepth >= Method.cutoff.modulation(isize);
        graphs.ValidPair(:,:,isize) = graphs.ValidPair(:,:,isize).*Min_Mod_;
        
      case 'dipole'
        Min_dipole_ = abs(PairwiseStatistics.Nuclear_Dipole) >= Method.cutoff.dipole(isize);
        graphs.ValidPair(:,:,isize) = graphs.ValidPair(:,:,isize).*Min_dipole_;
        
      case 'minimum-frequency'
        Min_Freq_ = abs(PairwiseStatistics.Frequency_Pair) >  Method.cutoff.minimum_frequency(isize);
        graphs.ValidPair(:,:,isize) = graphs.ValidPair(:,:,isize).*Min_Freq_;
        
      case 'hyperfine'
        Min_DeltaA_ = abs(PairwiseStatistics.DeltaHyperfine) > Method.cutoff.hyperfine_inf(isize);
        Max_DeltaA_ = abs(PairwiseStatistics.DeltaHyperfine) < Method.cutoff.hyperfine_sup(isize);
        graphs.ValidPair(:,:,isize) = graphs.ValidPair(:,:,isize).*Min_DeltaA_.*Max_DeltaA_;
    end
  end
end

end