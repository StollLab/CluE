function gofr = getRadialDistributionFunction(DistanceMatrix)

binWidth = 1e-11;
binEdges = 0:binWidth:max(DistanceMatrix(:));
binCenters = binEdges(1:end-1) + 0.5*binWidth; 
gofr = histcounts(DistanceMatrix(:),binEdges);

gofr(1) = 0;
gofr = 0.5*gofr./(binCenters*1e10)./(binCenters*1e10);
end