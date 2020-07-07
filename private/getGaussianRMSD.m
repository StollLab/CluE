function RMSD = getGaussianRMSD(modDepth,modFreq,TM)

x = TM.*modFreq*2*pi;
x2 = x.^2;

msd = (sqrt(2*pi)*(modDepth/8).^2).*(1/8*exp(-2*x2) - exp(-9/8*x2) + 7/2*exp(-1/2*x2) - 7*exp(-1/8*x2) +35/8);
RMSD = sqrt(msd);

end