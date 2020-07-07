% ========================================================================
% Calculate statistics for each bath spin pair. 
% ========================================================================

function pwstat = getPairwiseStatistics(System, Nuclei)

N = Nuclei.number;

% Initialize separation matrix.
pwstat.Coordinates = Nuclei.Coordinates;
pwstat.DistanceMatrix = zeros(N);
pwstat.Cylindrical_DistanceMatrix = zeros(N);
pwstat.ThetaMatrix = zeros(N);
pwstat.PhiMatrix = zeros(N); 

% Loop over all nuclear spins.
for ispin = 1:N
 
  % Loop over all nuclear spins with a higher index than ispin.
  for jspin = ispin+1:N
    
    deltaR_ = pwstat.Coordinates(ispin,:)-pwstat.Coordinates(jspin,:);
    
    % Set separation matrix entry.
    pwstat.DistanceMatrix(ispin,jspin) = norm(deltaR_);
    pwstat.Cylindrical_DistanceMatrix(ispin,jspin) = norm(deltaR_(1:2));
    
    cosTheta_ = deltaR_(3)/pwstat.DistanceMatrix(ispin,jspin);
    pwstat.ThetaMatrix(ispin,jspin) = acos(cosTheta_);
    
    cosPhi_ = deltaR_(1)/pwstat.Cylindrical_DistanceMatrix(ispin,jspin);
    pwstat.PhiMatrix(ispin,jspin) = acos(cosPhi_) + (1 - sign( deltaR_(2) ))*pi/2;
    
    % Use the symmetric nature of the separation matrix to set the ispin > jspin entries. 
    pwstat.DistanceMatrix(jspin,ispin) = pwstat.DistanceMatrix(ispin,jspin);
    pwstat.Cylindrical_DistanceMatrix(jspin,ispin) = pwstat.Cylindrical_DistanceMatrix(ispin,jspin);
    pwstat.ThetaMatrix(jspin,ispin) = pwstat.ThetaMatrix(ispin,jspin);
    pwstat.PhiMatrix(jspin,ispin) = pwstat.PhiMatrix(ispin,jspin);
  end
end

% Coordinates
pwstat.Distance = vecnorm(pwstat.Coordinates,2,2);
pwstat.Cylindrical_Distance = vecnorm(pwstat.Coordinates(:,1:2),2,2);
pwstat.theta = acos(pwstat.Coordinates(:,3)./pwstat.Distance);
cosTheta2 = (pwstat.Coordinates(:,3)./pwstat.Distance).^2;

cos_phi= pwstat.Coordinates(:,1)./pwstat.Cylindrical_Distance;
pwstat.phi = acos(cos_phi) + (1 - sign(pwstat.Coordinates(:,2)))*pi/2;


pwstat.gamma_e = -System.Electron.g*System.muB/System.hbar;
pwstat.gamma_n = Nuclei.Nuclear_g*System.muN/System.hbar;


% Hyperfine
pwstat.Hyperfine_perpendicular = Nuclei.FermiContact ...
  - System.hbar.*(System.mu0/4/pi).*pwstat.gamma_n.*pwstat.gamma_e.*pwstat.Distance.^-3; % Hz
pwstat.DeltaHyperfine_perpendicular = pwstat.Hyperfine_perpendicular - pwstat.Hyperfine_perpendicular';

pwstat.Hyperfine = pwstat.Hyperfine_perpendicular.*(1-3*cosTheta2);
pwstat.DeltaHyperfine = pwstat.Hyperfine - pwstat.Hyperfine';

% Nuclear Dipole-Dipole
b = 0.25*(System.mu0/4/pi).*(Nuclei.Nuclear_g'.*Nuclei.Nuclear_g).*System.muN^2; % J m^3.
b = -b./pwstat.DistanceMatrix.^3; % J.
b = b/(System.h); % Hz.
bp = b.*(1-3*cos(pwstat.ThetaMatrix).^2);
pwstat.Nuclear_Dipole_perpendicular = 4*b; % Hz.
pwstat.Nuclear_Dipole = 4*bp; % Hz.


% Hyperfine to nucler dipole-dipole ratio
cp = pwstat.DeltaHyperfine_perpendicular./pwstat.Nuclear_Dipole_perpendicular;
c = pwstat.DeltaHyperfine./pwstat.Nuclear_Dipole;

pwstat.DeltaHyperfine_over_Nuclear_Dipole_p = cp;
pwstat.DeltaHyperfine_over_Nuclear_Dipole = c;

% modulation depth
modDepth_p = 4*cp.^2./(1+cp.^2).^2;
modDepth = 4*c.^2./(1+c.^2).^2;
pwstat.Modulation_Depth_p = modDepth_p;
pwstat.Modulation_Depth = modDepth;

% modulation frequency
pwstat.Frequency_Pair_p =  2.*bp.*sqrt(1+cp.^2)./(2*pi);
pwstat.Frequency_Pair =  2.*b.*sqrt(1+c.^2)./(2*pi);

% delta_gm_gn
pwstat.Same_g = Nuclei.Nuclear_g == Nuclei.Nuclear_g';

TM = System.TMguess;
pwstat.GaussianRMSD_p = getGaussianRMSD(modDepth_p,pwstat.Frequency_Pair_p,TM);
pwstat.GaussianRMSD = getGaussianRMSD(modDepth,pwstat.Frequency_Pair,TM);

end



