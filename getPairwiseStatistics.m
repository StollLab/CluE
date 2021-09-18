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
if System.Methyl.include
  pwstat.Methyl_Data.ID = zeros(Nuclei.Methyl_Data.number_methyls,1);
  pwstat.Methyl_Data.number_methyls = 0;
end


% Methyl Groups
isMethylCarbon = strcmp(Nuclei.Type,'CH3');
isMethylHydron = strcmp(Nuclei.Type,'CH3_1H');

pwstat.Methyl_Data.number_methyls = sum(isMethylCarbon);
pwstat.Methyl_Data.ID = double(Nuclei.Index(isMethylCarbon));

deltaX = pwstat.Coordinates(:,1)-pwstat.Coordinates(:,1)';
deltaY = pwstat.Coordinates(:,2)-pwstat.Coordinates(:,2)';
deltaZ = pwstat.Coordinates(:,3)-pwstat.Coordinates(:,3)';

pwstat.DistanceMatrix = sqrt(deltaX.^2 + deltaY.^2 + deltaZ.^2);
pwstat.Cylindrical_DistanceMatrix = sqrt(deltaX.^2 + deltaY.^2);

pwstat.ThetaMatrix = acos(deltaZ./pwstat.DistanceMatrix);
pwstat.PhiMatrix = acos(deltaX./pwstat.Cylindrical_DistanceMatrix) ...
  + pi*(deltaY<0);

pwstat.PhiMatrix(deltaY==0) = 0;
 

% Coordinates
pwstat.Distance = vecnorm(pwstat.Coordinates,2,2);
pwstat.Cylindrical_Distance = vecnorm(pwstat.Coordinates(:,1:2),2,2);
pwstat.theta = acos(pwstat.Coordinates(:,3)./pwstat.Distance);
cosTheta2 = (pwstat.Coordinates(:,3)./pwstat.Distance).^2;

pwstat.maxDistance = max(pwstat.Distance,pwstat.Distance');

cos_phi= pwstat.Coordinates(:,1)./pwstat.Cylindrical_Distance;
pwstat.phi = acos(cos_phi) + (1 - sign(pwstat.Coordinates(:,2)))*pi/2;


pwstat.gamma_e = -System.Electron.g*System.muB/System.hbar;
pwstat.gamma_n = Nuclei.Nuclear_g*System.muN/System.hbar;


% Hyperfine
pwstat.Hyperfine_perpendicular = Nuclei.FermiContact ...
  - System.hbar^2.*(System.mu0/4/pi).*pwstat.gamma_n'.*pwstat.gamma_e .* ...
  pwstat.Distance.^-3; % J

pwstat.Hyperfine_perpendicular = pwstat.Hyperfine_perpendicular/(System.h); % Hz.
pwstat.DeltaHyperfine_perpendicular = pwstat.Hyperfine_perpendicular - pwstat.Hyperfine_perpendicular';

pwstat.Hyperfine = pwstat.Hyperfine_perpendicular.*(1-3*cosTheta2);
pwstat.DeltaHyperfine = pwstat.Hyperfine - pwstat.Hyperfine';


% Nuclear Dipole-Dipole
b_fliflopPerp = -(0.25*(System.mu0/4/pi).* ...
  (Nuclei.Nuclear_g'.*Nuclei.Nuclear_g).*System.muN^2) ...
  ./pwstat.DistanceMatrix.^3; % J.

b_fliflopPerp = b_fliflopPerp/(System.h); % Hz.
b = b_fliflopPerp.*(1-3*cos(pwstat.ThetaMatrix).^2);
pwstat.Nuclear_Dipole_perpendicular = 4*b_fliflopPerp; % Hz.
pwstat.Nuclear_Dipole = 4*b; % Hz.


pwstat.Nuclear_Dipole_x_iy_Z = zeroDiag(...
  pwstat.Nuclear_Dipole_perpendicular.* ...
  exp(1i*pwstat.PhiMatrix) .* ...
  sin(pwstat.ThetaMatrix) .* cos(pwstat.ThetaMatrix)...
);
if any(isnan(pwstat.Nuclear_Dipole_x_iy_Z(:)))
  error(['Error in getPairwiseStatistics(): ', ...
    'Nuclear_Dipole_x_iy_Z contains NANs.']);
end

% |bA_{mav}|^{1/2}.
pwstat.Amax = max( abs(pwstat.Hyperfine_perpendicular), abs(pwstat.Hyperfine_perpendicular'));
pwstat.bAmax = sqrt(pwstat.Amax.*abs(pwstat.Nuclear_Dipole));

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
pwstat.Frequency_Pair_p =  b_fliflopPerp.*sqrt(1+cp.^2); % Hz
pwstat.Frequency_Pair   =  b.*sqrt(1+c.^2); % Hz

% delta_gm_gn
pwstat.Same_g = Nuclei.Nuclear_g == Nuclei.Nuclear_g';

TM = System.TMguess;
pwstat.GaussianRMSD_p = getGaussianRMSD(modDepth_p,pwstat.Frequency_Pair_p,TM);
pwstat.GaussianRMSD = getGaussianRMSD(modDepth,pwstat.Frequency_Pair,TM);

  pwstat.Nuclear_Dipole = inf2nan(pwstat.Nuclear_Dipole);
  pwstat.Nuclear_Dipole_perpendicular = inf2nan(pwstat.Nuclear_Dipole_perpendicular);
  pwstat.Amax = inf2nan(pwstat.Amax);
  pwstat.bAmax = inf2nan(pwstat.bAmax);
  pwstat.Frequency_Pair_p = inf2nan(pwstat.Frequency_Pair_p);
  pwstat.Frequency_Pair = inf2nan(pwstat.Frequency_Pair);
  pwstat.Modulation_Depth_p = inf2nan(pwstat.Modulation_Depth_p);
  pwstat.Modulation_Depth = inf2nan(pwstat.Modulation_Depth);
  pwstat.GaussianRMSD_p = inf2nan(pwstat.GaussianRMSD_p);
  pwstat.GaussianRMSD = inf2nan(pwstat.GaussianRMSD);
%   pwstat.DistanceMatrix(pwstat.DistanceMatrix == 0) = nan;
  
  if System.Methyl.include
    methylList_ = pwstat.Methyl_Data.ID(pwstat.Methyl_Data.ID>0);
    for ispin  = methylList_
      hydrons_ = ispin + [1 ,2, 3];
      
      pwstat.Nuclear_Dipole(:,ispin) = maxabs(pwstat.Nuclear_Dipole(:,hydrons_)')';
      pwstat.Nuclear_Dipole_perpendicular(:,ispin) = maxabs( pwstat.Nuclear_Dipole_perpendicular(:,hydrons_)')';
      pwstat.Amax(:,ispin) = maxabs( pwstat.Amax(:,hydrons_)')';
      pwstat.bAmax(:,ispin) = maxabs( pwstat.bAmax(:,hydrons_)')';
      pwstat.Frequency_Pair_p(:,ispin) = maxabs( pwstat.Frequency_Pair_p(:,hydrons_)')';
      pwstat.Frequency_Pair(:,ispin) = maxabs( pwstat.Frequency_Pair(:,hydrons_)')';
      pwstat.Modulation_Depth_p(:,ispin) = maxabs( pwstat.Modulation_Depth_p(:,hydrons_)')';
      pwstat.Modulation_Depth(:,ispin) = maxabs( pwstat.Modulation_Depth(:,hydrons_)')';
      pwstat.GaussianRMSD_p(:,ispin) = maxabs( pwstat.GaussianRMSD_p(:,hydrons_)')';
      pwstat.GaussianRMSD(:,ispin) = maxabs( pwstat.GaussianRMSD(:,hydrons_)')';
      
      
      %       pwstat.DistanceMatrix(:,ispin) = minabs(pwstat.DistanceMatrix(:,hydrons_)')';
      %       pwstat.Cylindrical_DistanceMatrix(:,ispin) = minabs(pwstat.Cylindrical_DistanceMatrix(:,hydrons_)')';
      %       pwstat.ThetaMatrix(:,ispin) = mean(pwstat.ThetaMatrix(:,hydrons_)')';
      %       pwstat.PhiMatrix(:,ispin) = mean(pwstat.PhiMatrix(:,hydrons_)')';
      
      for jspin  = pwstat.Methyl_Data.ID'
        if jspin >= ispin, break; end
        hydrons_j = jspin + [1 ,2, 3];
        
        pwstat.Nuclear_Dipole(jspin,ispin) = supabs(pwstat.Nuclear_Dipole(hydrons_j,hydrons_));
        pwstat.Nuclear_Dipole_perpendicular(:,ispin) = supabs( pwstat.Nuclear_Dipole_perpendicular(hydrons_j,hydrons_));
        pwstat.Amax(jspin,ispin) = supabs( pwstat.Amax(hydrons_j,hydrons_));
        pwstat.bAmax(jspin,ispin) = supabs( pwstat.bAmax(hydrons_j,hydrons_));
        pwstat.Frequency_Pair_p(jspin,ispin) = supabs( pwstat.Frequency_Pair_p(hydrons_j,hydrons_));
        pwstat.Frequency_Pair(jspin,ispin) = supabs( pwstat.Frequency_Pair(hydrons_j,hydrons_));
        pwstat.Modulation_Depth_p(jspin,ispin) = supabs( pwstat.Modulation_Depth_p(hydrons_j,hydrons_));
        pwstat.Modulation_Depth(jspin,ispin) = supabs( pwstat.Modulation_Depth(hydrons_j,hydrons_));
        pwstat.GaussianRMSD_p(jspin,ispin) = supabs( pwstat.GaussianRMSD_p(hydrons_j,hydrons_));
        pwstat.GaussianRMSD(jspin,ispin) = supabs( pwstat.GaussianRMSD(hydrons_j,hydrons_));
      end
      
    end
    % Use the symmetric nature of the separation matrix to set the ispin > jspin entries.
    pwstat.Nuclear_Dipole(ispin,:) = pwstat.Nuclear_Dipole(:,ispin);
    pwstat.Nuclear_Dipole_perpendicular(ispin,:) = pwstat.Nuclear_Dipole_perpendicular(:,ispin);
    pwstat.Amax(ispin,:) = pwstat.Amax(:,ispin);
    pwstat.bAmax(ispin,:) = pwstat.bAmax(:,ispin);
    pwstat.Frequency_Pair_p(ispin,:) = pwstat.Frequency_Pair_p(:,ispin);
    pwstat.Frequency_Pair(ispin,:) = pwstat.Frequency_Pair(:,ispin);
    pwstat.Modulation_Depth_p(ispin,:) = pwstat.Modulation_Depth_p(:,ispin);
    pwstat.Modulation_Depth(ispin,:) = pwstat.Modulation_Depth(:,ispin);
    pwstat.GaussianRMSD_p(ispin,:) = pwstat.GaussianRMSD_p(:,ispin);
    pwstat.GaussianRMSD(ispin,:) = pwstat.GaussianRMSD(:,ispin);
    
%     pwstat.DistanceMatrix(ispin,:) = pwstat.DistanceMatrix(:,ispin);
%     pwstat.Cylindrical_DistanceMatrix(ispin,:) = pwstat.Cylindrical_DistanceMatrix(:,ispin);
%     pwstat.ThetaMatrix(ispin,:) =pwstat.ThetaMatrix(:,ispin);
%     pwstat.PhiMatrix(ispin,:) = pwstat.PhiMatrix(:,ispin);
    
end

end



