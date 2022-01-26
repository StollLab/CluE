% ========================================================================
% Calculate statistics for each bath spin pair.
% ========================================================================
%
% Input Arguments
%
% System = struct( ...
%  'Methyl', struct('include',bool), ...
%  'Electron' struct('g', double), ...
%  'muB', double, ...
%  'muN', double, ...
%  'mu0', double, ...
%  'hbar', double, ...
%  'h', double, ...
%  'TMguess', double);
%
% Method = struct(
%   'Criteria' , cell of char vectors, ...
%   'conserveMemory', bool);
%
% Nuclei = struct( 
%  'number', uint, ...
%  'Type', cell of char vectors, ...
%  'Index', array<uint>, ...
%  'Coordinates', array<double>
%  'Nuclear_g', array<double>, ...
%  'FermiContact', array<double>);
%
function pwstat = getPairwiseStatistics(System, Method, Nuclei)

ienum = 1;
STAT_SPATIAL = ienum; ienum = ienum +1;
STAT_HYPERFINE = ienum; ienum = ienum +1;
STAT_NUCLEAR_DIPOLE = ienum; ienum = ienum +1;
STAT_BAMAX = ienum; ienum = ienum +1;
STAT_ANALYTIC_HAHN= ienum; ienum = ienum +1;
STAT_NUM_ENUM = ienum-1;


if Method.conserveMemory
  stat_bools = false(STAT_NUM_ENUM,1);
else
  stat_bools = true(STAT_NUM_ENUM,1);
end
whichStatistics(Method);

N = Nuclei.number;

if System.Methyl.include
  % Methyl Groups
  isMethylCarbon = strcmp(Nuclei.Type,'CH3');
  %isMethylHydron = strcmp(Nuclei.Type,'CH3_1H');

  pwstat.Methyl_Data.number_methyls = sum(isMethylCarbon);
  pwstat.Methyl_Data.ID = double(Nuclei.Index(isMethylCarbon));

  methylList_ = pwstat.Methyl_Data.ID(pwstat.Methyl_Data.ID>0);
end


deltaX = Nuclei.Coordinates(:,1)-Nuclei.Coordinates(:,1)';
deltaY = Nuclei.Coordinates(:,2)-Nuclei.Coordinates(:,2)';
deltaZ = Nuclei.Coordinates(:,3)-Nuclei.Coordinates(:,3)';

pwstat.DistanceMatrix = sqrt(deltaX.^2 + deltaY.^2 + deltaZ.^2);
pwstat.Distance = vecnorm(Nuclei.Coordinates,2,2);

cosTheta2 = (Nuclei.Coordinates(:,3)./pwstat.Distance ).^2;

gamma_e = -System.Electron.g*System.muB/System.hbar;
gamma_n = Nuclei.Nuclear_g*System.muN/System.hbar;

% spatial
if stat_bools(STAT_SPATIAL)
  % Initialize separation matrix.
  pwstat.Coordinates = Nuclei.Coordinates;

  pwstat.Cylindrical_DistanceMatrix = sqrt(deltaX.^2 + deltaY.^2);

  pwstat.ThetaMatrix = acos(deltaZ./pwstat.DistanceMatrix);
  pwstat.PhiMatrix = acos(deltaX./pwstat.Cylindrical_DistanceMatrix) ...
    + pi*(deltaY<0);

  pwstat.PhiMatrix(deltaY==0) = 0;

  % Coordinates
  pwstat.Distance = vecnorm(Nuclei.Coordinates,2,2);
  pwstat.Cylindrical_Distance = vecnorm(Nuclei.Coordinates(:,1:2),2,2);
  pwstat.theta = acos(Nuclei.Coordinates(:,3)./pwstat.Distance);

  cos_phi= Nuclei.Coordinates(:,1)./pwstat.Cylindrical_Distance;
  pwstat.phi = acos(cos_phi) + (1 - sign(Nuclei.Coordinates(:,2)))*pi/2;

  pwstat.maxDistance = max(pwstat.Distance,pwstat.Distance');
end



% Hyperfine
if stat_bools(STAT_HYPERFINE)

  pwstat.Hyperfine_perpendicular = Nuclei.FermiContact ...
    - System.hbar^2.*(System.mu0/4/pi).*gamma_n'.*gamma_e .* ...
    pwstat.Distance.^-3; % J

  pwstat.Hyperfine_perpendicular ...
    = pwstat.Hyperfine_perpendicular/(System.h); % Hz.

  pwstat.DeltaHyperfine_perpendicular ...
    = pwstat.Hyperfine_perpendicular - pwstat.Hyperfine_perpendicular';

  pwstat.Hyperfine = pwstat.Hyperfine_perpendicular.*(1-3*cosTheta2);
  pwstat.DeltaHyperfine = pwstat.Hyperfine - pwstat.Hyperfine';
end

% Nuclear Dipole-Dipole
if stat_bools(STAT_NUCLEAR_DIPOLE)

  if isfield(pwstat,'ThetaMatrix')
    ThetaMatrix = pwstat.ThetaMatrix;
  else
    ThetaMatrix = acos(deltaZ./pwstat.DistanceMatrix);
  end

  if isfield(pwstat,'PhiMatrix')
    PhiMatrix = pwstat.PhiMatrix;
  else
    PhiMatrix = acos(deltaX./vecnorm(Nuclei.Coordinates(:,1:2),2,2)) ...
      + pi*(deltaY<0);

    PhiMatrix(deltaY==0) = 0;
  end

  b_fliflopPerp = -(0.25*(System.mu0/4/pi).* ...
    (Nuclei.Nuclear_g'.*Nuclei.Nuclear_g).*System.muN^2) ...
    ./pwstat.DistanceMatrix.^3; % J.

  b_fliflopPerp = b_fliflopPerp/(System.h); % Hz.

  b = b_fliflopPerp.*(1-3*cos(ThetaMatrix).^2);

  pwstat.Nuclear_Dipole_perpendicular = 4*b_fliflopPerp; % Hz.

  pwstat.Nuclear_Dipole = 4*b; % Hz.

  pwstat.Nuclear_Dipole_x_iy_Z = zeroDiag(...
    pwstat.Nuclear_Dipole_perpendicular.* ...
    exp(1i*PhiMatrix) .* ...
    sin(ThetaMatrix) .* cos(ThetaMatrix)...
    );
  if any(isnan(pwstat.Nuclear_Dipole_x_iy_Z(:)))
    error(['Error in getPairwiseStatistics(): ', ...
      'Nuclear_Dipole_x_iy_Z contains NANs.']);
  end

%   if System.Methyl.include 
%     methylList_ = pwstat.Methyl_Data.ID(pwstat.Methyl_Data.ID>0);
%     for ispin  = methylList_
%       hydrons_ = ispin + [1 ,2, 3];
% 
%       pwstat.Nuclear_Dipole(:,ispin) ...
%         = maxabs(pwstat.Nuclear_Dipole(:,hydrons_)')';
%       
%       pwstat.Nuclear_Dipole_perpendicular(:,ispin) ...
%         = maxabs( pwstat.Nuclear_Dipole_perpendicular(:,hydrons_)')';
% 
%       for jspin  = pwstat.Methyl_Data.ID'
%         if jspin >= ispin, break; end
%         hydrons_j = jspin + [1 ,2, 3];
% 
%         pwstat.Nuclear_Dipole(jspin,ispin) ...
%           = supabs(pwstat.Nuclear_Dipole(hydrons_j,hydrons_));
% 
%         pwstat.Nuclear_Dipole_perpendicular(:,ispin) ...
%         = supabs( pwstat.Nuclear_Dipole_perpendicular(hydrons_j,hydrons_));
% 
%       end
% 
%     end

    pwstat.Nuclear_Dipole = inf2nan(pwstat.Nuclear_Dipole);
    pwstat.Nuclear_Dipole_perpendicular ...
      = inf2nan(pwstat.Nuclear_Dipole_perpendicular);
  end

  % |bA_{mav}|^{1/2}.
  if stat_bools(STAT_BAMAX)
    pwstat.Amax = max( ...
      abs(pwstat.Hyperfine_perpendicular), ...
      abs(pwstat.Hyperfine_perpendicular'));

    pwstat.bAmax = sqrt(pwstat.Amax.*abs(pwstat.Nuclear_Dipole));

%     if System.Methyl.include
%       for ispin  = methylList_
%         hydrons_ = ispin + [1 ,2, 3];
% 
%         pwstat.Amax(:,ispin) = maxabs( pwstat.Amax(:,hydrons_)')';
%         pwstat.bAmax(:,ispin) = maxabs( pwstat.bAmax(:,hydrons_)')';
% 
%         for jspin  = pwstat.Methyl_Data.ID'
%           if jspin >= ispin, break; end
%           hydrons_j = jspin + [1 ,2, 3];
% 
%           pwstat.Amax(jspin,ispin) = supabs( pwstat.Amax(hydrons_j,hydrons_));
%           pwstat.bAmax(jspin,ispin) = supabs( pwstat.bAmax(hydrons_j,hydrons_));
%         end
% 
%       end
%     end

    pwstat.Amax = inf2nan(pwstat.Amax);
    pwstat.bAmax = inf2nan(pwstat.bAmax);
  end



  % modulation depth
  if stat_bools(STAT_ANALYTIC_HAHN)

    % Hyperfine to nucler dipole-dipole ratio
    cp = pwstat.DeltaHyperfine_perpendicular./...
      pwstat.Nuclear_Dipole_perpendicular;
    
    c = pwstat.DeltaHyperfine./pwstat.Nuclear_Dipole;

    pwstat.DeltaHyperfine_over_Nuclear_Dipole_p = cp;
    pwstat.DeltaHyperfine_over_Nuclear_Dipole = c;

    modDepth_p = 4*cp.^2./(1+cp.^2).^2;
    modDepth = 4*c.^2./(1+c.^2).^2;
    pwstat.Modulation_Depth_p = modDepth_p;
    pwstat.Modulation_Depth = modDepth;

    % modulation frequency
    pwstat.Frequency_Pair_p =  b_fliflopPerp.*sqrt(1+cp.^2); % Hz
    pwstat.Frequency_Pair   =  b.*sqrt(1+c.^2); % Hz
    
    TM = System.TMguess;
    
    pwstat.GaussianRMSD_p = getGaussianRMSD(...
      modDepth_p,pwstat.Frequency_Pair_p,TM);
    
    pwstat.GaussianRMSD = getGaussianRMSD(...
      modDepth,pwstat.Frequency_Pair,TM);

%     if System.Methyl.include
%       for ispin  = methylList_
%         hydrons_ = ispin + [1 ,2, 3];
% 
%         pwstat.Frequency_Pair_p(:,ispin) = maxabs( ...
%           pwstat.Frequency_Pair_p(:,hydrons_)')';
%         
%         pwstat.Frequency_Pair(:,ispin) = maxabs( ...
%           pwstat.Frequency_Pair(:,hydrons_)')';
%         
%         pwstat.Modulation_Depth_p(:,ispin) = maxabs( ...
%           pwstat.Modulation_Depth_p(:,hydrons_)')';
%         
%         pwstat.Modulation_Depth(:,ispin) = maxabs( ...
%           pwstat.Modulation_Depth(:,hydrons_)')';
%         
%         pwstat.GaussianRMSD_p(:,ispin) = maxabs( ...
%           pwstat.GaussianRMSD_p(:,hydrons_)')';
%         
%         pwstat.GaussianRMSD(:,ispin) = maxabs( ...
%           pwstat.GaussianRMSD(:,hydrons_)')';
% 
% 
%         for jspin  = pwstat.Methyl_Data.ID'
%           if jspin >= ispin, break; end
%           hydrons_j = jspin + [1 ,2, 3];
% 
%           pwstat.Frequency_Pair_p(jspin,ispin) = supabs( ...
%             pwstat.Frequency_Pair_p(hydrons_j,hydrons_));
%           pwstat.Frequency_Pair(jspin,ispin) = supabs( ...
%             pwstat.Frequency_Pair(hydrons_j,hydrons_));
%           pwstat.Modulation_Depth_p(jspin,ispin) = supabs( ...
%             pwstat.Modulation_Depth_p(hydrons_j,hydrons_));
%           pwstat.Modulation_Depth(jspin,ispin) = supabs( ...
%             pwstat.Modulation_Depth(hydrons_j,hydrons_));
%           pwstat.GaussianRMSD_p(jspin,ispin) = supabs( ...
%             pwstat.GaussianRMSD_p(hydrons_j,hydrons_));
%           pwstat.GaussianRMSD(jspin,ispin) = supabs( ...
%             pwstat.GaussianRMSD(hydrons_j,hydrons_));
%         end
% 
%       end

      pwstat.Frequency_Pair_p = inf2nan(pwstat.Frequency_Pair_p);
      pwstat.Frequency_Pair = inf2nan(pwstat.Frequency_Pair);
      pwstat.Modulation_Depth_p = inf2nan(pwstat.Modulation_Depth_p);
      pwstat.Modulation_Depth = inf2nan(pwstat.Modulation_Depth);
      pwstat.GaussianRMSD_p = inf2nan(pwstat.GaussianRMSD_p);
      pwstat.GaussianRMSD = inf2nan(pwstat.GaussianRMSD);
    end

    % delta_gm_gn
    pwstat.Same_g = Nuclei.Nuclear_g == Nuclei.Nuclear_g';


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
function whichStatistics(Method)
stat_bools(STAT_SPATIAL) = true;
stat_bools(STAT_HYPERFINE) = true;
stat_bools(STAT_NUCLEAR_DIPOLE) = true;

for ii = 1:numel(Method.Criteria)
  switch Method.Criteria{ii}

    case {'minimum-frequency','Gaussian RMSD','modulation'}
      stat_bools(STAT_ANALYTIC_HAHN) = true;

    case {'bAmax','minAmax','maxAmax','delta hyperfine'}
      stat_bools(STAT_BAMAX) = true;


  end
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end

