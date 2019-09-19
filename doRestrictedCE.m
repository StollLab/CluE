function [Signal,AuxiliarySignal,Order_n_Signal] = doRestrictedCE(System, Method, Nuclei, r0, verbose)

nNuclei = Nuclei.number;

% initialize 
Hyperfine = zeros(nNuclei,1);
ContributionSum = zeros(size(System.Time));
Order_n_Signal{1} = ones(size(System.Time));


if Method.conserveMemory
  AuxiliarySignal = 'The auxiliary signals are not saved in memory conservation mode.';
else
  AuxiliarySignal{nNuclei} = zeros(size(System.Time));
end

if ~all(Nuclei.Spin==1/2)
  warning('System has spins with I~=1/2. Restricted CE is limited to spin-1/2 nuclei and will ignore them.');
end

fprintf('%d of %d nuclei are spin-1/2.',sum(Nuclei.Spin==1/2),nNuclei);

verboseN = 0.01*nNuclei;

for inucleus = 1:nNuclei
  
  if ~Method.conserveMemory
    ln_AuxiliarySignal_ = zeros(size(System.Time));
  end
  
  if verbose && inucleus > verboseN
    verboseN = verboseN + 0.01*nNuclei;
    fprintf('Calculating: %d/%d.\n', inucleus,nNuclei);
  end
  
  if Nuclei.Spin(inucleus)~=1/2, continue; end
  
  % Calculating hyperfine coupling
  gamma_e = -System.Electron.g*System.muB/System.hbar;
  
  Rn = norm(System.Electron.Coordinates - Nuclei.Coordinates(inucleus,:) );
  
  gamma_n = Nuclei.Nuclear_g(inucleus)*System.muN/System.hbar;
  
  cosTheta2 = cosFieldAngle(System.Electron.Coordinates,Nuclei.Coordinates(inucleus,:))^2;
  
  Hyperfine(inucleus) = Nuclei.Hyperfine(inucleus)- ...
    (System.mu0/4/pi)*gamma_n*gamma_e*System.hbar*(1-3*cosTheta2)*Rn^-3;
  
  % Calculating bath coupling
  for jnucleus = 1:(inucleus-1)
    if Nuclei.Spin(jnucleus)~=1/2, continue; end
    
    % Calculate inter-nuclear distance and skip if larger than threshold
    Rij = norm(Nuclei.Coordinates(inucleus,:)-Nuclei.Coordinates(jnucleus,:));
    if Rij > r0, continue; end

    
    
    
    cosThetaSquared = (cosFieldAngle(Nuclei.Coordinates(inucleus,:),Nuclei.Coordinates(jnucleus,:)))^2;
    b = 0.25*(System.mu0/4/pi)*Nuclei.Nuclear_g(inucleus)*Nuclei.Nuclear_g(jnucleus)*System.muN^2; % J m^3
    r = norm(Nuclei.Coordinates(inucleus,:) - Nuclei.Coordinates(jnucleus,:));
    b = -b*(3*cosThetaSquared - 1)/r^3; % J
    b = b/System.hbar; % J -> rad/s
    
    c = ( Hyperfine(inucleus)-Hyperfine(jnucleus) )/(4*b); % unitless
    
    omega = 2*b*sqrt(1+c^2);
    
    if abs(c)>1e9
      error('c is larger than 1e9!');
    end
    
    % find the contribution from iNuc and jNuc to to ln( Signal )
    ln_AuxiliarySignal = -Nuclei.Abundance(inucleus)*Nuclei.Abundance(jnucleus)*( (c/(1+c^2)) * (cos(omega*System.Time) -1) ).^2;
    
    % add to running sum of ln( Signal ).
    ContributionSum = ContributionSum +ln_AuxiliarySignal;

    if ~Method.conserveMemory
      % add to running sum of ln( AuxiliarySignal )
      ln_AuxiliarySignal_ = ln_AuxiliarySignal_ + ln_AuxiliarySignal;
    end
    
  end
  
  if ~Method.conserveMemory
    % exponentiate to turn AuxiliarySignal into a true auxiliary signal
    AuxiliarySignal{inucleus} = exp(ln_AuxiliarySignal_);
  end
  
end

Signal = exp(ContributionSum);
Order_n_Signal{2} = Signal;

end

function cos_theta = cosFieldAngle(ParticleCoordinates_1,ParticleCoordinates_2)
R = ParticleCoordinates_2 - ParticleCoordinates_1;
cos_theta = R*[0;0;1]/norm(R);
end
