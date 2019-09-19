function [Signal,AuxiliarySignal,Order_n_Signal] = doRestrictedCCE(System, Method, Nuclei, r0,verbose)

% initialize
Hyperfine = zeros(Nuclei.number,1);
Signal = ones(size(System.Time));
Order_n_Signal{1} = Signal;
numberNuclei=Nuclei.number;

if Method.conserveMemory
  AuxiliarySignal = 'The auxiliary signals are not saved in memory conservation mode.';
end

if verbose
  % threshold and counter for how often to display info
  verboseThreshold = 0.05*numberNuclei;
  verboseCounter = 0;
  disp('Starting restricted CCE.');
end

if ~all(Nuclei.Spin==1/2)
  warning('System has spins with I~=1/2. Restricted CCE is limited to spin-1/2 nuclei and will ignore them.');
end

% loop over all nuclei
for inucleus = 1:numberNuclei
  if (abs(Nuclei.Spin(inucleus) -0.5)>1e-6)
    continue;
  end
  
  % initialize the nth auxiliary signal 
  if ~Method.conserveMemory
    AuxiliarySignal{inucleus} = ones(size(Signal));
  end
  % Calculating hyperfine coupling
  
  gamma_e = -System.Electron.g*System.muB/System.hbar;
  Rn = norm(System.Electron.Coordinates - Nuclei.Coordinates(inucleus,:) );
  gamma_n = Nuclei.Nuclear_g(inucleus)*System.muN/System.hbar;
  cosTheta2 = cosFieldAngle(System.Electron.Coordinates,Nuclei.Coordinates(inucleus,:));
  cosTheta2 = cosTheta2*cosTheta2;
  
  Hyperfine(inucleus) = Nuclei.Hyperfine(inucleus)-(System.mu0/4/pi)*gamma_n*gamma_e*System.hbar*(1-3*cosTheta2)*Rn^-3;
  % Calculating bath coupling
  for jnucleus = 1:inucleus-1
    
    % skip over I != 1/2
    if Nuclei.Spin(inucleus)~=1/2, continue; end
    
    % Calculate inter-nuclear distance and skip if larger than threshold
    Rij = norm(Nuclei.Coordinates(inucleus,:)-Nuclei.Coordinates(jnucleus,:));
    if Rij > r0, continue; end
    
    
    % calculate dipolar coupling
    cosThetaSquared = (cosFieldAngle(Nuclei.Coordinates(inucleus,:),Nuclei.Coordinates(jnucleus,:)))^2;
    b = 0.25*(System.mu0/4/pi)*Nuclei.Nuclear_g(inucleus)*Nuclei.Nuclear_g(jnucleus)*System.muN^2; % J m^3.
    r = norm(Nuclei.Coordinates(inucleus,:) - Nuclei.Coordinates(jnucleus,:));
    r3 = r^3;
    b = -b*(3*cosThetaSquared - 1)/r3; % J.
    b = b/(System.hbar); % 1/s.
    
    c = ( Hyperfine(inucleus)-Hyperfine(jnucleus) )/(4*b);
    w = 2*b*sqrt(1+c^2);
    
    %if abs(c)>1e9,  error('c is larger than 1e9!');  end
    
    % v(2t) = 1 - xi*xj*c^2/(1+c^2)^2{cos[ 2tb*sqrt(1+c^2) ] -1 }^2
    AuxiliarySignal_ = 1-Nuclei.Abundance(inucleus)*Nuclei.Abundance(jnucleus)*( (c/(1+c^2)) * (cos(w*System.Time) -1) ).^2;
    Signal = Signal.*AuxiliarySignal_;
    if ~Method.conserveMemory
      AuxiliarySignal{inucleus} =AuxiliarySignal{inucleus}.*AuxiliarySignal_;
    end
  end
  if verbose
    verboseCounter = verboseCounter + 1;
    
    if verboseCounter > verboseThreshold
      verboseCounter = 0;
      fprintf('Nuclei: %d/%d, (%s).\n',inucleus,numberNuclei,datetime);
    end
    
  end
  Order_n_Signal{2} = Signal;
end
end

function cos_theta = cosFieldAngle(ParticleCoordinates_1,ParticleCoordinates_2)
R = ParticleCoordinates_2 - ParticleCoordinates_1;
cos_theta = R*[0;0;1]/norm(R);
end