% ========================================================================
% Calculate statistics for each bath spin pair. 
% ========================================================================

function Nuclei = getPairwiseStatistics(System, Nuclei)

N = Nuclei.number;

Hyperfine = zeros(1,N);

Nuclei.ModulationDepth = zeros(N);
Nuclei.Hyperfine = zeros(N,1);
Nuclei.Nuclear_Dipole = zeros(N);
Nuclei.Frequency_Pair = zeros(N);
Nuclei.DeltaHyperfine = zeros(N);
Nuclei.SameType = eye(N);
Nuclei.Same_g = eye(N);

Nuclei.SpinSet = unique(Nuclei.Spin);

Nuclei.timescale.HF_A = 0;
Nuclei.timescale.HF_SzIxy = 0;
Nuclei.timescale.NucB = 0;
% loop over all nuclei
for inucleus = 1:N
  
  
  % Calculating hyperfine coupling
  
  gamma_e = -System.Electron.g*System.muB/System.hbar;
  Rn = norm(System.Electron.Coordinates - Nuclei.Coordinates(inucleus,:) );
  gamma_n = Nuclei.Nuclear_g(inucleus)*System.muN/System.hbar;
  
  R = System.Electron.Coordinates - Nuclei.Coordinates(inucleus,:);
  cosTheta2 = (R*[0;0;1]/norm(R))^2;
  
  Hyperfine(inucleus) = Nuclei.FermiContact(inucleus)*(2*pi)-(System.mu0/4/pi)*gamma_n*gamma_e*System.hbar*(1-3*cosTheta2)*Rn^-3;
  Nuclei.Hyperfine(inucleus) = Hyperfine(inucleus)/(2*pi); % Hz;
  % Calculating bath coupling
  for jnucleus = 1:inucleus-1
    
    if strcmp(Nuclei.Type{inucleus},Nuclei.Type{jnucleus})
      Nuclei.SameType(inucleus,jnucleus) = true;
      Nuclei.SameType(jnucleus,inucleus) = true;
    end
    
    if abs(Nuclei.Nuclear_g(inucleus)-Nuclei.Nuclear_g(jnucleus)) < 1e-9
      Nuclei.Same_g(inucleus,jnucleus) = true;
      Nuclei.Same_g(jnucleus,inucleus) = true;
    end
    
    % calculate dipolar coupling
    delta_r = Nuclei.Coordinates(inucleus,:) -Nuclei.Coordinates(jnucleus,:);
    cosThetaSquared = ( delta_r(3)/norm(delta_r) )^2;
    b = 0.25*(System.mu0/4/pi)*Nuclei.Nuclear_g(inucleus)*Nuclei.Nuclear_g(jnucleus)*System.muN^2; % J m^3.
    r = norm(Nuclei.Coordinates(inucleus,:) - Nuclei.Coordinates(jnucleus,:));
    r3 = r^3;
    b = -b/r3; % J.
%     b = b*(3*cosThetaSquared - 1); % J.
    b = b/(System.hbar); % rad/s.
    
    
    c = ( Hyperfine(inucleus)-Hyperfine(jnucleus) )/(4*b);
    
    mod_amp = 4*c^2/(1+c^2)^2;
    
    Nuclei.Frequency_Pair(inucleus,jnucleus) = 2*b*sqrt(1+c^2)/(2*pi);
    Nuclei.Frequency_Pair(jnucleus,inucleus)=    Nuclei.Frequency_Pair(inucleus,jnucleus);
    
    Nuclei.ModulationDepth(inucleus,jnucleus) = mod_amp;
    Nuclei.ModulationDepth(jnucleus,inucleus) = Nuclei.ModulationDepth(inucleus,jnucleus);
    
    Nuclei.Nuclear_Dipole(inucleus,jnucleus) = 4*b/(2*pi); % Hz.
    Nuclei.Nuclear_Dipole(jnucleus,inucleus) = Nuclei.Nuclear_Dipole(inucleus,jnucleus);
    
    Nuclei.DeltaHyperfine(inucleus,jnucleus) = (Nuclei.Hyperfine(inucleus) - Nuclei.Hyperfine(jnucleus)); % Hz
    Nuclei.DeltaHyperfine(jnucleus,inucleus) = -Nuclei.DeltaHyperfine(inucleus,jnucleus);
    
    Nuclei.timescale.HF_A = Nuclei.timescale.HF_A + abs(Nuclei.DeltaHyperfine(jnucleus,inucleus));
    Nuclei.timescale.NucB = Nuclei.timescale.NucB + abs(Nuclei.Nuclear_Dipole(jnucleus,inucleus));
  end
  
end

Nuclei.timescale.HF_SzIxy = 1/mean(abs(Nuclei.Hyperfine));
numPair = double(N*(N+1))/2;
Nuclei.timescale.NucB = numPair/Nuclei.timescale.NucB;
Nuclei.timescale.HF_A = numPair/Nuclei.timescale.HF_A;

end