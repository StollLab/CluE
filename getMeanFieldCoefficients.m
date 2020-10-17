% Calculate < psi | H | psi > for a given | psi >,
% where  | psi > = Kron_{n} |n>.
function [h,hTr] = getMeanFieldCoefficients(Nuclei,System)

% Only operators proportional Iz are needed.
theory = System.theory;
% useEZ       = theory(1);
useNZ       = false; %  theory(2);
useHF_SzIz  = false; %  theory(3);
% useHF_SzIxy = theory(4);
useNucA     = theory(5);
% useNucB     = theory(6);
useNucCD    = false; % theory(7);
% useNucEF    = theory(8);
useNQ       = false; %  theory(9);

nStates = max(Nuclei.nStates);
h_ = zeros(Nuclei.number,Nuclei.number,4);
h = zeros(Nuclei.number,Nuclei.number,4,nStates);
hTr = zeros(1,nStates);
% ENUM
E = 1; Z = 2; RAISE = 3; SZ = 4;


% Constants
km = System.mu0/(4*pi);
muB = System.muB;
muN = System.muN;
hbar = System.hbar;
ge = System.Electron.g;

for istate = 1:nStates
  for iSpin = 1:Nuclei.number
    
    psi_i = Nuclei.ZeemanStates(iSpin);
    I = Nuclei.Spin(iSpin);
    MI = psi_i - I - 1;
    
    % Nuclear Zeeman
    if useNZ
      
      NuclearZeeman = MI*Nuclei.Nuclear_g(iSpin)*System.magneticField*System.muN; % J.
      NuclearZeeman = NuclearZeeman/(2*pi*hbar);
      h_(iSpin,iSpin,E) = h_(iSpin,iSpin,E) + NuclearZeeman;
      
    end
    
    % Hyperfine
    if useHF_SzIz
      
      
      if strcmp(Nuclei.Type{iSpin},'e')
        muN = -muB;
      end
      
      gni = Nuclei.Nuclear_g(iSpin);
      r = Nuclei.Coordinates(iSpin,:);
      if size(r,2)==3
        r=r';
      end
      n = r/norm(r);
      nnt = n*n';
      r3 = norm(r)^3;
      dd = km*ge*muB*gni*muN/r3*(eye(3)-3*nnt);
      Hyperfine = -dd/(2*pi*hbar); % Hz.
      h_(iSpin,iSpin,SZ) = h_(iSpin,iSpin,SZ) + Hyperfine(3,3)*MI;
      
    end
    
    % Quadrupole
    
    if useNQ && Nuclei.StateMultiplicity(iSpin) > 2
      Q_ = Nuclei.Qtensor(:,:,iSpin);
      H_nuclear_quadrupole = MI*Q_(3,3)*MI ...
        + (Q_(1,1) + Q_(2,2))*(I*(I + 1)  -1/2*MI^2);
      
      h_(iSpin,iSpin,E) = h_(iSpin,iSpin,E) ...
        + H_nuclear_quadrupole;
      
    end
    
    for jSpin = 1:iSpin-1
      
      psi_j = Nuclei.ZeemanStates(jSpin);
      J = Nuclei.Spin(jSpin);
      MJ = psi_j - J -1;
      
      if useNucA || useNucCD
        
        
        muNi = System.muN;
        if strcmp(Nuclei.Type{iSpin},'e')
          muNi = -muB;
        end
        muNj = System.muN;
        if strcmp(Nuclei.Type{jSpin},'e')
          muNj = -muB;
        end
        
        gni = Nuclei.Nuclear_g(iSpin);
        gnj = Nuclei.Nuclear_g(jSpin);
        r = Nuclei.Coordinates(iSpin,:)-Nuclei.Coordinates(jSpin,:);
        if size(r,2)==3
          r=r';
        end
        n = r/norm(r);
        nnt = n*n';
        r3 = norm(r)^3;
        dd = -km*gni*gnj*muNi*muNj/r3*(eye(3)-3*nnt);
        dd = dd/(2*pi*hbar); % Hz.
        
        % nucleus-nucleus secular (A term)
        if useNucA
          
          h_(iSpin,jSpin,E) = h_(iSpin,jSpin,E) + MI*dd(3,3)*MJ;
          h_(jSpin,iSpin,E) = h_(iSpin,jSpin,E); % + MJ*dd(3,3)*MI;
          
          h_(iSpin,jSpin,Z) = h_(iSpin,jSpin,Z) + MI*dd(3,3);
          h_(jSpin,iSpin,Z) = h_(iSpin,jSpin,Z);% + MJ*dd(3,3);
        end
        
        % nucleus-nucleus dipolar C and D terms
        if useNucCD
          
          cd = 1/2*(dd(1,3) - 1i*dd(2,3));
          
          h_(iSpin,jSpin,RAISE) = h_(iSpin,jSpin,RAISE) + MI*cd;
          h_(jSpin,iSpin,RAISE) = h_(jSpin,iSpin,RAISE) + MJ*cd;
          
        end
      end
    end
  end
  
  for iSpin = 1:Nuclei.number
    h_(iSpin,iSpin,Z) = sum( h_(:,iSpin,Z));
    h_(iSpin,iSpin,RAISE) = sum( h_(:,iSpin,RAISE));
  end
  hTr(istate) = trace(h_(:,:,E));
  
h(:,:,:,istate) = h_;
end
 
end