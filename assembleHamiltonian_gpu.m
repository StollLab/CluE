
% CHANGES
% electron manifold = ms
% Hamiltonian: cell(N+1,N+1) --> array(3,3,N+1,N+1)
%System
%{
Nuclei_Spin = Nuclei.Spin;
Nuclei_NumberStates = Nuclei.NumberStates;
System_full_Sz_Hyperfine = System.full_Sz_Hyperfine;
System_nuclear_dipole_A  = System.nuclear_dipole_A;
System_nuclear_dipole_B  = System.nuclear_dipole_B;
System_nuclear_dipole_CD = System.nuclear_dipole_CD;
System_nuclear_dipole_EF = System.nuclear_dipole_EF;

System_theory = ...
[System.full_Sz_Hyperfine,System.fullDipoleTensor,System.nuclear_dipole_A,System.nuclear_dipole_B,System.nuclear_dipole_CD,System.nuclear_dipole_EF];
%}
function H_out = assembleHamiltonian_gpu(Hamiltonian,SpinOp,Cluster,NumberStates,System_full_Sz_Hyperfine, ms,zeroIndex,clusterSize)
%--------------------------------------------------------------------------
% ENUM 1-Clusters: E O
%--------------------------------------------------------------------------
% E (1)  : E
% O (2-4): z + -

E = 1;
Z = 2; R = 3; L = 4;
%--------------------------------------------------------------------------
% 2-Clusters: EE EO OE OO
%--------------------------------------------------------------------------
% EE (1)   : EE
% EO (2-4) : Ez E+ E-
% OE (5-7) : zE +E -E
% OO (8-10): zz +- -+

EE = 1;
EZ = 2; ER = 3; EL =  4;
ZE = 5; RE = 6; LE =  7;
ZZ = 8; RL = 9; LR = 10;
%--------------------------------------------------------------------------
% ENUM 3-Clusters: EEE EEO EOE OEE EOO OEO OOE
%--------------------------------------------------------------------------
% EEE (1)    : EEE
% EEO (2-4)  : EEz EE+ EE-
% EOE (5-7)  : EzE E+E E-E
% OEE (8-10) : zEE +EE -EE
% EOO (11-13): Ezz E+- E-+
% OEO (14-16): zEz +E- -E+
% OOE (17-19): zzE +-E -+E

EEE =  1;
EEZ =  2; EER =  3; EEL =  4;
EZE =  5; ERE =  6; ELE =  7;
ZEE =  8; REE =  9; LEE = 10;
EZZ = 11; ERL = 12; ELR = 13;
ZEZ = 14; REL = 15; LER = 16;
ZZE = 17; RLE = 18; LRE = 19;

%--------------------------------------------------------------------------
% ENUM 4-Clusters: EEEE EEEO EEOE EOEE OEEE EEOO EOEO OEEO EOOE OEOE OOEE
%--------------------------------------------------------------------------
% EEEE (1)    :  EEEE
% EEEO (2-4)  :  EEEz EEE+ EEE-
% EEOE (5-7)  :  EEzE EE+E EE-E
% EOEE (8-10) :  EzEE E+EE E-EE
% OEEE (11-13):  zEEE +EEE -EEE
% EEOO (14-16):  EEzz EE+- EE-+
% EOEO (17-19):  EzEz E+E- -E+E
% OEEO (20-22):  zEEz +EE- -EE+
% EOOE (23-25):  EzzE E+-E E-+E
% OEOE (26-28):  zEzE +E-E -E+E
% OOEE (29-31):  zzEE +-EE -+EE


EEEE =  1;
EEEZ =  2; EEER =  3; EEEL =  4;
EEZE =  5; EERE =  6; EELE =  7;
EZEE =  8; EREE =  9; ELEE = 10;
ZEEE = 11; REEE = 12; LEEE = 13;
EEZZ = 14; EERL = 15; EELR = 16;
EZEZ = 17; EREL = 18; ELER = 19;
ZEEZ = 20; REEL = 21; LEER = 22;
EZZE = 23; ERLE = 24; ELRE = 25;
ZEZE = 26; RELE = 27; LERE = 28;
ZZEE = 29; RLEE = 30; LREE = 31;


Cluster = sort(unique(Cluster));

if abs(double(zeroIndex) + 1 -double(Cluster(1)))>=1
  error('Cluster reference failure.');
end

n_cluster = length(Cluster);

if clusterSize ~= n_cluster
  error('Cluster reference failure.');
end

% initialize output cluster Hamiltonian
switch clusterSize
  case 1
    H_out = SpinOp(:,:,E);
  case 2
    H_out = SpinOp(:,:,EE);
  case 3
    H_out = SpinOp(:,:,EEE);
  case 4
    H_out = SpinOp(:,:,EEEE);
end

% electron Zeeman splitting
eZeeman = ms*Hamiltonian(:,:,1,1); % Hz

% electron Zeeman frequency
H_out = H_out*eZeeman(3,3); % Hz


% icluster is a index for the cluster
inucleus = 1;
for icluster = 1:n_cluster
  
  % the ith nuclear spin in the system
  % ispin = Cluster(icluster) ;
  
  % the ith nuclear spin in the input Hamiltonian
  % inucleus = ispin - zeroIndex + 1; % since the electron is given position 1;
  inucleus = inucleus + 1;
  
  % number of identity matrices preceding the ith nucleus
  pre_operators = icluster-1;

  %------------------------------------------------------------------------
  % One Nucleus Spin Hamiltonian
  %------------------------------------------------------------------------
  Hhf = Hamiltonian(:,:,1,inucleus)+Hamiltonian(:,:,inucleus,1);
  
  switch clusterSize
    case 1
      nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,Z);  
      hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,Z);
      hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,R) + SpinOp(:,:,L) )/2;
      hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,R) - SpinOp(:,:,L) )/2i;
    case 2
      
      switch pre_operators
        case 0
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,ZE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,ZE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,RE) + SpinOp(:,:,LE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,RE) - SpinOp(:,:,LE) )/2i;
        case 1
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EZ);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EZ);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,ER) + SpinOp(:,:,EL) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,ER) - SpinOp(:,:,EL) )/2i;
      end
      
    case 3
      
      switch pre_operators
        case 0
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,ZEE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,ZEE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,REE) + SpinOp(:,:,LEE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,REE) - SpinOp(:,:,LEE) )/2i;
        case 1
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EZE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EZE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,ERE) + SpinOp(:,:,ELE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,ERE) - SpinOp(:,:,ELE) )/2i;
        case 2
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EEZ);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EEZ);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EER) + SpinOp(:,:,EEL) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EER) - SpinOp(:,:,EEL) )/2i;
          
      end
      
    case 4
      switch pre_operators
        case 0
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,ZEEE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,ZEEE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,REEE) + SpinOp(:,:,LEEE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,REEE) - SpinOp(:,:,LEEE) )/2i;
        case 1
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EZEE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EZEE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EREE) + SpinOp(:,:,ELEE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EREE) - SpinOp(:,:,ELEE) )/2i;
        case 2
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EEZE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EEZE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EERE) + SpinOp(:,:,EELE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EERE) - SpinOp(:,:,EELE) )/2i;
        case 3
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EEEZ);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EEEZ);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EEER) + SpinOp(:,:,EEEL) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EEER) - SpinOp(:,:,EEEL) )/2i;
          
      end
      
  end
  
  H_out = H_out + nuclear_Zeeman_Sz + hyperfine_SzIz;
  if System_full_Sz_Hyperfine
    H_out = H_out + hyperfine_SzIx + hyperfine_SzIy;
  end
  
  %------------------------------------------------------------------------
  % Loop over all nuclei with index greater than the ith nucleus.
  %------------------------------------------------------------------------
  jnucleus = inucleus;
  for jcluster = icluster+1:n_cluster
    
    % the jth nuclear spin in the system
    %jspin = Cluster(jcluster);
    
    % the jth nuclear spin in the input Hamiltonian
    %jnucleus = jspin - zeroIndex + 1; % since the electron is given position 1;
    
    jnucleus = jnucleus + 1;
    % number of identity matrices following the jth nucleus
    post_operators = clusterSize - jcluster;
    % get H_ij
    Hdd = Hamiltonian(:,:,inucleus,jnucleus) + Hamiltonian(:,:,jnucleus,inucleus);
    
    switch clusterSize
      
      case 2
        
        nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZZ);
        nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,RL)+SpinOp(:,:,LR));
        
      case 3
        
        switch pre_operators
          case 0
            if post_operators==1
              nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZZE);
              nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,RLE)+SpinOp(:,:,LRE));
            else
              nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEZ);
              nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,REL)+SpinOp(:,:,LER));
            end
          case 1
            nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZZ);
            nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,ERL)+SpinOp(:,:,ELR));
            
        end
        
      case 4
        switch pre_operators
          case 0
            switch post_operators
              case 0
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEEZ);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,REEL)+SpinOp(:,:,LEER));
              case 1
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEZE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,RELE)+SpinOp(:,:,LERE));
              case 2
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZZEE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,RLEE)+SpinOp(:,:,LREE));
            end
          case 1
            switch post_operators
              case 0
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZEZ);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EREL)+SpinOp(:,:,ELER));
              case 1
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZZE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,ERLE)+SpinOp(:,:,ELRE));
            end
          case 2
            nuclear_bath = Hdd(3,3)*SpinOp(:,:,EEZZ);
            nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EERL)+SpinOp(:,:,EELR));
        end
        
    end
    H_out = H_out + nuclear_bath + nuclear_flipflop;
    
    
  end
  
  
  
end

if max(max(abs( H_out - H_out') ) )  >1e-12

  if max(max(abs(nuclear_Zeeman_Sz -nuclear_Zeeman_Sz'))) > 1e-12
    disp('The Hamiltonian component nuclear_Zeeman_Sz is not Hermitian.')
  end
  
  if max(max(abs(hyperfine_SzIz -hyperfine_SzIz'))) > 1e-12
    disp('The Hamiltonian component hyperfine_SzIz is not Hermitian.')
  end
  
  if max(max(abs(hyperfine_SzIx -hyperfine_SzIx'))) > 1e-12
    disp('The Hamiltonian component hyperfine_SzIx is not Hermitian.')
  end
  
  if max(max(abs(hyperfine_SzIy -hyperfine_SzIy'))) > 1e-12
    disp('The Hamiltonian component hyperfine_SzIy is not Hermitian.')
  end
  
  if max(max(abs( nuclear_bath - nuclear_bath' ))) > 1e-12
    disp('The Hamiltonian component nuclear_bath is not Hermitian.')
  end
  
  if max(max(abs( nuclear_flipflop - nuclear_flipflop' ))) > 1e-12
    disp('The Hamiltonian component nuclear_flipflop is not Hermitian.')
  end
  
  error('Cluster Hamiltonian is not Hermitian.');
end

if max(max(abs( H_out - H_out') ) )>1e-12
  error('Cluster Hamiltonian is not Hermitian.');
end

end



function Iz = spinZ(spin)
Iz = eye(2*spin+1);
for ii = 1:(2*spin+1)
  Iz(ii,ii) = spin+ 1 - ii;
end
end

function Iminus = spinLower(spin)
Iminus = zeros(2*spin+1);
for ii = 1:(2*spin+1)-1
  Iminus(ii+1,ii) = sqrt(spin*(spin+1)-(spin - ii)*(spin+ 1 - ii));
end
end

function Iplus = spinRaise(spin)
Iplus = zeros(2*spin+1);
for ii = 1:(2*spin+1)-1
  Iplus(ii,ii+1) = sqrt(spin*(spin+1)-(spin - ii)*(spin+ 1 - ii));
end
end

function Ix = spinX(spin)
Ix = (spinRaise(spin) + spinLower(spin))/2;
end

function Iy = spinY(spin)
Iy = 1i*(-spinRaise(spin) + spinLower(spin))/2;
end



