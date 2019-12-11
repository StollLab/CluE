
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
function H_out = assembleHamiltonian_gpu(Hamiltonian,SpinOp,SpinXiXjOp,Cluster,NumberStates,System_full_Sz_Hyperfine, ms,zeroIndex,clusterSize,MethylID,methyl_number,HNQ,state_multiplicity)
switch clusterSize
  case 1
    %--------------------------------------------------------------------------
    % ENUM 1-Clusters: E O
    %--------------------------------------------------------------------------
    % E (1)  : E
    % O (2-4): z + -
    
    E = 1;
    Z = 2; R = 3; L = 4;
    
    
  case 2
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
    
    ZR = 11;  ZL = 12;  RZ = 13;  LZ = 14;
    RR = 15;  LL = 16;
    
  case 3
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

    EZR = 20;  EZL = 21;  ERZ = 22;  ELZ = 23;
    ZER = 24;  ZEL = 25;  REZ = 26;  LEZ = 27;
    ZRE = 28;  ZLE = 29;  RZE = 30;  LZE = 31;
    
    ERR = 32;  ELL = 33;
    RER = 34;  LEL = 35;
    RRE = 36;  LLE = 37;
    
    
  case 4
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
    
    EEZR = 31;  EEZL = 32;  EERZ = 33;  EELZ = 34;
    EZER = 35;  EZEL = 36;  EREZ = 37;  ELEZ = 38;
    ZEER = 39;  ZEEL = 40;  REEZ = 41;  LEEZ = 42;
    EZRE = 43;  EZLE = 44;  ERZE = 45;  ELZE = 46;
    ZERE = 47;  ZELE = 48;  REZE = 49;  LEZE = 50;
    ZREE = 51;  ZLEE = 52;  RZEE = 53;  LZEE = 54;
    
    EERR = 55;  EELL = 56;
    ERER = 57;  ELEL = 58;
    REER = 59;  LEEL = 60;
    ERRE = 61;  ELLE = 62;
    RERE = 63;  LELE = 64;
    RREE = 65;  LLEE = 66;
    
    
  case 5
    % EEEEE (1)    :  EEEEE
    % EEEEO (2-4)  :  EEEEz EEEE+ EEEE-
    % EEEOE (5-7)  :  EEEzE EEE+E EEE-E
    % EEOEE (8-10) :  EEzEE EE+EE EE-EE
    % EOEEE (11-13):  EzEEE E+EEE E-EEE
    % OEEEE (14-16):  zEEEE +EEEE -EEEE
    % EEEOO (17-19):  EEEzz EEE+- EEE-+
    % EEOEO (20-22):  EEzEz EE+E- EE-E+
    % EOEEO (23-25):  EzEEz E+EE- E-EE+
    % OEEEO (26-28):  zEEEz +EEE- -EEE+
    % EEOOE (29-31):  EEzzE EE+-E EE-+E
    % EOEOE (32-34):  EzEzE E+E-E E-E+E
    % EOEOE (35-37):  zEEzE +EE-E -EE+E
    % EOOEE (38-40):  EzzEE E+-EE E-+EE
    % OEOEE (41-43):  zEzEE +E-EE -E+EE
    % OOEEE (44-46):  zzEEE +-EEE -+EEE
    
    EEEEE =  1;
    EEEEZ =  2; EEEER =  3; EEEEL =  4;
    EEEZE =  5; EEERE =  6; EEELE =  7;
    EEZEE =  8; EEREE =  9; EELEE = 10;
    EZEEE = 11; EREEE = 12; ELEEE = 13;
    ZEEEE = 14; REEEE = 15; LEEEE = 16;
    EEEZZ = 17; EEERL = 18; EEELR = 19;
    EEZEZ = 20; EEREL = 21; EELER = 22;
    EZEEZ = 23; EREEL = 24; ELEER = 25;
    ZEEEZ = 26; REEEL = 27; LEEER = 28;
    EEZZE = 29; EERLE = 30; EELRE = 31;
    EZEZE = 32; ERELE = 33; ELERE = 34;
    ZEEZE = 35; REELE = 36; LEERE = 37;
    EZZEE = 38; ERLEE = 39; ELREE = 40;
    ZEZEE = 41; RELEE = 42; LEREE = 43;
    ZZEEE = 44; RLEEE = 45; LREEE = 46;
    
  case 6
    % EEEEEE (1)    :  EEEEEE
    % EEEEEO (2-4)  :  EEEEEz EEEEE+ EEEEE-
    % EEEEOE (5-7)  :  EEEEzE EEEE+E EEEE-E
    % EEEOEE (8-10) :  EEEzEE EEE+EE EEE-EE
    % EEOEEE (11-13):  EEzEEE EE+EEE EE-EEE
    % EOEEEE (14-16):  EzEEEE +EEEEE -EEEEE
    % OEEEEE (17-19):  zEEEEE +EEEEE -EEEEE
    % EEEEOO (20-22):  EEEEzz EEEE+- EEEE-+
    % EEEOEO (23-25):  EEEzEz EEE+E- EEE-E+
    % EEOEEO (26-28):  EEzEEz EE+EE- EE-EE+
    % EOEEEO (29-31):  EzEEEz E+EEE- E-EEE+
    % OEEEEO (32-34):  zEEEEz +EEEE- -EEEE+
    % EEEOOE (35-37):  EEEzzE EEE+-E EEE-+E
    % EEOEOE (38-40):  EEzEzE EE+E-E EE-E+E
    % EOEEOE (41-43):  EzEEzE E+EE-E E-EE+E
    % OEEEOE (44-46):  zEEEzE +EEE-E -EEE+E
    % EEOOEE (47-49):  EEzzEE EE+-EE EE-+EE
    % EOEOEE (50-51):  EzEzEE E+E-EE E-E+EE
    % OEEOEE (53-55):  zEEzEE +EE-EE -EE+EE
    % EOOEEE (56-58):  EzzEEE E+-EEE E-+EEE
    % OEOEEE (59-61):  zEzEEE +E-EEE -E+EEE
    % OOEEEE (62-64):  zzEEEE +-EEEE -+EEEE
    
    EEEEEE =  1;
    EEEEEZ =  2; EEEEER =  3; EEEEEL =  4;
    EEEEZE =  5; EEEERE =  6; EEEELE =  7;
    EEEZEE =  8; EEEREE =  9; EEELEE = 10;
    EEZEEE = 11; EEREEE = 12; EELEEE = 13;
    EZEEEE = 14; EREEEE = 15; ELEEEE = 16;
    ZEEEEE = 17; REEEEE = 18; LEEEEE = 19;
    EEEEZZ = 20; EEEERL = 21; EEEELR = 22;
    EEEZEZ = 23; EEEREL = 24; EEELER = 25;
    EEZEEZ = 26; EEREEL = 27; EELEER = 28;
    EZEEEZ = 29; EREEEL = 30; ELEEER = 31;
    ZEEEEZ = 32; REEEEL = 33; LEEEER = 34;
    EEEZZE = 35; EEERLE = 36; EEELRE = 37;
    EEZEZE = 38; EERELE = 39; EELERE = 40;
    EZEEZE = 41; EREELE = 42; ELEERE = 43;
    ZEEEZE = 44; REEELE = 45; LEEERE = 46;
    EEZZEE = 47; EERLEE = 48; EELREE = 49;
    EZEZEE = 50; ERELEE = 51; ELEREE = 52;
    ZEEZEE = 53; REELEE = 54; LEEREE = 55;
    EZZEEE = 56; ERLEEE = 57; ELREEE = 58;
    ZEZEEE = 59; RELEEE = 60; LEREEE = 61;
    ZZEEEE = 62; RLEEEE = 63; LREEEE = 64;
    
end

% ENUM
X_ = 1;  Y_ = 2; Z_ = 3;
XX_ = 4; XY_ = 5;  XZ_ = 6;
YX_ = 7; YY_ = 8;  YZ_ = 9;
ZX_ = 10; ZY_ = 11;  ZZ_ = 12;

Cluster = sort(unique(Cluster));

if methyl_number==0 && abs(double(zeroIndex) + 1 - double(Cluster(1)) )>=1
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
    H_nuclear_quadrupole = 0;
  case 2
    H_out = SpinOp(:,:,EE);
    H_nuclear_quadrupole = SpinOp(:,:,EE);
  case 3
    H_out = SpinOp(:,:,EEE);
    H_nuclear_quadrupole = 0;
  case 4
    H_out = SpinOp(:,:,EEEE);
    H_nuclear_quadrupole = 0;
  case 5
    H_out = SpinOp(:,:,EEEEE);
    H_nuclear_quadrupole = 0;
  case 6
    H_out = SpinOp(:,:,EEEEEE);
    H_nuclear_quadrupole = 0;
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
          
  if state_multiplicity(Cluster(inucleus-1)) > 2
    H_nuclear_quadrupole = ...
      HNQ(1,1,Cluster(inucleus-1))*SpinXiXjOp(:,:,XX_+pre_operators*ZZ_) + HNQ(1,2,Cluster(inucleus-1))*SpinXiXjOp(:,:,XY_+pre_operators*ZZ_) + HNQ(1,3,Cluster(inucleus-1))*SpinXiXjOp(:,:,XZ_+pre_operators*ZZ_) + ...
      HNQ(2,1,Cluster(inucleus-1))*SpinXiXjOp(:,:,YX_+pre_operators*ZZ_) + HNQ(2,2,Cluster(inucleus-1))*SpinXiXjOp(:,:,YY_+pre_operators*ZZ_) + HNQ(2,3,Cluster(inucleus-1))*SpinXiXjOp(:,:,YZ_+pre_operators*ZZ_) + ...
      HNQ(3,1,Cluster(inucleus-1))*SpinXiXjOp(:,:,ZX_+pre_operators*ZZ_) + HNQ(3,2,Cluster(inucleus-1))*SpinXiXjOp(:,:,ZY_+pre_operators*ZZ_) + HNQ(3,3,Cluster(inucleus-1))*SpinXiXjOp(:,:,ZZ_+pre_operators*ZZ_);
  else
    H_nuclear_quadrupole = 0;
  end
  
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
    case 5
      
      switch pre_operators
        case 0
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,ZEEEE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,ZEEEE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,REEEE) + SpinOp(:,:,LEEEE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,REEEE) - SpinOp(:,:,LEEEE) )/2i;
          
        case 1
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EZEEE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EZEEE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EREEE) + SpinOp(:,:,ELEEE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EREEE) - SpinOp(:,:,ELEEE) )/2i; 
          
        case 2
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EEZEE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EEZEE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EEREE) + SpinOp(:,:,EELEE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EEREE) - SpinOp(:,:,EELEE) )/2i; 
          
        case 3
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EEEZE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EEEZE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EEERE) + SpinOp(:,:,EEELE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EEERE) - SpinOp(:,:,EEELE) )/2i; 
         
        case 4
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EEEEZ);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EEEEZ);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EEEER) + SpinOp(:,:,EEEEL) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EEEER) - SpinOp(:,:,EEEEL) )/2i; 
          
      end
      
    case 6
      
      switch pre_operators
        case 0
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,ZEEEEE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,ZEEEEE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,REEEEE) + SpinOp(:,:,LEEEEE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,REEEEE) - SpinOp(:,:,LEEEEE) )/2i; 
          
        case 1
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EZEEEE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EZEEEE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EREEEE) + SpinOp(:,:,ELEEEE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EREEEE) - SpinOp(:,:,ELEEEE) )/2i; 
          
        case 2
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EEZEEE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EEZEEE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EEREEE) + SpinOp(:,:,EELEEE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EEREEE) - SpinOp(:,:,EELEEE) )/2i; 
          
        case 3
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EEEZEE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EEEZEE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EEEREE) + SpinOp(:,:,EEELEE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EEEREE) - SpinOp(:,:,EEELEE) )/2i; 
          
        case 4
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EEEEZE);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EEEEZE);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EEEERE) + SpinOp(:,:,EEEELE) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EEEERE) - SpinOp(:,:,EEEELE) )/2i; 
          
        case 5 
          nuclear_Zeeman_Sz = -Hamiltonian(3,3,inucleus,inucleus)*SpinOp(:,:,EEEEEZ);
          hyperfine_SzIz =  -Hhf(3,3)*ms*SpinOp(:,:,EEEEEZ);
          hyperfine_SzIx = -Hhf(1,3)*ms*(SpinOp(:,:,EEEEER) + SpinOp(:,:,EEEEEL) )/2;
          hyperfine_SzIy = -Hhf(2,3)*ms*(SpinOp(:,:,EEEEER) - SpinOp(:,:,EEEEEL) )/2i;
          
      end
      
  end
  
  H_out = H_out + nuclear_Zeeman_Sz + hyperfine_SzIz + 1/2*(H_nuclear_quadrupole+H_nuclear_quadrupole');
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
        if useCD
          nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                       *(SpinOp(:,:,ZR)+SpinOp(:,:,RZ)) ...
                     + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                      *(SpinOp(:,:,ZL)+SpinOp(:,:,LZ));
        end
        if useEF
          nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                        *SpinOp(:,:,RR)...
                      + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                        *SpinOp(:,:,LL);
        end
        
      case 3
        
        switch pre_operators
          case 0 % OOE OR OEO
            if post_operators==1 % OOE
              nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZZE);
              nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,RLE)+SpinOp(:,:,LRE));
              
              if useCD
                nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2))*(SpinOp(:,:,ZRE)+SpinOp(:,:,RZE));
                nuclear_CD = nuclear_CD + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))*(SpinOp(:,:,ZLE)+SpinOp(:,:,LZE));
              end
              if useEF
                nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                  *SpinOp(:,:,RRE)...
                  + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                  *SpinOp(:,:,LLE);
              end
              
            else % OEO
              nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEZ);
              nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,REL)+SpinOp(:,:,LER));
              
              if useCD
                nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2))*(SpinOp(:,:,ZER)+SpinOp(:,:,REZ));
                nuclear_CD = nuclear_CD + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))*(SpinOp(:,:,ZEL)+SpinOp(:,:,LEZ));
              end
              if useEF
                nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                  *SpinOp(:,:,RER)...
                  + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                  *SpinOp(:,:,LEL);
              end
            end
            
          case 1 % EOO
            nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZZ);
            nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,ERL)+SpinOp(:,:,ELR));
            
            if useCD
              nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2))*(SpinOp(:,:,EZR)+SpinOp(:,:,ERZ));
              nuclear_CD = nuclear_CD + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))*(SpinOp(:,:,EZL)+SpinOp(:,:,ELZ));
            end
            if useEF
              nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                *SpinOp(:,:,ERR)...
                + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                *SpinOp(:,:,ELL);
            end
            
        end
        
      case 4
        switch pre_operators
          case 0
            switch post_operators
              case 0 % OEEO
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEEZ);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,REEL)+SpinOp(:,:,LEER));
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,ZEER)+SpinOp(:,:,REEZ)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,ZEEL)+SpinOp(:,:,LEEZ));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,REER)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,LEEL);
                end
                
              case 1 % OEOE
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEZE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,RELE)+SpinOp(:,:,LERE));
                
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,ZERE)+SpinOp(:,:,REZE)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,ZELE)+SpinOp(:,:,LEZE));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,RERE)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,LELE);
                end
                
              case 2 % OOEE
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZZEE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,RLEE)+SpinOp(:,:,LREE));
                
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,ZREE)+SpinOp(:,:,RZEE)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,ZLEE)+SpinOp(:,:,LZEE));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,RREE)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,LLEE);
                end
                
            end
            
          case 1
            switch post_operators
              case 0 % EOEO
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZEZ);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EREL)+SpinOp(:,:,ELER));
           
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,EZER)+SpinOp(:,:,EREZ)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,EZEL)+SpinOp(:,:,ELEZ));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,ERER)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,ELEL);
                end
                
              case 1 % EOOE
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZZE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,ERLE)+SpinOp(:,:,ELRE));
                
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,EZRE)+SpinOp(:,:,ERZE)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,EZLE)+SpinOp(:,:,ELZE));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,ERRE)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,ELLE);
                end
            end
            
          case 2 % EEOO
            nuclear_bath = Hdd(3,3)*SpinOp(:,:,EEZZ);
            nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EERL)+SpinOp(:,:,EELR));
            
            if useCD
              nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                *(SpinOp(:,:,EEZR)+SpinOp(:,:,EERZ)) ...
                + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                *(SpinOp(:,:,EZL)+SpinOp(:,:,EELZ));
            end
            if useEF
              nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                *SpinOp(:,:,EERR)...
                + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                *SpinOp(:,:,EELL);
            end
        end
        
      case 5
        switch pre_operators
          case 0
            switch post_operators
              case 0 % OEEEO
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEEEZ);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,REEEL)+SpinOp(:,:,LEEER));
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,ZEEER)+SpinOp(:,:,REEEZ)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,ZEEEL)+SpinOp(:,:,LEEEZ));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,REEER)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,LEEEL);
                end
                
              
              case 1 % OEEOE
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEEZE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,REELE)+SpinOp(:,:,LEERE));
                
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,ZEERE)+SpinOp(:,:,REEZE)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,ZEELE)+SpinOp(:,:,LEEZE));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,REERE)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,LEELE);
                end
                
              case 2 % OEOEE
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEZEE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,RELEE)+SpinOp(:,:,LEREE));
                
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,ZEREE)+SpinOp(:,:,REZEE)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,ZELEE)+SpinOp(:,:,LEZEE));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,REREE)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,LELEE);
                end
                
                
              case 3 % OOEEE
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZZEEE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,RLEEE)+SpinOp(:,:,LREEE));
                
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,ZREEE)+SpinOp(:,:,RZEEE)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,ZLEEE)+SpinOp(:,:,LZEEE));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,RREEE)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,LLEEE);
                end
                
                
            end
          case 1
            switch post_operators
              case 0 % EOEEO
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZEEZ);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EREEL)+SpinOp(:,:,ELEER));
        
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,EZEER)+SpinOp(:,:,EREEZ)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,EZEEL)+SpinOp(:,:,ELEEZ));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,EREER)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,ELEEL);
                end
                
              
              case 1 % EOEOE
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZEZE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,ERELE)+SpinOp(:,:,ELERE));
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,EZERE)+SpinOp(:,:,EREZE)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,EZELE)+SpinOp(:,:,ELEZE));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,ERERE)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,ELELE);
                end
                
              
              case 2 % EOOEE
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZZEE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,ERLEE)+SpinOp(:,:,ELREE));
                
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,EZREE)+SpinOp(:,:,ERZEE)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,EZLEE)+SpinOp(:,:,ELZEE));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,ERREE)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,ELLEE);
                end
                
            end
            
          case 2
            switch post_operators
              case 0 % EEOEO
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EEZEZ);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EEREL)+SpinOp(:,:,EELER));
                
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,EEZER)+SpinOp(:,:,EEREZ)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,EEZEL)+SpinOp(:,:,EELEZ));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,EERER)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,EELEL);
                end
                
              case 1 % EEOOE
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EEZZE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EERLE)+SpinOp(:,:,EELRE));
            
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,EEZRE)+SpinOp(:,:,EERZE)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,EEZLE)+SpinOp(:,:,EELZE));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,EERRE)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,EELLE);
                end
                
            end
            
          case 3 % EEEOO
            
            nuclear_bath = Hdd(3,3)*SpinOp(:,:,EEEZZ);
            nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EEERL)+SpinOp(:,:,EEELR));
            
            if useCD
              nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                *(SpinOp(:,:,EEEZR)+SpinOp(:,:,EEERZ)) ...
                + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                *(SpinOp(:,:,EEEZL)+SpinOp(:,:,EEELZ));
            end
            if useEF
              nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                *SpinOp(:,:,EEERR)...
                + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                *SpinOp(:,:,EEELL);
            end
            
        end
        
        case 6
        switch pre_operators
          case 0
            switch post_operators
              case 0 % OEEEEO
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEEEEZ);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,REEEEL)+SpinOp(:,:,LEEEER));
                
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,ZEEEER)+SpinOp(:,:,REEEEZ)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,ZEEEEL)+SpinOp(:,:,LEEEEZ));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,REEEER)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,LEEEEL);
                end
                
              case 1 % OEEEOE
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEEEZE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,REEELE)+SpinOp(:,:,LEEERE));
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,ZEEERE)+SpinOp(:,:,REEEZE)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,ZEEELE)+SpinOp(:,:,LEEEZE));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,REEERE)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,LEEELE);
                end
                
                
              case 2 % OEEOEE
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEEZEE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,REELEE)+SpinOp(:,:,LEEREE));
              
                if useCD
                  nuclear_CD = 1/2*(Hdd(1,3) - 1i*Hdd(2,2)) ...
                    *(SpinOp(:,:,ZEEREE)+SpinOp(:,:,REEZEE)) ...
                    + 1/2*(Hdd(1,3) + 1i*Hdd(2,2))...
                    *(SpinOp(:,:,ZEELEE)+SpinOp(:,:,LEEZEE));
                end
                if useEF
                  nuclear_EF =  1/2*(Hdd(1,1) - Hdd(1,1) - 1i*Hdd(1,2) - 1i*Hdd(2,1))...
                    *SpinOp(:,:,REEREE)...
                    + 1/2*(Hdd(1,1) - Hdd(1,1) + 1i*Hdd(1,2) + 1i*Hdd(2,1))...
                    *SpinOp(:,:,LEELEE);
                end
                
              case 3
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZEZEEE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,RELEEE)+SpinOp(:,:,LEREEE));
              case 4
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,ZZEEEE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,RLEEEE)+SpinOp(:,:,LREEEE));
                
            end
          case 1
            switch post_operators
              case 0
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZEEEZ);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EREEEL)+SpinOp(:,:,ELEEER));
              case 1
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZEEZE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EREELE)+SpinOp(:,:,ELEERE)); 
              case 2
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZEZEE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,ERELEE)+SpinOp(:,:,ELEREE));
              case 3
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EZZEEE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,ERLEEE)+SpinOp(:,:,ELREEE));
                
            end
            
          case 2
            switch post_operators
              case 0
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EEZEEZ);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EEREEL)+SpinOp(:,:,EELEER));
              case 1
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EEZEZE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EERELE)+SpinOp(:,:,EELERE));
              case 2
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EEZZEE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EERLEE)+SpinOp(:,:,EELREE));
            end
            
          case 3
            switch post_operators
              case 0
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EEEZEZ);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EEEREL)+SpinOp(:,:,EEELER));
              case 1
                nuclear_bath = Hdd(3,3)*SpinOp(:,:,EEEZZE);
                nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EEERLE)+SpinOp(:,:,EEELRE));
            end

          case 4
            nuclear_bath = Hdd(3,3)*SpinOp(:,:,EEEEZZ);
            nuclear_flipflop = 0.25*(Hdd(1,1) + Hdd(2,2))*(SpinOp(:,:,EEEERL)+SpinOp(:,:,EEEELR));
        end
        
    end
    
    H_out = H_out + nuclear_bath + nuclear_flipflop;
    if useCD, H_out = H_out + nuclear_CD; end 
    if useEF, H_out = H_out + nuclear_EF; end
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



