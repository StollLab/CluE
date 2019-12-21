%{
CHANGES
Hamiltonian: cell(N+1,N+1) --> array(3,3,N+1,N+1)
System
Nuclei_Spin = Nuclei.Spin;
Nuclei_NumberStates = Nuclei.NumberStates;
System_full_Sz_Hyperfine = System.full_Sz_Hyperfine;
System_nuclear_dipole_A  = System.nuclear_dipole_A;
System_nuclear_dipole_B  = System.nuclear_dipole_B;
System_nuclear_dipole_CD = System.nuclear_dipole_CD;
System_nuclear_dipole_EF = System.nuclear_dipole_EF;

theory = ...
 [electron_Zeeman,...
  nuclear_Zeeman,...
  secular_Hyperfine, System.full_Sz_Hyperfine, ...
  nuclear_dipole_A, nuclear_dipole_B,nuclear_dipole_CD,nuclear_dipole_EF, ...
  nuclear_quadrupole]
%}
function [H_alpha,H_beta] = assembleHamiltonian_gpu(state_multiplicity,tensors,SpinOp,Qtensors,SpinXiXjOp,...
  theory,zeroIndex,methyl_number)

useEZ       = theory(1);
useNZ       = theory(2);
useHF_A     = theory(3);
useHF_SxIxy = theory(4);
useNucA     = theory(5);
useNucB     = theory(6);
useNucCD    = theory(7);
useNucEF    = theory(8);
useNQ       = theory(9);

clustersize = numel(state_multiplicity);

%Cluster = sort(unique(Cluster));

%if methyl_number==0 && abs(double(zeroIndex) + 1 - double(Cluster(1)) )>=1
%  error('Cluster reference failure.');
%end

%if clusterSize ~= length(Cluster)
%  error('Cluster reference failure.');
%end

Hnuc = 0;
Hhf = 0;

% iSpin is the index of the nuclear spin in the cluster
for iSpin = 1:clustersize
  
  % the ith nuclear spin in the input Hamiltonian
  % inucleus = ispin - zeroIndex + 1; % since the electron is given position 1
  
  %------------------------------------------------------------------------
  % One Nucleus Spin Hamiltonian
  %------------------------------------------------------------------------
  % Calculate nuclear quadrupole Hamiltonian
  if useNQ && state_multiplicity(iSpin) > 2
    Q_ = Qtensors(:,:,iSpin);
    [xx,xy,xz,yx,yy,yz,zx,zy,zz] = spinopidx_nq(iSpin);
    H_nuclear_quadrupole = ...
      Q_(1,1)*SpinXiXjOp(:,:,xx) + ...
      Q_(1,2)*SpinXiXjOp(:,:,xy) + ...
      Q_(1,3)*SpinXiXjOp(:,:,xz) + ...
      Q_(2,1)*SpinXiXjOp(:,:,yx) + ...
      Q_(2,2)*SpinXiXjOp(:,:,yy) + ...
      Q_(2,3)*SpinXiXjOp(:,:,yz) + ...
      Q_(3,1)*SpinXiXjOp(:,:,zx) + ...
      Q_(3,2)*SpinXiXjOp(:,:,zy) + ...
      Q_(3,3)*SpinXiXjOp(:,:,zz);
  else
    H_nuclear_quadrupole = 0;
  end
  
  % Calculate nuclear Zeeman and hyperfine Hamiltonians
  [z,r,l] = spinopidx(clusterSize,iSpin);
  Iz = SpinOp(:,:,z);
  Ix = (SpinOp(:,:,r) + SpinOp(:,:,l) )/2;
  Iy = (SpinOp(:,:,r) - SpinOp(:,:,l) )/2i;
  
  if useNZ
    H_nuclear_Zeeman_Iz = -tensors(3,3,iSpin+1,iSpin+1)*Iz;
  else
    H_nuclear_Zeeman_Iz = 0;
  end
  A = tensors(:,:,1,iSpin+1) + tensors(:,:,iSpin+1,1);
  if useHF_A
    H_hyperfine_SzIz = -A(3,3)*Iz;
  else
    H_hyperfine_SzIz = 0;
  end
  if useHF_SxIxy
    H_hyperfine_SzIx = -A(1,3)*Ix;
    H_hyperfine_SzIy = -A(2,3)*Iy;
  else
    H_hyperfine_SzIx = 0;
    H_hyperfine_SzIy = 0;
  end
  
  % Assemble single-nucleus terms in nuclear Hamiltonian
  Hnuc = Hnuc + H_nuclear_Zeeman_Iz + H_nuclear_quadrupole;
  Hhf = Hhf + H_hyperfine_SzIz + H_hyperfine_SzIx + H_hyperfine_SzIy;

  
  %------------------------------------------------------------------------
  % Loop over all nuclei with index greater than the ith nucleus.
  %------------------------------------------------------------------------
  for jSpin = iSpin+1:clustersize
        
    % the jth nuclear spin in the input Hamiltonian
    %jnucleus = jspin - zeroIndex + 1; % since the electron is given position 1
    
    % get dipolar coupling tensor
    dd = tensors(:,:,iSpin+1,jSpin+1) + tensors(:,:,jSpin+1,iSpin+1);
    
    [zz,rl,lr,zr,zl,rz,lz,rr,ll] = spinopidx2(clusterSize,iSpin,jSpin);
    
    if useNucA
      IzJz = SpinOp(:,:,zz);
      Hnn_A = dd(3,3)*IzJz;
    else
      Hnn_A = 0;
    end
    
    if useNucB
      IrJl = SpinOp(:,:,rl);
      IlJr = SpinOp(:,:,lr);
      Hnn_B = 0.25*(dd(1,1) + dd(2,2))*(IrJl+IlJr);
    else
      Hnn_B = 0;
    end
    
    if useNucCD
      IzJr = SpinOp(:,:,zr);
      IrJz = SpinOp(:,:,rz);
      IzJl = SpinOp(:,:,zl);
      IlJz = SpinOp(:,:,lz);
      cd = 1/2*(dd(1,3) - 1i*dd(2,3));
      Hnn_CD = cd*(IzJr+IrJz) + cd'*(IzJl+IlJz);
    else
      Hnn_CD = 0;
    end
    
    if useNucEF
      IrJr = SpinOp(:,:,rr);
      IlJl = SpinOp(:,:,ll);
      ef = 1/4*(dd(1,1) - dd(2,2) - 1i*(dd(1,2) + dd(2,1)));
      Hnn_EF = ef*IrJr + ef'*IlJl;
    else
      Hnn_EF = 0;
    end
    
    Hnuc = Hnuc + Hnn_A + Hnn_B + Hnn_CD + Hnn_EF;
    
  end
  
end

% Calculate electron Zeeman Hamiltonian
if useEZ
  eZeeman = tensors(:,:,1,1); % Hz
  HEZ = eZeeman(3,3)*eye(size(Hnuc)); % Hz
else
  HEZ = 0;
end

% Calculate total nuclear Hamiltonians for alpha and beta electron manifolds
H_alpha = +1/2*(HEZ + Hhf) + Hnuc;
H_beta  = -1/2*(HEZ + Hhf) + Hnuc;

% Check Hermitianity
threshold = 1e-12;
isHermA = isHermitian(H_alpha,threshold);
isHermB = isHermitian(H_beta,threshold);
if ~isHermA || ~isHermB
  error('Cluster Hamiltonian is not Hermitian.');
end

% Hermitianize
H_alpha = (H_alpha+H_alpha')/2;
H_beta = (H_beta+H_beta')/2;

end

function [ishermitian,nonHermiticity] = isHermitian(H,threshold)
mma = @(A)max(max(abs(A)));
nonHermiticity = mma(H-H')/mma(H);
ishermitian = nonHermiticity<=threshold;
end

%{
--------------------------------------------------------------------------
ENUM 1-Clusters: E O
--------------------------------------------------------------------------
E (1)  : E
O (2-4): z + -
--------------------------------------------------------------------------
2-Clusters: EE EO OE OO
--------------------------------------------------------------------------
EE (1)   : EE
EO (2-4) : Ez E+ E-
OE (5-7) : zE +E -E
OO (8-10): zz +- -+
--------------------------------------------------------------------------
ENUM 3-Clusters: EEE EEO EOE OEE EOO OEO OOE
--------------------------------------------------------------------------
EEE (1)    : EEE
EEO (2-4)  : EEz EE+ EE-
EOE (5-7)  : EzE E+E E-E
OEE (8-10) : zEE +EE -EE
EOO (11-13): Ezz E+- E-+
OEO (14-16): zEz +E- -E+
OOE (17-19): zzE +-E -+E
--------------------------------------------------------------------------
ENUM 4-Clusters: EEEE EEEO EEOE EOEE OEEE EEOO EOEO OEEO EOOE OEOE OOEE
--------------------------------------------------------------------------
EEEE (1)    :  EEEE
EEEO (2-4)  :  EEEz EEE+ EEE-
EEOE (5-7)  :  EEzE EE+E EE-E
EOEE (8-10) :  EzEE E+EE E-EE
OEEE (11-13):  zEEE +EEE -EEE
EEOO (14-16):  EEzz EE+- EE-+
EOEO (17-19):  EzEz E+E- -E+E
OEEO (20-22):  zEEz +EE- -EE+
EOOE (23-25):  EzzE E+-E E-+E
OEOE (26-28):  zEzE +E-E -E+E
OOEE (29-31):  zzEE +-EE -+EE
--------------------------------------------------------------------------
ENUM 5-Clusters
--------------------------------------------------------------------------
EEEEE (1)    :  EEEEE
EEEEO (2-4)  :  EEEEz EEEE+ EEEE-
EEEOE (5-7)  :  EEEzE EEE+E EEE-E
EEOEE (8-10) :  EEzEE EE+EE EE-EE
EOEEE (11-13):  EzEEE E+EEE E-EEE
OEEEE (14-16):  zEEEE +EEEE -EEEE
EEEOO (17-19):  EEEzz EEE+- EEE-+
EEOEO (20-22):  EEzEz EE+E- EE-E+
EOEEO (23-25):  EzEEz E+EE- E-EE+
OEEEO (26-28):  zEEEz +EEE- -EEE+
EEOOE (29-31):  EEzzE EE+-E EE-+E
EOEOE (32-34):  EzEzE E+E-E E-E+E
EOEOE (35-37):  zEEzE +EE-E -EE+E
EOOEE (38-40):  EzzEE E+-EE E-+EE
OEOEE (41-43):  zEzEE +E-EE -E+EE
OOEEE (44-46):  zzEEE +-EEE -+EEE
--------------------------------------------------------------------------
ENUM 6-Clusters
--------------------------------------------------------------------------
EEEEEE (1)    :  EEEEEE
EEEEEO (2-4)  :  EEEEEz EEEEE+ EEEEE-
EEEEOE (5-7)  :  EEEEzE EEEE+E EEEE-E
EEEOEE (8-10) :  EEEzEE EEE+EE EEE-EE
EEOEEE (11-13):  EEzEEE EE+EEE EE-EEE
EOEEEE (14-16):  EzEEEE +EEEEE -EEEEE
OEEEEE (17-19):  zEEEEE +EEEEE -EEEEE
EEEEOO (20-22):  EEEEzz EEEE+- EEEE-+
EEEOEO (23-25):  EEEzEz EEE+E- EEE-E+
EEOEEO (26-28):  EEzEEz EE+EE- EE-EE+
EOEEEO (29-31):  EzEEEz E+EEE- E-EEE+
OEEEEO (32-34):  zEEEEz +EEEE- -EEEE+
EEEOOE (35-37):  EEEzzE EEE+-E EEE-+E
EEOEOE (38-40):  EEzEzE EE+E-E EE-E+E
EOEEOE (41-43):  EzEEzE E+EE-E E-EE+E
OEEEOE (44-46):  zEEEzE +EEE-E -EEE+E
EEOOEE (47-49):  EEzzEE EE+-EE EE-+EE
EOEOEE (50-51):  EzEzEE E+E-EE E-E+EE
OEEOEE (53-55):  zEEzEE +EE-EE -EE+EE
EOOEEE (56-58):  EzzEEE E+-EEE E-+EEE
OEOEEE (59-61):  zEzEEE +E-EEE -E+EEE
OOEEEE (62-64):  zzEEEE +-EEEE -+EEEE
%}

function [z,r,l] = spinopidx(clusterSize,iSpin)

switch clusterSize
  case 1
    E = 1;
    Z = 2; R = 3; L = 4;
    zrl = [Z R L];
  case 2
    EE = 1;
    EZ = 2; ER = 3; EL =  4;
    ZE = 5; RE = 6; LE =  7;
    switch iSpin
      case 1, zrl = [ZE RE LE];
      case 2, zrl = [EZ ER EL];
    end
  case 3
    EEE =  1;
    EEZ =  2; EER =  3; EEL =  4;
    EZE =  5; ERE =  6; ELE =  7;
    ZEE =  8; REE =  9; LEE = 10;
    switch iSpin
      case 3, zrl = [EEZ EER EEL];
      case 2, zrl = [EZE ERE ELE];
      case 1, zrl = [ZEE REE LEE];
    end
  case 4
    EEEE =  1;
    EEEZ =  2; EEER =  3; EEEL =  4;
    EEZE =  5; EERE =  6; EELE =  7;
    EZEE =  8; EREE =  9; ELEE = 10;
    ZEEE = 11; REEE = 12; LEEE = 13;
    switch iSpin
      case 4, zrl = [EEEZ EEER EEEL];
      case 3, zrl = [EEZE EERE EELE];
      case 2, zrl = [EZEE EREE ELEE];
      case 1, zrl = [ZEEE REEE LEEE];
    end
  case 5
    EEEEE =  1;
    EEEEZ =  2; EEEER =  3; EEEEL =  4;
    EEEZE =  5; EEERE =  6; EEELE =  7;
    EEZEE =  8; EEREE =  9; EELEE = 10;
    EZEEE = 11; EREEE = 12; ELEEE = 13;
    ZEEEE = 14; REEEE = 15; LEEEE = 16;
    switch iSpin
      case 5, zrl = [EEEEZ EEEER EEEEL];
      case 4, zrl = [EEEZE EEERE EEELE];
      case 3, zrl = [EEZEE EEREE EELEE];
      case 2, zrl = [EZEEE EREEE ELEEE];
      case 1, zrl = [ZEEEE REEEE LEEEE];
    end
  case 6
    EEEEEE =  1;
    EEEEEZ =  2; EEEEER =  3; EEEEEL =  4;
    EEEEZE =  5; EEEERE =  6; EEEELE =  7;
    EEEZEE =  8; EEEREE =  9; EEELEE = 10;
    EEZEEE = 11; EEREEE = 12; EELEEE = 13;
    EZEEEE = 14; EREEEE = 15; ELEEEE = 16;
    ZEEEEE = 17; REEEEE = 18; LEEEEE = 19;
    switch iSpin
      case 6, zrl = [EEEEEZ EEEEER EEEEEL];
      case 5, zrl = [EEEEZE EEEERE EEEELE];
      case 4, zrl = [EEEZEE EEEREE EEELEE];
      case 3, zrl = [EEZEEE EEREEE EELEEE];
      case 2, zrl = [EZEEEE EREEEE ELEEEE];
      case 1, zrl = [ZEEEEE REEEEE LEEEEE];
    end
end

z = zrl(1);
r = zrl(2);
l = zrl(3);
end

function [zz,rl,lr,zr,zl,rz,lz,rr,ll] = spinopidx2(clusterSize,iSpin,jSpin)
switch clusterSize
  
  case 2
    ZZ = 8; RL = 9; LR = 10;
    ZR = 11;  ZL = 12;  RZ = 13;  LZ = 14;
    RR = 15;  LL = 16;
    
    IJ = [ZZ RL LR ZR ZL RZ LZ RR LL];
    
  case 3
    EZZ = 11; ERL = 12; ELR = 13;
    ZEZ = 14; REL = 15; LER = 16;
    ZZE = 17; RLE = 18; LRE = 19;
    
    EZR = 20;  EZL = 21;  ERZ = 22;  ELZ = 23;
    ZER = 24;  ZEL = 25;  REZ = 26;  LEZ = 27;
    ZRE = 28;  ZLE = 29;  RZE = 30;  LZE = 31;
    
    ERR = 32;  ELL = 33;
    RER = 34;  LEL = 35;
    RRE = 36;  LLE = 37;
    switch iSpin
      case 1 % OOE or OEO
        switch jSpin
          case 2, IJ = [ZEZ REL LER ZER ZEL REZ LEZ RER LEL]; % OEO
          case 1, IJ = [ZZE RLE LRE ZRE ZLE RZE LZE RRE LLE]; % OOE
        end
      case 2, IJ = [EZZ ERL ELR EZR EZL ERZ ELZ ERR ELL]; % EOO
    end
    
  case 4
    EEZZ = 14; EERL = 15; EELR = 16;
    EZEZ = 17; EREL = 18; ELER = 19;
    ZEEZ = 20; REEL = 21; LEER = 22;
    EZZE = 23; ERLE = 24; ELRE = 25;
    ZEZE = 26; RELE = 27; LERE = 28;
    ZZEE = 29; RLEE = 30; LREE = 31;
    EEZR = 31; EEZL = 32; EERZ = 33; EELZ = 34;
    EZER = 35; EZEL = 36; EREZ = 37; ELEZ = 38;
    ZEER = 39; ZEEL = 40; REEZ = 41; LEEZ = 42;
    EZRE = 43; EZLE = 44; ERZE = 45; ELZE = 46;
    ZERE = 47; ZELE = 48; REZE = 49; LEZE = 50;
    ZREE = 51; ZLEE = 52; RZEE = 53; LZEE = 54;
    
    EERR = 55; EELL = 56;
    ERER = 57; ELEL = 58;
    REER = 59; LEEL = 60;
    ERRE = 61; ELLE = 62;
    RERE = 63; LELE = 64;
    RREE = 65; LLEE = 66;
    switch iSpin
      case 1
        switch jSpin
          case 4, IJ = [ZEEZ REEL LEER ZEER ZEEL REEZ LEEZ REER LEEL]; % OEEO
          case 3, IJ = [ZEZE RELE LERE ZERE ZELE REZE LEZE RERE LELE]; % OEOE
          case 2, IJ = [ZZEE RLEE LREE ZREE ZLEE RZEE LZEE RREE LLEE]; % OOEE
        end
      case 2
        switch jSpin
          case 4, IJ = [EZEZ EREL ELER EZER EZEL EREZ ELEZ ERER ELEL]; % EOEO
          case 3, IJ = [EZZE ERLE ELRE EZRE EZLE ERZE ELZE ERRE ELLE]; % EOOE
        end
      case 3, IJ = [EEZZ EERL EELR EEZR EEZL EERZ EELZ EERR EELL]; % EEOO
    end
    
  case 5
    
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
    switch iSpin
      case 1
        switch jSpin
          case 5, IJ = [ZEEEZ REEEL LEEER ZEEER ZEEEL REEEZ LEEEZ REEER LEEEL]; % OEEEO
          case 4, IJ = [ZEEZE REELE LEERE ZEERE ZEELE REEZE LEEZE REERE LEELE]; % OEEOE
          case 3, IJ = [ZEZEE RELEE LEREE ZEREE ZELEE REZEE LEZEE REREE LELEE]; % OEOEE
          case 2, IJ = [ZZEEE RLEEE LREEE ZREEE ZLEEE RZEEE LZEEE RREEE LLEEE]; % OOEEE
        end
      case 2
        switch jSpin
          case 5, IJ = [EZEEZ EREEL ELEER EZEER EZEEL EREEZ ELEEZ EREER ELEEL]; % EOEEO
          case 4, IJ = [EZEZE ERELE ELERE EZERE EZELE EREZE ELEZE ERERE ELELE]; % EOEOE
          case 3, IJ = [EZZEE ERLEE ELREE EZREE EZLEE ERZEE ELZEE ERREE ELLEE]; % EOOEE
        end
      case 3
        switch jSpin
          case 5, IJ = [EEZEZ EEREL EELER EEZER EEZEL EEREZ EELEZ EERER EELEL]; % EEOEO
          case 4, IJ = [EEZZE EERLE EELRE EEZRE EEZLE EERZE EELZE EERRE EELLE]; % EEOOE
        end
      case 4, IJ = [EEEZZ EEERL EEELR EEEZR EEEZL EEERZ EEELZ EEERR EEELL]; % EEEOO
    end
    
  case 6
    
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
    switch iSpin
      case 1
        switch jSpin
          case 6, IJ = [ZEEEEZ REEEEL LEEEER ZEEEER ZEEEEL REEEEZ LEEEEZ REEEER LEEEEL]; % OEEEEO
          case 5, IJ = [ZEEEZE REEELE LEEERE ZEEERE ZEEELE REEEZE LEEEZE REEERE LEEELE]; % OEEEOE
          case 4, IJ = [ZEEZEE REELEE LEEREE ZEEREE ZEELEE REEZEE LEEZEE REEREE LEELEE]; % OEEOEE
          case 3, IJ = [ZEZEEE RELEEE LEREEE ZEREEE ZELEEE REZEEE LEZEEE REREEE LELEEE]; % OEOEEE
          case 2, IJ = [ZZEEEE RLEEEE LREEEE ZREEEE ZLEEEE RZEEEE LZEEEE RREEEE LLEEEE]; % OOEEEE
        end
      case 2
        switch jSpin
          case 6, IJ = [EZEEEZ EREEEL ELEEER EZEEER EZEEEL EREEEZ ELEEEZ EREEER ELEEEL]; % EOEEEO
          case 5, IJ = [EZEEZE EREELE ELEERE EZEERE EZEELE EREEZE ELEEZE EREERE ELEELE]; % EOEEOE
          case 4, IJ = [EZEZEE ERELEE ELEREE EZEREE EZELEE EREEZE ELEZEE EREREE ELELEE]; % EOEOEE
          case 3, IJ = [EZZEEE ERLEEE ELREEE EZREEE EZLEEE ERZEEE ELZEEE ERREEE ELLEEE]; % EOOEEE
        end
      case 3
        switch jSpin
          case 6, IJ = [EEZEEZ EEREEL EELEER EEZEER EEZEEL EEREEZ EELEEZ EEREER EELEEL]; % EEOEEO
          case 5, IJ = [EEZEZE EERELE EELERE EEZERE EEZELE EEREZE EELEZE EERERE EELELE]; % EEOEOE
          case 4, IJ = [EEZZEE EERLEE EELREE EEZREE EEZLEE EERZEE EELZEE EERREE EELLEE]; % EEOOEE
        end
      case 4
        switch jSpin
          case 6, IJ = [EEEZEZ EEEREL EEELER EEEZER EEEZEL EEEREZ EEELEZ EEERER EEELEL]; % EEEOEO
          case 5, IJ = [EEEZZE EEERLE EEELRE EEEZRE EEEZLE EEERZE EEELZE EEERRE EEELLE]; % EEEOOE
        end
      case 5, IJ = [EEEEZZ EEEERL EEEELR EEEEZR EEEEZL EEEERZ EEEELZ EEEERR EEEELL]; % EEEEOO
    end
end
%   1  2  3  4  5  6  7  8  9
% [ZZ RL LR ZR ZL RZ LZ RR LL]
zz = IJ(1);
rl = IJ(2); lr = IJ(3);
zr = IJ(4); zl = IJ(5);
rz = IJ(6); lz = IJ(7);
rr = IJ(8); ll = IJ(9);
end

function [xx,xy,xz,yx,yy,yz,zx,zy,zz] = spinopidx_nq(iSpin)

% ENUM for spin-operator indices.
XX_ = 1;  XY_ = 4;  XZ_ = 7;
YX_ = 2;  YY_ = 5;  YZ_ = 8;
ZX_ = 3;  ZY_ = 6;  ZZ_ = 9;

off = (iSpin-1)*9;

xx = XX_ + off;
xy = XY_ + off;
xz = XZ_ + off;
yx = YX_ + off;
yy = YY_ + off;
yz = YZ_ + off;
zx = ZX_ + off;
zy = ZY_ + off;
zz = ZZ_ + off;

end
