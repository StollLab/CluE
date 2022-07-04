function [H_alpha,H_beta] = assembleHamiltonian(state_multiplicity,...
  tensors,SpinOp,Qtensors,SpinXiXjOp,...
  theory,isMethyl,methylMethod, methylTunnelingSplitting, methylID)



% debug test
% test_spinopidx();
% test_spinopidx2();

useEZ       = theory(1);
useNZ       = theory(2);
useHF_SzIz  = theory(3);
useHF_SzIxy = theory(4);
useNucA     = theory(5);
useNucB     = theory(6);
useNucCD    = theory(7);
useNucEF    = theory(8);
useNQ       = theory(9);

% Methyl version 2
useMethyl = methylMethod>0;

clusterSize = numel(state_multiplicity);

%Cluster = sort(unique(Cluster));

%if methyl_number==0 && abs(double(zeroIndex) + 1 - double(Cluster(1)) )>=1
%  error('Cluster reference failure.');
%end

%if clusterSize ~= length(Cluster)
%  error('Cluster reference failure.');
%end
% I0 = SpinOp(:,:,1);
Hnuc = 0;
Hhf = 0;

% E = 1; Z = 2; RAISE = 3; SZ = 4;
Hmf = 0;
Hmf0 = 0; 
Hmf0_ = 0;


% iSpin is the index of the nuclear spin in the cluster
for iSpin = 1:clusterSize
  
  % the ith nuclear spin in the input Hamiltonian
  % inucleus = ispin - zeroIndex + 1; % since the electron is given position 1
  
  %------------------------------------------------------------------------
  % One Nucleus Spin Hamiltonian
  %------------------------------------------------------------------------

  % Calculate nuclear Zeeman Hamiltonian
  [ze,pe,me] = spinopidx(clusterSize,iSpin);
  Iz = SpinOp(:,:,ze);
  Ix = (SpinOp(:,:,pe) + SpinOp(:,:,me) )/2;
  Iy = (SpinOp(:,:,pe) - SpinOp(:,:,me) )/2i;
  
  if useNZ
    H_nuclear_Zeeman = tensors(3,3,iSpin+1,iSpin+1)*Iz ...
    + tensors(1,1,iSpin+1,iSpin+1)*Ix + tensors(2,2,iSpin+1,iSpin+1)*Iy;% ...
    %+ tensors(2,1,iSpin+1,iSpin+1)*Ix + tensors(3,1,iSpin+1,iSpin+1)*Ix ...
    %+ tensors(1,2,iSpin+1,iSpin+1)*Iy + tensors(3,2,iSpin+1,iSpin+1)*Iy ...
    %+ tensors(1,3,iSpin+1,iSpin+1)*Iz + tensors(2,3,iSpin+1,iSpin+1)*Iz;
  else
    H_nuclear_Zeeman = 0;
  end
  
  % Calculate hyperfine Hamiltonian
  A = tensors(:,:,1,iSpin+1) + tensors(:,:,iSpin+1,1);
  if useHF_SzIz
    H_hyperfine_SzIz = A(3,3)*Iz;
  else
    H_hyperfine_SzIz = 0;
  end
  if useHF_SzIxy
    H_hyperfine_SzIx = A(3,1)*Ix;
    H_hyperfine_SzIy = A(3,2)*Iy;
  else
    H_hyperfine_SzIx = 0;
    H_hyperfine_SzIy = 0;
  end
  
  % Calculate nuclear quadrupole Hamiltonian
  if useNQ && state_multiplicity(iSpin) > 2
    Q_ = Qtensors(:,:,iSpin);
    [xx,xy,xz,yx,yy,yz,zx,zy,zz] = spinopidx_nq(iSpin,SpinXiXjOp);
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
  
  
  
  % Assemble single-nucleus terms in nuclear Hamiltonian
  Hnuc = Hnuc + H_nuclear_Zeeman + H_nuclear_quadrupole;
  Hhf = Hhf + H_hyperfine_SzIz + H_hyperfine_SzIx + H_hyperfine_SzIy;

  
  %------------------------------------------------------------------------
  % Loop over all nuclei with index greater than the ith nucleus.
  %------------------------------------------------------------------------
  for jSpin = iSpin+1:clusterSize
        
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
    
    % nucleus-nucleus pseudo-secular (B term)      
    if useNucB
      IpJm = SpinOp(:,:,rl);
      ImJp = SpinOp(:,:,lr);
      Hnn_B = 0.25*(dd(1,1) + dd(2,2))*(IpJm+ImJp);
    else
      Hnn_B = 0;
    end
    
    % nucleus-nucleus dipolar C and D terms
    if useNucCD
      IzJp = SpinOp(:,:,zr);
      IpJz = SpinOp(:,:,rz);
      IzJm = SpinOp(:,:,zl);
      ImJz = SpinOp(:,:,lz);
      cd = 1/2*(dd(1,3) - 1i*dd(2,3));
      Hnn_CD = cd*(IzJp+IpJz) + cd'*(IzJm+ImJz);
    else
      Hnn_CD = 0;
    end
    
    % nucleus-nucleus dipolar E and F terms
    if useNucEF
      IpJp = SpinOp(:,:,rr);
      ImJm = SpinOp(:,:,ll);
      ef = 1/4*(dd(1,1) - dd(2,2) - 1i*(dd(1,2) + dd(2,1)));
      Hnn_EF = ef*IpJp + ef'*ImJm;
    else
      Hnn_EF = 0;
    end
    
    if useMethyl && isMethyl(iSpin) && isMethyl(jSpin) ...
        && methylID(iSpin) == methylID(jSpin) %&& sum(methylID==methylID(iSpin))==3
      IzJz = SpinOp(:,:,zz);
      IpJm = SpinOp(:,:,rl);
      ImJp = SpinOp(:,:,lr);
      Jmethyl = -2/3*full( methylTunnelingSplitting(iSpin,jSpin) );
      H_methyl = 0;
      %H_methyl = Jmethyl*(IzJz + 0.5*IpJm + 0.5*ImJp);
      if useNucA
        H_methyl = H_methyl + Jmethyl*IzJz;
      end
      if useNucB
        H_methyl = H_methyl + Jmethyl*(0.5*IpJm + 0.5*ImJp);
      end
    else
      H_methyl = 0;
    end
    
   
    
    Hmf0 = Hmf0 + Hmf0_;
    Hnuc = Hnuc + Hnn_A + Hnn_B + Hnn_CD + Hnn_EF + H_methyl;
    
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
H_alpha = +1/2*(HEZ + Hhf + Hmf) + Hnuc + Hmf0;
H_beta  = -1/2*(HEZ + Hhf + Hmf) + Hnuc + Hmf0;
  
% Check Hermitianity
threshold = 1e-12;
[isHermA,nonHermiticityA] = isHermitian(H_alpha,threshold);
[isHermB,nonHermiticityB] = isHermitian(H_beta,threshold);
if ~isHermA || ~isHermB
  hline = '--------------------------------------------------------------';
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
  fprintf('Non-Hermiticity = {%d,%d}.\n',nonHermiticityA,nonHermiticityB);
  disp(hline);
  
  [isHerm,nonHermiticity] = isHermitian(HEZ,threshold);
  fprintf('HEZ non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('pass = %d.\n',isHerm); disp(hline);
  
  [isHerm,nonHermiticity] = isHermitian(Hhf,threshold);
  fprintf('Hhf non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('pass = %d.\n',isHerm); disp(hline);
  
  [isHerm,nonHermiticity] = isHermitian(H_hyperfine_SzIz,threshold);
  fprintf('  H_hyperfine_SzIz non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);

  [isHerm,nonHermiticity] = isHermitian(H_hyperfine_SzIy,threshold);
  fprintf('  H_hyperfine_SzIy non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);
  
  [isHerm,nonHermiticity] = isHermitian(H_hyperfine_SzIx,threshold);
  fprintf('  H_hyperfine_SzIx non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);

    [isHerm,nonHermiticity] = isHermitian(Hmf0,threshold);
  fprintf('  H_meanField_1 non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);

    [isHerm,nonHermiticity] = isHermitian(Hmf,threshold);
  fprintf('  H_meanField_Sz non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('  pass = %d.\n',isHerm); disp(hline);

  [isHerm,nonHermiticity] = isHermitian(Hnuc,threshold);
  fprintf('Hnuc non-Hermiticity = %d.\n',nonHermiticity);
  fprintf('pass = %d.\n',isHerm); disp(hline);
  if clusterSize > 1
    [isHerm,nonHermiticity] = isHermitian(Hnn_A,threshold);
    fprintf('  Hnn_A non-Hermiticity = %d.\n',nonHermiticity);
    fprintf('  pass = %d.\n',isHerm); disp(hline);
    
    [isHerm,nonHermiticity] = isHermitian(Hnn_B,threshold);
    fprintf('  Hnn_B non-Hermiticity = %d.\n',nonHermiticity);
    fprintf('  pass = %d.\n',isHerm); disp(hline);
    
    [isHerm,nonHermiticity] = isHermitian(Hnn_CD,threshold);
    fprintf('  Hnn_CD non-Hermiticity = %d.\n',nonHermiticity);
    fprintf('  pass = %d.\n',isHerm); disp(hline);
    
    [isHerm,nonHermiticity] = isHermitian(Hnn_EF,threshold);
    fprintf('  Hnn_EF non-Hermiticity = %d.\n',nonHermiticity);
    fprintf('  pass = %d.\n',isHerm); disp(hline);
    
    [isHerm,nonHermiticity] = isHermitian(H_nuclear_quadrupole,threshold);
    fprintf('  H_nuclear_quadrupole non-Hermiticity = %d.\n',nonHermiticity);
    fprintf('  pass = %d.\n',isHerm); disp(hline);
    
    [isHerm,nonHermiticity] = isHermitian(H_methyl,threshold);
    fprintf('  H_methyl non-Hermiticity = %d.\n',nonHermiticity);
    fprintf('  pass = %d.\n',isHerm); disp(hline);
  end
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
  test_spinopidx();
  test_spinopidx2()
  test_spinops(SpinOp,clusterSize);
  error('Cluster Hamiltonian is not Hermitian.');
end

% Hermitianize
H_alpha = (H_alpha+H_alpha')/2;
H_beta = (H_beta+H_beta')/2;

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

zrl = [2,3,4] + 3*(clusterSize - iSpin);

z = zrl(1);
r = zrl(2);
l = zrl(3);

end

function [z,r,l] = spinopidx_v0(clusterSize,iSpin)
switch clusterSize
  case 1
%     E = 1;
    Z = 2; R = 3; L = 4;
    zrl = [Z R L];
  case 2
%     EE = 1;
    EZ = 2; ER = 3; EL = 4;
    ZE = 5; RE = 6; LE = 7;
    switch iSpin
      case 1, zrl = [ZE RE LE];
      case 2, zrl = [EZ ER EL];
    end
  case 3
%     EEE = 1;
    EEZ = 2; EER = 3; EEL = 4;
    EZE = 5; ERE = 6; ELE = 7;
    ZEE = 8; REE = 9; LEE = 10;
    switch iSpin
      case 3, zrl = [EEZ EER EEL];
      case 2, zrl = [EZE ERE ELE];
      case 1, zrl = [ZEE REE LEE];
    end
  case 4
%     EEEE = 1;
    EEEZ = 2; EEER = 3; EEEL = 4;
    EEZE = 5; EERE = 6; EELE = 7;
    EZEE = 8; EREE = 9; ELEE = 10;
    ZEEE = 11; REEE = 12; LEEE = 13;
    switch iSpin
      case 4, zrl = [EEEZ EEER EEEL];
      case 3, zrl = [EEZE EERE EELE];
      case 2, zrl = [EZEE EREE ELEE];
      case 1, zrl = [ZEEE REEE LEEE];
    end
  case 5
%     EEEEE = 1;
    EEEEZ = 2; EEEER = 3; EEEEL = 4;
    EEEZE = 5; EEERE = 6; EEELE = 7;
    EEZEE = 8; EEREE = 9; EELEE = 10;
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
%     EEEEEE = 1;
    EEEEEZ = 2; EEEEER = 3; EEEEEL = 4;
    EEEEZE = 5; EEEERE = 6; EEEELE = 7;
    EEEZEE = 8; EEEREE = 9; EEELEE = 10;
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

function test_spinopidx()

for clusterSize = 1:6
  for iSpin = 1:clusterSize
    [z0,r0,l0] = spinopidx_v0(clusterSize,iSpin);
    zrl0 = [z0,r0,l0];
    
    [z,r,l] = spinopidx(clusterSize,iSpin);
    zrl1 = [z,r,l];
    
    if(max( abs(zrl1-zrl0))>0)
      fprintf(['\n\nError in spinopidx(): ', ...
        'for cluster size %d, single spin indices are \n'],clusterSize);
      disp(zrl1);
      disp('when they should be');
      disp(zrl0);
      error('Error above.');
      
    end
  end
end
end


function [zz,rl,lr,zr,zl,rz,lz,rr,ll] = spinopidx2(clusterSize,iSpin,jSpin)


iParticle = clusterSize + 1 -iSpin;
jParticle = clusterSize + 1 -jSpin;
% ZZ RL LR ------------------
offset = 3*clusterSize -1;
m = iParticle-jParticle;
n = jParticle-1;

zzIncrementation = m + clusterSize*n -(n^2 + n)/2;

zz = 3*zzIncrementation +  offset;
rl = zz + 1;
lr = rl + 1;

% ZR ZL RZ LZ------------------
m = clusterSize - (clusterSize - 1);
n = clusterSize - 2;
zzOffset = m + clusterSize*n -(n^2 + n)/2;

offset = offset + 3*zzOffset + 3;

zr = 4*( zzIncrementation -1) + offset;
zl = zr + 1;
rz = zl + 1;
lz = rz + 1;

% RR LL ------------------
offset = offset + 4*zzOffset;
rr = 2*( zzIncrementation -1) + offset;
ll = rr + 1;
end

function [zz,rl,lr,zr,zl,rz,lz,rr,ll] = spinopidx2_v0(clusterSize,iSpin,jSpin)

switch clusterSize
  
  case 2    
    ZZ = 8; RL = 9; LR = 10;
    ZR = 11; ZL = 12; RZ = 13;
    LZ = 14; RR = 15; LL = 16;
    
    IJ = [ZZ RL LR ZR ZL RZ LZ RR LL];
    
  case 3

    % ZZ RL LR ------------------
    EZZ = 11; ERL = 12; ELR = 13;
    ZEZ = 14; REL = 15; LER = 16;
    ZZE = 17; RLE = 18; LRE = 19;
    % ZR ZL RZ LZ------------------
    EZR = 20; EZL = 21; ERZ = 22; ELZ = 23; 
    ZER = 24; ZEL = 25; REZ = 26; LEZ = 27; 
    ZRE = 28; ZLE = 29; RZE = 30; LZE = 31;
    % RR LL ------------------
    ERR = 32; ELL = 33; 
    RER = 34; LEL = 35; 
    RRE = 36; LLE = 37;
    
    switch iSpin
      case 1 % OOE or OEO
        switch jSpin
          case 3, IJ = [ZEZ REL LER ZER ZEL REZ LEZ RER LEL]; % OEO
          case 2, IJ = [ZZE RLE LRE ZRE ZLE RZE LZE RRE LLE]; % OOE
        end
      case 2, IJ = [EZZ ERL ELR EZR EZL ERZ ELZ ERR ELL]; % EOO
    end
    
  case 4
    
    % ZZ RL LR ------------------
    EEZZ = 14; EERL = 15; EELR = 16;
    EZEZ = 17; EREL = 18; ELER = 19;
    ZEEZ = 20; REEL = 21; LEER = 22;
    EZZE = 23; ERLE = 24; ELRE = 25;
    ZEZE = 26; RELE = 27; LERE = 28;
    ZZEE = 29; RLEE = 30; LREE = 31; 
    
    % ZR ZL RZ LZ------------------
    EEZR = 32; EEZL = 33; EERZ = 34; EELZ = 35; 
    EZER = 36; EZEL = 37; EREZ = 38; ELEZ = 39; 
    ZEER = 40; ZEEL = 41; REEZ = 42; LEEZ = 43;
    EZRE = 44; EZLE = 45; ERZE = 46; ELZE = 47; 
    ZERE = 48; ZELE = 49; REZE = 50; LEZE = 51; 
    ZREE = 52; ZLEE = 53; RZEE = 54; LZEE = 55;
    
    % RR LL ------------------
    EERR = 56; EELL = 57; 
    ERER = 58; ELEL = 59; 
    REER = 60; LEEL = 61; 
    ERRE = 62; ELLE = 63; 
    RERE = 64; LELE = 65; 
    RREE = 66; LLEE = 67;
    
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
    EEEZR = 47; EEEZL = 48; EEERZ = 49;
    EEELZ = 50; EEZER = 51; EEZEL = 52;
    EEREZ = 53; EELEZ = 54; EZEER = 55;
    EZEEL = 56; EREEZ = 57; ELEEZ = 58;
    ZEEER = 59; ZEEEL = 60; REEEZ = 61;
    LEEEZ = 62; EEZRE = 63; EEZLE = 64;
    EERZE = 65; EELZE = 66; EZERE = 67;
    EZELE = 68; EREZE = 69; ELEZE = 70;
    ZEERE = 71; ZEELE = 72; REEZE = 73;
    LEEZE = 74; EZREE = 75; EZLEE = 76;
    ERZEE = 77; ELZEE = 78; ZEREE = 79;
    ZELEE = 80; REZEE = 81; LEZEE = 82;
    ZREEE = 83; ZLEEE = 84; RZEEE = 85;
    LZEEE = 86; EEERR = 87; EEELL = 88;
    EERER = 89; EELEL = 90; EREER = 91;
    ELEEL = 92; REEER = 93; LEEEL = 94;
    EERRE = 95; EELLE = 96; ERERE = 97;
    ELELE = 98; REERE = 99; LEELE = 100;
    ERREE = 101; ELLEE = 102; REREE = 103;
    LELEE = 104; RREEE = 105; LLEEE = 106;
    
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
    % ZZ RL LR ------------------
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
    
    % ZR ZL RZ LZ------------------
    EEEEZR = 65; EEEEZL = 66; EEEERZ = 67; EEEELZ = 68; 
    EEEZER = 69; EEEZEL = 70; EEEREZ = 71; EEELEZ = 72; 
    EEZEER = 73; EEZEEL = 74; EEREEZ = 75; EELEEZ = 76;
    EZEEER = 77; EZEEEL = 78; EREEEZ = 79; ELEEEZ = 80; 
    ZEEEER = 81; ZEEEEL = 82; REEEEZ = 83; LEEEEZ = 84; 
    EEEZRE = 85; EEEZLE = 86; EEERZE = 87; EEELZE = 88;
    EEZERE = 89; EEZELE = 90; EEREZE = 91; EELEZE = 92; 
    EZEERE = 93; EZEELE = 94; EREEZE = 95; ELEEZE = 96; 
    ZEEERE = 97; ZEEELE = 98; REEEZE = 99; LEEEZE = 100;
    EEZREE = 101; EEZLEE = 102; EERZEE = 103; EELZEE = 104; 
    EZEREE = 105; EZELEE = 106; EREZEE = 107; ELEZEE = 108; 
    ZEEREE = 109; ZEELEE = 110; REEZEE = 111; LEEZEE = 112;
    EZREEE = 113; EZLEEE = 114; ERZEEE = 115; ELZEEE = 116; 
    ZEREEE = 117; ZELEEE = 118; REZEEE = 119; LEZEEE = 120; 
    ZREEEE = 121; ZLEEEE = 122; RZEEEE = 123; LZEEEE = 124;
    
    % RR LL ------------------
    EEEERR = 125; EEEELL = 126; 
    EEERER = 127; EEELEL = 128; 
    EEREER = 129; EELEEL = 130;
    EREEER = 131; ELEEEL = 132; 
    REEEER = 133; LEEEEL = 134; 
    EEERRE = 135; EEELLE = 136;
    EERERE = 137; EELELE = 138; 
    EREERE = 139; ELEELE = 140; 
    REEERE = 141; LEEELE = 142;
    EERREE = 143; EELLEE = 144; 
    EREREE = 145; ELELEE = 146; 
    REEREE = 147; LEELEE = 148;
    ERREEE = 149; ELLEEE = 150; 
    REREEE = 151; LELEEE = 152; 
    RREEEE = 153; LLEEEE = 154;
    
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
          case 4, IJ = [EZEZEE ERELEE ELEREE EZEREE EZELEE EREZEE ELEZEE EREREE ELELEE]; % EOEOEE
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

function test_spinopidx2()

for clusterSize = 2:6
  for iSpin = 1:clusterSize-1
    for jSpin = iSpin+1:clusterSize
      [zz0,rl0,lr0,zr0,zl0,rz0,lz0,rr0,ll0] = ...
        spinopidx2_v0(clusterSize,iSpin,jSpin);
      zrl0 = [zz0,rl0,lr0,zr0,zl0,rz0,lz0,rr0,ll0];
      
      [zz,rl,lr,zr,zl,rz,lz,rr,ll] = spinopidx2(clusterSize,iSpin,jSpin);
      zrl1 = [zz,rl,lr,zr,zl,rz,lz,rr,ll];
      
      if(max( abs(zrl1-zrl0) )>0)
        fprintf(['\n\nError in spinopidx2(): ', ...
          'for cluster size %d, double spin indices are \n'],clusterSize);
        disp(zrl1)
        disp('when they should be');
        disp(zrl0)
        error('Error above.');
      end
    end
  end
end
end

function [xx,xy,xz,yx,yy,yz,zx,zy,zz] = spinopidx_nq(iSpin,SpinXiXjOp)

% ENUM for spin-operator indices.
XX_ = 1;  XY_ = 4;  XZ_ = 7;
YX_ = 2;  YY_ = 5;  YZ_ = 8;
ZX_ = 3;  ZY_ = 6;  ZZ_ = 9;

if size(SpinXiXjOp,3)==9
  off = 0;
else
  off = (iSpin-1)*9;
end
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


function test_spinops(SpinOp,clusterSize)


for iSpin = 1:clusterSize
  [ze,pe,me] = spinopidx(clusterSize,iSpin);
  for jSpin = iSpin+1:clusterSize


    [ez,ep,em] = spinopidx(clusterSize,jSpin);
    [zz,rl,lr,zr,zl,rz,lz,~,~] = spinopidx2(clusterSize,iSpin,jSpin);

    assert(max(max(abs( ...
      SpinOp(:,:,ze) - 2*commutator(SpinOp(:,:,rz),SpinOp(:,:,lz) ) ...
      )))==0);
    
    assert(max(max(abs( ...
      SpinOp(:,:,ez) - 2*commutator(SpinOp(:,:,zr),SpinOp(:,:,zl) ) ...
      )))==0);

    assert(max(max(abs( ...
      2*SpinOp(:,:,ze) - commutator(SpinOp(:,:,pe),SpinOp(:,:,me) ) ...
      ) ))==0);

    assert(max(max(abs( ...
      2*SpinOp(:,:,ez) - commutator(SpinOp(:,:,ep),SpinOp(:,:,em) ) ...
      ) ))==0);

    assert(max(max(abs( ...
      SpinOp(:,:,zz) -SpinOp(:,:,ze)*SpinOp(:,:,ez) ...
      )))==0);

    assert(max(max(abs( ...
      SpinOp(:,:,rl) -SpinOp(:,:,pe)*SpinOp(:,:,em) ...
      )))==0)

    assert(max(max(abs( ...
      SpinOp(:,:,lr) -SpinOp(:,:,ep)*SpinOp(:,:,me) ...
      )))==0)

    assert(max(max(abs( ...
      SpinOp(:,:,ep) -SpinOp(:,:,pe) ...
      ))) ~=0 );

    assert(max(max(abs( ...
      SpinOp(:,:,em) -SpinOp(:,:,me) ...
      ))) ~=0 );

    assert(max(max(abs( ...
      SpinOp(:,:,ez) -SpinOp(:,:,ze) ...
      ))) ~=0 );
  end
end
end
