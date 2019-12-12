function output = generateSpinOperators(spin)

% enum
E = 1;  Z = 2; RAISE = 3; LOWER = 4;
OO = [E,Z,LOWER,RAISE];

% Set spin and spin multiplicity;
multiplicity = 2*spin+1;

%--------------------------------------------------------------------------
% 1-Clusters: E O
%--------------------------------------------------------------------------
% E (1)  : E
% O (2-4): z + -

% Initialize Spin multiplicity Cluster size arrays.
maxSize = 6;
NumSpinOps = 3/2*((1:maxSize).^2 + (1:maxSize)) +1; 
SpinOp_1 = zeros(multiplicity,multiplicity,NumSpinOps(1));
SpinOp_2 = zeros(multiplicity^2,multiplicity^2,NumSpinOps(2));
SpinOp_3 = zeros(multiplicity^3,multiplicity^3,NumSpinOps(3));
SpinOp_4 = zeros(multiplicity^4,multiplicity^4,NumSpinOps(4));
SpinOp_5 = zeros(multiplicity^5,multiplicity^5,NumSpinOps(5));
SpinOp_6 = zeros(multiplicity^6,multiplicity^6,NumSpinOps(6));

% Assign single spin operators.
SpinOp_1(:,:,E) = eye(multiplicity);
SpinOp_1(:,:,Z) = spinZ(spin);
SpinOp_1(:,:,RAISE) = spinRaise(spin);
SpinOp_1(:,:,LOWER) = spinLower(spin);

%--------------------------------------------------------------------------
% 2-Clusters: EE EO OE OO
%--------------------------------------------------------------------------
% EE (1)   : EE
% EO (2-4) : Ez E+ E-
% OE (5-7) : zE +E -E
% OO (8-10): zz +- -+

% OO (11-14): z+ z- +z -z
% OO (15-16): ++ -- 

EE =  1;
EZ =  2;  ER =  3;  EL =  4;
ZE =  5;  RE =  6;  LE =  7;
ZZ =  8;  RL =  9;  LR = 10;

ZR = 11;  ZL = 12;  RZ = 13;  LZ = 14;
RR = 15;  LL = 16;

% Assign 2-cluster operators that contain an identity.
% EE

spin_index = EE;
SpinOp_2(:,:,EE) = eye(multiplicity^2);

for iop =Z:LOWER
  % EO
  spin_index = spin_index + 1;
  SpinOp_2(:,:,spin_index) = kron(SpinOp_1(:,:,E),SpinOp_1(:,:,iop));
end

for iop =Z:LOWER
  % OE
  spin_index = spin_index + 1;
  SpinOp_2(:,:,spin_index) = kron(SpinOp_1(:,:,iop),SpinOp_1(:,:,E));
end

% Assign 2-cluster operators that do not contain an identity.
%OO
SpinOp_2(:,:,ZZ) = kron(SpinOp_1(:,:,Z),SpinOp_1(:,:,Z));
SpinOp_2(:,:,RL) = kron(SpinOp_1(:,:,RAISE),SpinOp_1(:,:,LOWER));
SpinOp_2(:,:,LR) = kron(SpinOp_1(:,:,LOWER),SpinOp_1(:,:,RAISE));

SpinOp_2(:,:,ZR) = kron(SpinOp_1(:,:,Z),SpinOp_1(:,:,RAISE));
SpinOp_2(:,:,ZL) = kron(SpinOp_1(:,:,Z),SpinOp_1(:,:,LOWER));
SpinOp_2(:,:,RZ) = kron(SpinOp_1(:,:,RAISE),SpinOp_1(:,:,Z));
SpinOp_2(:,:,LZ) = kron(SpinOp_1(:,:,LOWER),SpinOp_1(:,:,Z));


SpinOp_2(:,:,RR) = kron(SpinOp_1(:,:,RAISE),SpinOp_1(:,:,RAISE));
SpinOp_2(:,:,LL) = kron(SpinOp_1(:,:,LOWER),SpinOp_1(:,:,LOWER));
%--------------------------------------------------------------------------
% 3-Clusters: EEE EEO EOE OEE EOO OEO OOE
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

% Assign 3-cluster operators that contain one non-identity operator.
SpinOp_3(:,:,E) = eye(multiplicity^3);
spin_index = E;
for iop =Z:LOWER
  % EEO
  spin_index = spin_index + 1;
  SpinOp_3(:,:,spin_index) = kron(SpinOp_2(:,:,E),SpinOp_1(:,:,iop));
end

for iop =Z:LOWER  
  % EOE
  spin_index = spin_index + 1;
  SpinOp_3(:,:,spin_index) = kron(  kron(SpinOp_1(:,:,E),SpinOp_1(:,:,iop))  ,SpinOp_1(:,:,E));
end

for iop =Z:LOWER  
  % OEE
  spin_index = spin_index + 1;
  SpinOp_3(:,:,spin_index) = kron(SpinOp_1(:,:,iop),SpinOp_2(:,:,E));
  
end

% Assign 3-cluster operators that contain two non-identity operators.
% EOO
for iop =8:10
  spin_index = spin_index + 1;
  SpinOp_3(:,:,spin_index)= kron(SpinOp_1(:,:,E),SpinOp_2(:,:,iop));
end

% OEO

spin_index = spin_index + 1;
SpinOp_3(:,:,spin_index) = kron(  kron(SpinOp_1(:,:,Z)     ,SpinOp_1(:,:,E))  ,SpinOp_1(:,:,Z));

spin_index = spin_index + 1;
SpinOp_3(:,:,spin_index) = kron(  kron(SpinOp_1(:,:,RAISE) ,SpinOp_1(:,:,E))  ,SpinOp_1(:,:,LOWER));

spin_index = spin_index + 1;
SpinOp_3(:,:,spin_index) = kron( kron(SpinOp_1(:,:,LOWER) ,SpinOp_1(:,:,E))  ,SpinOp_1(:,:,RAISE));

% OOE

for iop =8:10
  spin_index = spin_index + 1;
  SpinOp_3(:,:,spin_index)= kron(SpinOp_2(:,:,iop),SpinOp_1(:,:,E));
end

SpinOp_3(:,:,EZR) = kron(SpinOp_1(:,:,E), SpinOp_2(:,:,ZR));
SpinOp_3(:,:,EZL) = kron(SpinOp_1(:,:,E), SpinOp_2(:,:,ZL));
SpinOp_3(:,:,ERZ) = kron(SpinOp_1(:,:,E), SpinOp_2(:,:,RZ));
SpinOp_3(:,:,ELZ) = kron(SpinOp_1(:,:,E), SpinOp_2(:,:,LZ));

SpinOp_3(:,:,ZER) = kron(SpinOp_1(:,:,Z),     SpinOp_2(:,:,EZ));
SpinOp_3(:,:,ZEL) = kron(SpinOp_1(:,:,Z),     SpinOp_2(:,:,EL));
SpinOp_3(:,:,REZ) = kron(SpinOp_1(:,:,RAISE), SpinOp_2(:,:,EZ));
SpinOp_3(:,:,LEZ) = kron(SpinOp_1(:,:,LOWER), SpinOp_2(:,:,EZ));

SpinOp_3(:,:,ZRE) = kron(SpinOp_2(:,:,ZR), SpinOp_1(:,:,E));
SpinOp_3(:,:,ZLE) = kron(SpinOp_2(:,:,ZL), SpinOp_1(:,:,E));
SpinOp_3(:,:,RZE) = kron(SpinOp_2(:,:,RZ), SpinOp_1(:,:,E));
SpinOp_3(:,:,LZE) = kron(SpinOp_2(:,:,LZ), SpinOp_1(:,:,E));


SpinOp_3(:,:,ERR) = kron(SpinOp_1(:,:,E), SpinOp_2(:,:,RR));
SpinOp_3(:,:,ELL) = kron(SpinOp_1(:,:,E), SpinOp_2(:,:,LL));

SpinOp_3(:,:,RER) = kron(SpinOp_1(:,:,RAISE), SpinOp_2(:,:,ER));
SpinOp_3(:,:,LEL) = kron(SpinOp_1(:,:,LOWER), SpinOp_2(:,:,EL));

SpinOp_3(:,:,RRE) = kron(SpinOp_2(:,:,RR), SpinOp_1(:,:,E));
SpinOp_3(:,:,LLE) = kron(SpinOp_2(:,:,LL) ,SpinOp_1(:,:,E));
%--------------------------------------------------------------------------
% 4-Clusters: EEEE EEEO EEOE EOEE OEEE EEOO EOEO OEEO EOOE OEOE OOEE
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
% Assign 4-cluster operators that contain one non-identity operator.
SpinOp_4(:,:,E) = eye(multiplicity^4);
spin_index = E;
for iop =Z:LOWER
  % EEEO
  spin_index = spin_index + 1;
  SpinOp_4(:,:,spin_index) = kron(SpinOp_3(:,:,E),SpinOp_1(:,:,iop));
end

for iop =Z:LOWER
  % EEOE
  spin_index = spin_index + 1;
  SpinOp_4(:,:,spin_index) = kron(  kron(SpinOp_2(:,:,E),SpinOp_1(:,:,iop))  ,SpinOp_1(:,:,E));
end

for iop =Z:LOWER
  % EOEE
  spin_index = spin_index + 1;
  SpinOp_4(:,:,spin_index) = kron(  kron(SpinOp_1(:,:,E),SpinOp_1(:,:,iop))  ,SpinOp_2(:,:,E));
end

for iop =Z:LOWER
  % OEEE
  spin_index = spin_index + 1;
  SpinOp_4(:,:,spin_index) = kron(SpinOp_1(:,:,iop),SpinOp_3(:,:,E));
  
end

% Assign 4-cluster operators that contain two non-identity operators.
% EEOO
for iop =8:10
  spin_index = spin_index + 1;
  SpinOp_4(:,:,spin_index)= kron(SpinOp_2(:,:,E),SpinOp_2(:,:,iop));
end

% EOEO
spin_index = spin_index + 1;
SpinOp_4(:,:,spin_index) = kron(SpinOp_1(:,:,E), kron(  kron(SpinOp_1(:,:,Z)     ,SpinOp_1(:,:,E))  ,SpinOp_1(:,:,Z)));

spin_index = spin_index + 1;
SpinOp_4(:,:,spin_index) = kron(SpinOp_1(:,:,E), kron(  kron(SpinOp_1(:,:,RAISE) ,SpinOp_1(:,:,E))  ,SpinOp_1(:,:,LOWER)));

spin_index = spin_index + 1;
SpinOp_4(:,:,spin_index) = kron(SpinOp_1(:,:,E), kron( kron(SpinOp_1(:,:,LOWER) ,SpinOp_1(:,:,E))  ,SpinOp_1(:,:,RAISE)));

% OEEO
spin_index = spin_index + 1;
SpinOp_4(:,:,spin_index) = kron(  kron(SpinOp_1(:,:,Z)     ,SpinOp_2(:,:,E))  ,SpinOp_1(:,:,Z));

spin_index = spin_index + 1;
SpinOp_4(:,:,spin_index) = kron(  kron(SpinOp_1(:,:,RAISE) ,SpinOp_2(:,:,E))  ,SpinOp_1(:,:,LOWER));

spin_index = spin_index + 1;
SpinOp_4(:,:,spin_index) = kron( kron(SpinOp_1(:,:,LOWER) ,SpinOp_2(:,:,E))  ,SpinOp_1(:,:,RAISE));


% EOOE

for iop =8:10
  spin_index = spin_index + 1;
  SpinOp_4(:,:,spin_index)=  kron(SpinOp_1(:,:,E), kron(SpinOp_2(:,:,iop),SpinOp_1(:,:,E)));
end

% OEOE
spin_index = spin_index + 1;
SpinOp_4(:,:,spin_index) = kron( kron(  kron(SpinOp_1(:,:,Z)     ,SpinOp_1(:,:,E))  ,SpinOp_1(:,:,Z)),SpinOp_1(:,:,E));

spin_index = spin_index + 1;
SpinOp_4(:,:,spin_index) = kron( kron(  kron(SpinOp_1(:,:,RAISE) ,SpinOp_1(:,:,E))  ,SpinOp_1(:,:,LOWER)),SpinOp_1(:,:,E));

spin_index = spin_index + 1;
SpinOp_4(:,:,spin_index) = kron( kron( kron(SpinOp_1(:,:,LOWER) ,SpinOp_1(:,:,E))  ,SpinOp_1(:,:,RAISE)) ,SpinOp_1(:,:,E));

% OOEE

for iop =8:10
  spin_index = spin_index + 1;
  SpinOp_4(:,:,spin_index)= kron(SpinOp_2(:,:,iop),SpinOp_2(:,:,E));
end

SpinOp_4(:,:,EEZR) = kron(SpinOp_2(:,:,EE), SpinOp_2(:,:,ZR));
SpinOp_4(:,:,EEZL) = kron(SpinOp_2(:,:,EE), SpinOp_2(:,:,ZL));
SpinOp_4(:,:,EERZ) = kron(SpinOp_2(:,:,EE), SpinOp_2(:,:,RZ));
SpinOp_4(:,:,EELZ) = kron(SpinOp_2(:,:,EE), SpinOp_2(:,:,LZ));

SpinOp_4(:,:,EZER) = kron(SpinOp_2(:,:,EZ), SpinOp_2(:,:,EZ));
SpinOp_4(:,:,EZEL) = kron(SpinOp_2(:,:,EZ), SpinOp_2(:,:,EL));
SpinOp_4(:,:,EREZ) = kron(SpinOp_2(:,:,EZ), SpinOp_2(:,:,EZ));
SpinOp_4(:,:,ELEZ) = kron(SpinOp_2(:,:,EZ), SpinOp_2(:,:,EZ));

SpinOp_4(:,:,ZEER) = kron(SpinOp_2(:,:,ZE), SpinOp_2(:,:,ER));
SpinOp_4(:,:,ZEEL) = kron(SpinOp_2(:,:,ZE), SpinOp_2(:,:,EL));
SpinOp_4(:,:,REEZ) = kron(SpinOp_2(:,:,RE), SpinOp_2(:,:,EZ));
SpinOp_4(:,:,LEEZ) = kron(SpinOp_2(:,:,LE), SpinOp_2(:,:,EZ));

SpinOp_4(:,:,EZRE) = kron(SpinOp_2(:,:,EZ), SpinOp_2(:,:,RE));
SpinOp_4(:,:,EZLE) = kron(SpinOp_2(:,:,EZ), SpinOp_2(:,:,LE));
SpinOp_4(:,:,ERZE) = kron(SpinOp_2(:,:,ER), SpinOp_2(:,:,ZE));
SpinOp_4(:,:,ELZE) = kron(SpinOp_2(:,:,EL), SpinOp_2(:,:,ZE));

SpinOp_4(:,:,ZERE) = kron(SpinOp_2(:,:,ZE), SpinOp_2(:,:,RE));
SpinOp_4(:,:,ZELE) = kron(SpinOp_2(:,:,ZE), SpinOp_2(:,:,LE));
SpinOp_4(:,:,REZE) = kron(SpinOp_2(:,:,RE), SpinOp_2(:,:,ZE));
SpinOp_4(:,:,LEZE) = kron(SpinOp_2(:,:,LE), SpinOp_2(:,:,ZE));

SpinOp_4(:,:,ZREE) = kron(SpinOp_2(:,:,ZR), SpinOp_2(:,:,EE));
SpinOp_4(:,:,ZLEE) = kron(SpinOp_2(:,:,ZL), SpinOp_2(:,:,EE));
SpinOp_4(:,:,RZEE) = kron(SpinOp_2(:,:,RZ), SpinOp_2(:,:,EE));
SpinOp_4(:,:,LZEE) = kron(SpinOp_2(:,:,LZ), SpinOp_2(:,:,EE));

SpinOp_4(:,:,EERR) = kron(SpinOp_2(:,:,EE), SpinOp_2(:,:,RR));
SpinOp_4(:,:,EELL) = kron(SpinOp_2(:,:,EE), SpinOp_2(:,:,LL));

SpinOp_4(:,:,ERER) = kron(SpinOp_2(:,:,ER), SpinOp_2(:,:,ER));
SpinOp_4(:,:,ELEL) = kron(SpinOp_2(:,:,EL), SpinOp_2(:,:,EL));

SpinOp_4(:,:,REER) = kron(SpinOp_2(:,:,RE), SpinOp_2(:,:,ER));
SpinOp_4(:,:,LEEL) = kron(SpinOp_2(:,:,LE) ,SpinOp_2(:,:,EL));

SpinOp_4(:,:,ERRE) = kron(SpinOp_2(:,:,ER), SpinOp_2(:,:,RE));
SpinOp_4(:,:,ELLE) = kron(SpinOp_2(:,:,EL), SpinOp_2(:,:,LE));

SpinOp_4(:,:,RERE) = kron(SpinOp_2(:,:,RE), SpinOp_2(:,:,RE));
SpinOp_4(:,:,LELE) = kron(SpinOp_2(:,:,LE), SpinOp_2(:,:,LE));

SpinOp_4(:,:,RREE) = kron(SpinOp_2(:,:,RR), SpinOp_2(:,:,EE));
SpinOp_4(:,:,LLEE) = kron(SpinOp_2(:,:,LL) ,SpinOp_2(:,:,EE));

%--------------------------------------------------------------------------
% 5-Clusters: EEEEE EEEEO EEEOE EEOEE EOEEE OEEEE EEEOO EEOEO EOEEO OEEEO
% EEOOE EOEOE OEEOE EOOEE OEOEE OOEEE
%--------------------------------------------------------------------------
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

% EEEEE (1)    :  EEEEE
spin_index = 1;
SpinOp_5(:,:,spin_index) = eye(multiplicity^5);

% EEEEO (2-4)  :  EEEEz EEEE+ EEEE-
for iop =Z:LOWER
  % EEEO
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron(SpinOp_4(:,:,E),SpinOp_1(:,:,iop));
end

% EEEOE (5-7)  :  EEEzE EEE+E EEE-E
for iop =Z:LOWER
  % EEOE
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron(SpinOp_3(:,:,E),SpinOp_1(:,:,iop)),SpinOp_1(:,:,E));
end

% EEOEE (8-10) :  EEzEE EE+EE EE-EE
for iop =Z:LOWER
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron(SpinOp_2(:,:,E),SpinOp_1(:,:,iop)),SpinOp_2(:,:,E));
end

% EOEEE (11-13):  EzEEE E+EEE E-EEE
for iop =Z:LOWER
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron(SpinOp_1(:,:,E),SpinOp_1(:,:,iop)),SpinOp_3(:,:,E));
end

% OEEEE (14-16):  zEEEE +EEEE -EEEE
for iop =Z:LOWER
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron(SpinOp_1(:,:,iop),SpinOp_4(:,:,E));
end

% EEEOO (17-19):  EEEzz EEE+- EEE-+
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron(SpinOp_3(:,:,E),SpinOp_1(:,:,iop)),SpinOp_1(:,:,jop));
end

% EEOEO (20-22):  EEzEz EE+E- EE-E+
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron( kron(...
    SpinOp_2(:,:,E), ...
    SpinOp_1(:,:,iop)), ...
    SpinOp_1(:,:,E) ), ...
    SpinOp_1(:,:,jop) );
end

% EOEEO (23-25):  EzEEz E+EE- E-EE+
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron( kron(...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,iop)), ...
    SpinOp_2(:,:,E) ), ...
    SpinOp_1(:,:,jop) );
end

% OEEEO (26-28):  zEEEz +EEE- -EEE+ 
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron( kron(...
    SpinOp_1(:,:,iop), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_2(:,:,E) ), ...
    SpinOp_1(:,:,jop) );
end

% EEOOE (29-31):  EEzzE EE+-E EE-+E
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron( kron(...
    SpinOp_2(:,:,E), ...
    SpinOp_1(:,:,iop)), ...
    SpinOp_1(:,:,jop) ), ...
    SpinOp_1(:,:,E) );
end

% EOEOE (32-34):  EzEzE E+E-E E-E+E
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron( kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,iop)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,jop)) , ...
    SpinOp_1(:,:,E) );
end

% EOEOE (35-37):  zEEzE +EE-E -EE+E
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron( kron( kron( ...
    SpinOp_1(:,:,iop), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,jop)) , ...
    SpinOp_1(:,:,E) );
end

% EOOEE (38-40):  EzzEE E+-EE E-+EE
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron( kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,iop)), ...
    SpinOp_1(:,:,jop)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E) );
end

% OEOEE (41-43):  zEzEE +E-EE -E+EE
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron( kron( kron( ...
    SpinOp_1(:,:,iop), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,jop)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E) );
end

% OOEEE (44-46):  zzEEE +-EEE -+EEE
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_5(:,:,spin_index) = kron( kron( kron( kron( ...
    SpinOp_1(:,:,iop), ...
    SpinOp_1(:,:,jop)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E) );
end


%--------------------------------------------------------------------------
% 6-Clusters: EEEEEE EEEEEO EEEEOE EEEOEE EEOEEE EOEEEE OEEEEE EEEEOO
% EEEOEO EEOEEO EOEEEO OEEEEO EEEOOE EEOEOE EOEEOE OEEEOE EEOOEE EOEOEE OEEOEE EOOEEE OOEEEE
%--------------------------------------------------------------------------
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

% EEEEEE (1)    :  EEEEEE
spin_index = 1;
SpinOp_6(:,:,spin_index) = eye(multiplicity^6);

% EEEEEO (2-4)  :  EEEEEz EEEEE+ EEEEE-
for iop =Z:LOWER
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( ...
    SpinOp_4(:,:,E), ...
    SpinOp_1(:,:,E) ), ...
    SpinOp_1(:,:,iop) );
end

% EEEEOE (5-7)  :  EEEEzE EEEE+E EEEE-E
for iop =Z:LOWER
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( ...
    SpinOp_4(:,:,E), ...
    SpinOp_1(:,:,iop) ), ...
    SpinOp_1(:,:,E) );
end

% EEEOEE (8-10) :  EEEzEE EEE+EE EEE-EE
for iop =Z:LOWER
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( ...
    SpinOp_3(:,:,E), ...
    SpinOp_1(:,:,iop) ), ...
    SpinOp_2(:,:,E) );
end

% EEOEEE (11-13):  EEzEEE EE+EEE EE-EEE
for iop =Z:LOWER
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( ...
    SpinOp_2(:,:,E), ...
    SpinOp_1(:,:,iop) ), ...
    SpinOp_3(:,:,E) );
end

% EOEEEE (14-16):  EzEEEE +EEEEE -EEEEE
for iop =Z:LOWER
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,iop) ), ...
    SpinOp_4(:,:,E) );
end

% OEEEEE (17-19):  zEEEEE +EEEEE -EEEEE
for iop =Z:LOWER
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( ...
    SpinOp_1(:,:,iop), ...
    SpinOp_1(:,:,E) ), ...
    SpinOp_4(:,:,E) );
end

% EEEEOO (20-22):  EEEEzz EEEE+- EEEE-+
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,iop)) , ...
    SpinOp_1(:,:,jop) );
end


% EEEOEO (23-25):  EEEzEz EEE+E- EEE-E+
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,iop)) , ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,jop) );
end


% EEOEEO (26-28):  EEzEEz EE+EE- EE-EE+
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,iop)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,jop) );
end


% EOEEEO (29-31):  EzEEEz E+EEE- E-EEE+ 
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,iop)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,jop) );
end


% OEEEEO (32-34):  zEEEEz +EEEE- -EEEE+ 
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,iop), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,jop) );
end


% EEEOOE (35-37):  EEEzzE EEE+-E EEE-+E
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,iop)) , ...
    SpinOp_1(:,:,jop)) , ...
    SpinOp_1(:,:,E) );
end


% EEOEOE (38-40):  EEzEzE EE+E-E EE-E+E
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,iop)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,jop)) , ...
    SpinOp_1(:,:,E) );
end


% EOEEOE (41-43):  EzEEzE E+EE-E E-EE+E
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,iop)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,jop)) , ...
    SpinOp_1(:,:,E) );
end


% OEEEOE (44-46):  zEEEzE +EEE-E -EEE+E
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,iop), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,jop)) , ...
    SpinOp_1(:,:,E) );
end


% EEOOEE (47-49):  EEzzEE EE+-EE EE-+EE
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,iop)), ...
    SpinOp_1(:,:,jop)) , ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E) );
end


% EOEOEE (50-51):  EzEzEE E+E-EE E-E+EE
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,iop)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,jop)) , ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E) );
end


% OEEOEE (53-55):  zEEzEE +EE-EE -EE+EE
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,iop), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,jop)) , ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E) );
end


% EOOEEE (56-58):  EzzEEE E+-EEE E-+EEE
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,E), ...
    SpinOp_1(:,:,iop)), ...
    SpinOp_1(:,:,jop)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E) );
end


% OEOEEE (59-61):  zEzEEE +E-EEE -E+EEE
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,iop), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,jop)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E) );
end


% OOEEEE (62-64):  zzEEEE +-EEEE -+EEEE
for iop =Z:LOWER
  jop = OO(iop);
  spin_index = spin_index + 1;
  SpinOp_6(:,:,spin_index) = kron( kron( kron( kron( kron( ...
    SpinOp_1(:,:,iop), ...
    SpinOp_1(:,:,jop)), ...
    SpinOp_1(:,:,E)), ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E)) , ...
    SpinOp_1(:,:,E) );
end




%--------------------------------------------------------------------------
% 7-Clusters: 
%--------------------------------------------------------------------------
% EEEEEEE (1)    :  EEEEEEE
% EEEEEEO (2-4)  :  EEEEEEz EEEEEE+ EEEEEE-
% EEEEEOE (5-7)  :  EEEEEzE EEEEE+E EEEEE-E
% EEEEOEE (8-10) :  EEEEzEE EEEE+EE EEEE-EE
% EEEOEEE (11-13):  EEEzEEE EEE+EEE EEE-EEE
% EEOEEEE (14-19):  EEzEEEE EE+EEEE EE-EEEE
% EOEEEEE (17-22):  EzEEEEE E+EEEEE E-EEEEE
% OEEEEEE (20-25):  zEEEEEE +EEEEEE -EEEEEE
% EEEEEOO (20-28):  EEEEEzz EEEEE+- EEEEE-+
% EEEEOEO (23-31):  EEEEzEz EEEE+E- EEEE-E+
% EEEOEEO (26-34):  EEEzEEz EEE+EE- EEE-EE+
% EEOEEEO (29-37):  EEzEEEz EE+EEE- EE-EEE+ 
% EOEEEEO (32-40):  EzEEEEz E+EEEE- E-EEEE+ 
% OEEEEEO (32-43):  zEEEEEz +EEEEE- -EEEEE+ 
% EEEEOOE (35-46):  EEEEzzE EEEE+-E EEEE-+E
% EEEOEOE (38-49):  EEEzEzE EEE+E-E EEE-E+E
% EEOEEOE (41-52):  EEzEEzE EE+EE-E EE-EE+E
% EOEEEOE (44-55):  EzEEEzE E+EEE-E E-EEE+E
% OEEEEOE (44-58):  zEEEEzE +EEEE-E -EEEE+E
% EEEOOEE (47-61):  EEEzzEE EEE+-EE EEE-+EE
% EEOEOEE (50-64):  EEzEzEE EE+E-EE EE-E+EE
% EOEEOEE (53-67):  EzEEzEE E+EE-EE E-EE+EE
% OEEEOEE (53-70):  zEEEzEE +EEE-EE -EEE+EE
% EEOOEEE (56-73):  EEzzEEE EE+-EEE EE-+EEE
% EOEOEEE (59-76):  EzEzEEE E+E-EEE E-E+EEE
% OEEOEEE (59-79):  zEEzEEE +EE-EEE -E+EEE
% EOOEEEE (62-82):  EzzEEEE E+-EEEE E-+EEEE
% OEOEEEE (62-85):  zEzEEEE +E-EEEE -E+EEEE
% OOEEEEE (62-88):  zzEEEEE +-EEEEE -+EEEEE

%--------------------------------------------------------------------------
% 8-Clusters: 
%--------------------------------------------------------------------------
% EEEEEEEE (1)    :  EEEEEEE
% EEEEEEEO (2-4)  :  EEEEEEz EEEEEE+ EEEEEE-
% EEEEEEOE (5-7)  :  EEEEEzE EEEEE+E EEEEE-E
% EEEEEOEE (8-10) :  EEEEzEE EEEE+EE EEEE-EE
% EEEEOEEE (11-13):  EEEzEEE EEE+EEE EEE-EEE
% EEEOEEEE (14-16):  EEzEEEE EE+EEEE EE-EEEE
% EEOEEEEE (17-19):  EzEEEEE E+EEEEE E-EEEEE
% EOEEEEEE (17-19):  zEEEEEE +EEEEEE -EEEEEE
% OEEEEEEE (17-19):  zEEEEEE +EEEEEE -EEEEEE
% EEEEEEOO (20-22):  EEEEEzz EEEEE+- EEEEE-+
% EEEEEOEO (23-25):  EEEEzEz EEEE+E- EEEE-E+
% EEEEOEEO (26-28):  EEEzEEz EEE+EE- EEE-EE+
% EEEOEEEO (29-31):  EEzEEEz EE+EEE- EE-EEE+ 
% EEOEEEEO (32-34):  EzEEEEz E+EEEE- E-EEEE+ 
% EOEEEEEO (32-34):  zEEEEEz +EEEEE- -EEEEE+ 
% OEEEEEEO (32-34):  zEEEEEz +EEEEE- -EEEEE+ 
% EEEEEOOE (35-37):  EEEEzzE EEEE+-E EEEE-+E
% EEEEOEOE (38-40):  EEEzEzE EEE+E-E EEE-E+E
% EEEOEEOE (41-43):  EEzEEzE EE+EE-E EE-EE+E
% EEOEEEOE (44-46):  EzEEEzE E+EEE-E E-EEE+E
% EOEEEEOE (44-46):  zEEEEzE +EEEE-E -EEEE+E
% EOEEEEOE (44-46):  zEEEEzE +EEEE-E -EEEE+E
% OEEEEEOE (44-46):  zEEEEzE +EEEE-E -EEEE+E
% EEEEOOEE (47-49):  EEEzzEE EEE+-EE EEE-+EE
% EEEOEOEE (50-51):  EEzEzEE EE+E-EE EE-E+EE
% EEOEEOEE (53-55):  EzEEzEE E+EE-EE E-EE+EE
% EOEEEOEE (53-55):  zEEEzEE +EEE-EE -EEE+EE
% OEEEEOEE (53-55):  zEEEzEE +EEE-EE -EEE+EE
% EEEOOEEE (56-58):  EEzzEEE EE+-EEE EE-+EEE
% EEOEOEEE (59-61):  EzEzEEE E+E-EEE E-E+EEE
% EOEEOEEE (59-61):  zEEzEEE +EE-EEE -E+EEE
% OEEEOEEE (59-61):  zEEzEEE +EE-EEE -E+EEE
% EEOOEEEE (62-64):  EzzEEEE E+-EEEE E-+EEEE
% EOEOEEEE (62-64):  zEzEEEE +E-EEEE -E+EEEE
% OEEOEEEE (62-64):  zEzEEEE +E-EEEE -E+EEEE
% EOOEEEEE (62-64):  zzEEEEE +-EEEEE -+EEEEE
% OEOEEEEE (62-64):  zzEEEEE +-EEEEE -+EEEEE
% OOEEEEEE (62-112):  zzEEEEE +-EEEEE -+EEEEE

output = {SpinOp_1,SpinOp_2,SpinOp_3,SpinOp_4,SpinOp_5,SpinOp_6};

end




