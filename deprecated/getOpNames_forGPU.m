
%==========================================================================
% Generate a list of strings specifying the cluster spin-operator.
%==========================================================================

function opNames = getOpNames_forGPU(clusterSize)

% Get number of operators for the given cluster size.
numOps = 1 + 3*clusterSize + 9*NchooseK(clusterSize,2);

% Initialize output.
opNames = cell(1,numOps);

% Begin the identity string.
opNames{1} = 'E';

% Get the 1 identity operator EE....EE.
for iSpin = 2:clusterSize
  opNames{1} = [opNames{1},'E'];
end

% Set index.
idx = 1;

Op = {'Z','R','L'};
% Add 3*clusterSize operators of the form EE...O...EE,
% where O in {Z,R,L}.
% Loop over spins.
for iSpin = 1:clusterSize
  
  % Adjust for placement order.
  ispin_index = clusterSize-iSpin + 1;
  
  % Loop over the operator set indices.
  for iop = 1:3
    idx = idx + 1;
    % Initialize operator string with the identity string.
    opNames{idx} = opNames{1};
    % Set the non-identity element.
    opNames{idx}(ispin_index) = Op{iop};
  end
end

OpOp = {'ZZ','RL','LR'};

% Add 3*NchooseK(clusterSize,2) operators of the form EE...O1...E...O2...EE,
% where O1,2 in OpOp{}(1,2).
% Loop over spins.
for jSpin = 1:clusterSize-1
  
  % Adjust for placement order.
  jspin_index = clusterSize-jSpin + 1;
  
  % Loop over spins.
  for iSpin = jSpin+1:clusterSize
    
    % Adjust for placement order.
    ispin_index = clusterSize-iSpin + 1;
    
    % Loop over the operator set indices.
    for iop = 1:3
      
      % Increment index.
      idx = idx + 1;
      % Initialize operator string with the identity string.
      opNames{idx} = opNames{1};
      % Set the non-identity elements.
      opNames{idx}(ispin_index) = OpOp{iop}(1);
      opNames{idx}(jspin_index) = OpOp{iop}(2);
    end
  end
end

%

Op_CD = {'ZR', 'ZL', 'RZ', 'LZ'};

% Add 4*NchooseK(clusterSize,2) operators of the form EE...O1...E...O2...EE,
% where O1,2 in Op_CD{}(1,2).
% Loop over spins.
for jSpin = 1:clusterSize-1
  
  % Adjust for placement order.
  jspin_index = clusterSize-jSpin + 1;
  
  % Loop over spins.
  for iSpin = jSpin+1:clusterSize
    
    % Adjust for placement order.
    ispin_index = clusterSize-iSpin + 1;
    
    % Loop over the operator set indices.
    for iop = 1:4
      
      % Increment index.
      idx = idx + 1;
      % Initialize operator string with the identity string.
      opNames{idx} = opNames{1};
      % Set the non-identity elements.
      opNames{idx}(ispin_index) = Op_CD{iop}(1);
      opNames{idx}(jspin_index) = Op_CD{iop}(2);
    end
  end
end

% RR and LL
% Add 4*NchooseK(clusterSize,2) operators of the form EE...O...E...O...EE,
% where O1,2 in {R,L}.
% Loop over spins.
for jSpin = 1:clusterSize-1
  
  % Adjust for placement order.
  jspin_index = clusterSize-jSpin + 1;
  
  % Loop over spins.
  for iSpin = jSpin+1:clusterSize
    
    % Adjust for placement order.
    ispin_index = clusterSize-iSpin + 1;
   
    % Loop over the operator set indices.
    for iop = 2:3
      
      % Increment index.
      idx = idx + 1;
      % Initialize operator string with the identity string.
      opNames{idx} = opNames{1};
      % Set the non-identity elements.
      opNames{idx}(ispin_index) = Op{iop};
      opNames{idx}(jspin_index) = Op{iop};
    end
  end
end

end

% Information on Operator Indices

%--------------------------------------------------------------------------
% 1-Clusters: E O
%--------------------------------------------------------------------------
% E (1)  : E
% O (2-4): z + -

%--------------------------------------------------------------------------
% 2-Clusters: EE EO OE OO
%--------------------------------------------------------------------------

% EE (1)   : EE
% EO (2-4) : Ez E+ E-
% OE (5-7) : zE +E -E
% OO (8-10): zz +- -+

% OO (11-14): z+ z- +z -z
% OO (15-16): ++ -- 

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

% EEE =  1;
% EEZ =  2; EER =  3; EEL =  4;
% EZE =  5; ERE =  6; ELE =  7;
% ZEE =  8; REE =  9; LEE = 10;
% EZZ = 11; ERL = 12; ELR = 13;
% ZEZ = 14; REL = 15; LER = 16;
% ZZE = 17; RLE = 18; LRE = 19;
% 
% EZR = 20;  EZL = 21;  ERZ = 22;  ELZ = 23;
% ZER = 24;  ZEL = 25;  REZ = 26;  LEZ = 27;
% ZRE = 28;  ZLE = 29;  RZE = 30;  LZE = 31;
% 
% ERR = 32;  ELL = 33;
% RER = 34;  LEL = 35;
% RRE = 36;  LLE = 37;

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

% EEEE =  1;
% EEEZ =  2; EEER =  3; EEEL =  4;
% EEZE =  5; EERE =  6; EELE =  7;
% EZEE =  8; EREE =  9; ELEE = 10;
% ZEEE = 11; REEE = 12; LEEE = 13;
% EEZZ = 14; EERL = 15; EELR = 16;
% EZEZ = 17; EREL = 18; ELER = 19;
% ZEEZ = 20; REEL = 21; LEER = 22;
% EZZE = 23; ERLE = 24; ELRE = 25;
% ZEZE = 26; RELE = 27; LERE = 28;
% ZZEE = 29; RLEE = 30; LREE = 31;
% 
% 
% EEZR = 31;  EEZL = 32;  EERZ = 33;  EELZ = 34;
% EZER = 35;  EZEL = 36;  EREZ = 37;  ELEZ = 38;
% ZEER = 39;  ZEEL = 40;  REEZ = 41;  LEEZ = 42;
% EZRE = 43;  EZLE = 44;  ERZE = 45;  ELZE = 46;
% ZERE = 47;  ZELE = 48;  REZE = 49;  LEZE = 50;
% ZREE = 51;  ZLEE = 52;  RZEE = 53;  LZEE = 54;
% 
% EERR = 55;  EELL = 56;
% ERER = 57;  ELEL = 58;
% REER = 59;  LEEL = 60;
% ERRE = 61;  ELLE = 62;
% RERE = 63;  LELE = 64;
% RREE = 65;  LLEE = 66;


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