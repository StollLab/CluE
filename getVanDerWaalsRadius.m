% Data from Mathematica 12.1 ElementData[].
function R = getVanDerWaalsRadius(type)
pm = 1e-12;
switch type
  case {'H','D'}
    R = 120*pm;
  case 'He'
    R = 140*pm;
  case 'Li'
    R = 182*pm;
  case 'C'
    R = 170*pm;
  case 'N'
    R = 155*pm;
  case 'O'
    R = 152*pm;
  case 'F'
    R = 147*pm;
  case 'Ne'
    R = 154*pm;
  case 'Na'
    R = 227*pm;
  case 'Mg'
    R = 173*pm;
  case 'Si'
    R = 210*pm;
  case 'P'
    R = 180*pm;
  case 'S'
    R = 180*pm;
  case 'Cl'
    R = 175*pm;
  case 'Ar'
    R = 188*pm;
  case 'K'
    R = 275*pm;
  case 'Ni'
    R = 163*pm;
  case 'Cu'
    R = 140*pm;
  case 'Zn'
    R = 139*pm;
  case 'Ga'
    R = 187*pm;
  case 'As'
    R = 185*pm;
  case 'Se'
    R = 190*pm;
  case 'Br'
    R = 185*pm;
  case 'Kr'
    R = 202*pm;
  case {'M','MW','LP'}
    R = 0;
  otherwise
    err= ['Could not find Van der Waals radius for ', type,'.'];
    error(err)
end
end