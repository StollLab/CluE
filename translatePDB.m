function pdb = translatePDB(pdb,System)

Coordinates = [pdb.x,pdb.y,pdb.z];

originVec = = getRotationOrigin(System,Coordinates);

Coordinates = Coordinates - originVec;

R = rotateZYZ(System.pdbAlpha,System.pdbBeta,System.pdbGamma);

Coordinates = R*Coordinates*R';

Coordinates = Coordinates + originVec;

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function Rotation_Coordinates = getRotationOrigin(System,Coordinates)
 
Rotation_Coordinates = System.pdbTranslation;

if iscell(Rotation_Coordinates)
  % find nuclei to average over
  replaceNuclei = [Rotation_Coordinates{:}];
  ReplaceNuclei = zeros(size(replaceNuclei));
  for irep = 1:length(replaceNuclei)
    ReplaceNuclei(irep) = find(pdb.serial==replaceNuclei(irep));
  end
  % place the electron at the mean coordinates
  System.Electron.Coordinates = mean( Coordinates(ReplaceNuclei,:),1);
  
  % set initial electron coordinates to a 3-vector
  Rotation_Coordinates = System.Electron.Coordinates;
  
  
  
elseif length(System.Electron.Coordinates) == 1
  
  replaceNucleus = System.Electron.Coordinates;
  System.Electron.Coordinates = Coordinates(replaceNucleus,:);
  Rotation_Coordinates = Coordinates(replaceNucleus,:);
  
elseif length(System.Electron.Coordinates) == 2
  
  replaceNucleus1 = System.Electron.Coordinates(1);
  replaceNucleus2 = System.Electron.Coordinates(2);
  System.Electron.Coordinates = ...
    0.5*( Coordinates(replaceNucleus1,:) +Coordinates(replaceNucleus2,:));
  Rotation_Coordinates = System.Electron.Coordinates;
  
end
if all( size(Rotation_Coordinates)==[3,1]) 
  return;
elseif all( size(Rotation_Coordinates)==[1,3])
  Rotation_Coordinates = Rotation_Coordinates';
else 
  error(['Error in getElectronCoordinates():', ...
    'could not get electron coordinates.']);
end

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end