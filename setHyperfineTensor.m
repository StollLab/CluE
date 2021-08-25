function Atensor_L = setHyperfineTensor(Azz,Afc,zA,xA,iNuc,Nuclei)
if norm(zA)==0
  warning('Failed to set quadrupole tensor orientation.')
  return;
end
if norm(xA)==0
  while xA*zA==0
    warning('Failed to set quadrupole tensor orientation; using a random direction.')
    xA = rand(1,3);
    xA = xA/norm(xA);
  end
end

zA = zA/norm(zA);
xA = xA - (xA*zA')*zA;
xA = xA/norm(xA);
yA = cross(zA,xA);
yA = yA/norm(yA);
R_A2L = [xA; yA; zA]; % rotation matrix from Q to lab frame

I = Nuclei.Spin(iNuc);

Atensor_A = eye(3)*Afc + diag(Azz/2*[-1,-1,2]);
Atensor_L = R_A2L*Atensor_A*R_A2L';

% Nuclei.hyperfine2lab(:,:,iNuc) = R_A2L;
% Nuclei.Atensor(:,iNuc) = Atensor_L(:);
end