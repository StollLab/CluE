function Nuclei = setQuadrupoleTensor(e2qQh_,eta_,zQ,xQ,iNuc,Nuclei)
if norm(zQ)==0
  disp('Failed to set quadrupole tensor orientation.')
  return;
end
if norm(xQ)==0
  while xQ*zQ'==0
    disp('Failed to set quadrupole tensor orientation; using a random direction.')
    xQ = rand(1,3);
    xQ = xQ/norm(xQ);
  end
end

zQ = zQ/norm(zQ);
xQ = xQ - (xQ*zQ')*zQ;
xQ = xQ/norm(xQ);
yQ = cross(zQ,xQ);
yQ = yQ/norm(yQ);
R_Q2L = [xQ; yQ; zQ]; % rotation matrix from Q to lab frame

I = Nuclei.Spin(iNuc);

Qtensor_Q = e2qQh_/4/I/(2*I-1)*diag([-1+eta_, -1-eta_, 2]);
Qtensor_L = R_Q2L*Qtensor_Q*R_Q2L';

Nuclei.quadrupole2lab(:,:,iNuc) = R_Q2L;
Nuclei.Qtensor(:,:,iNuc) = Qtensor_L;
end