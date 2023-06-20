function save_tensors(Nuclei, System, Method,Data)

% Extract physical constants.
ge = System.ge;
geff = System.g(end);
muB = System.muB;
muN = System.muN;
mu0 = System.mu0;
hbar = System.hbar;

magneticField = System.magneticField;

Theory = System.Theory;
theory = Theory(Method.order,:);

thisCluster = Nuclei.Index;

[tensors,~] = pairwisetensors(Nuclei.Nuclear_g, ...
  Nuclei.Coordinates,thisCluster,Nuclei.Atensor,magneticField,ge,geff,...
  muB,muN,mu0,hbar,theory,0,0,0);

tensor_indices = [0, thisCluster];
for ii = tensor_indices
  idx = 0;
  if ii >= 1, idx = Nuclei.pdbID(ii); end
  fprintf('T[%i] = [%f, %f, %f];\n', idx,...
    tensors(1,1,ii+1,ii+1), tensors(2,2,ii+1,ii+1),tensors(3,3,ii+1,ii+1));
end

 
for ii = tensor_indices

  idx = 0;
  if ii >= 1, idx = Nuclei.pdbID(ii); end

  for jj = tensor_indices
    if jj < ii, continue; end

    jdx = 0;
    if jj >= 1, jdx = Nuclei.pdbID(jj); end

    fprintf('T[%i,%i] = [%f, %f, %f, %f, %f, %f];\n',...
      idx,jdx,...
      tensors(1,1,ii+1,jj+1), tensors(1,2,ii+1,jj+1),tensors(1,3,ii+1,jj+1),...
      tensors(2,2,ii+1,jj+1), tensors(2,3,ii+1,jj+1),tensors(3,3,ii+1,jj+1)...
      );
  end
end

if ii > 1
  ten = Nuclei.Qtensor(:,:,ii);
  if max(abs(ten)) > 1e-12
    fprintf('T[%i,%i] = [%f, %f, %f, %f, %f, %f];\n',...
      Nuclei.pdbID(ii),Nuclei.pdbID(ii),...
      ten(1,1), ten(1,2),ten(1,3),...
      ten(2,2), ten(2,3),ten(3,3)...
      );
  end
end
end