function B = zeroDiag(A)
  N = min(size(A));
  B = A;
  for ii =1:N
    B(ii,ii) = 0;
  end

end