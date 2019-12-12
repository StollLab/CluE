function NCK = NchooseK(N,K)
nN = numel(N);
nK = numel(K);
NCK = zeros(nN,nK);

for n=1:nN
  for k=1:nK
    
    if K(k) > N(n) || K(k) < 0
      NCK(n,k) = 0;
    else
      NCK(n,k) = nchoosek(N(n),K(k));
    end
    
  end
end

end
