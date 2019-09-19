function n = NchooseK(N,K)

if K > N || K < 0
  n = 0;
else 
  n = nchoosek(N,K);
end
end
