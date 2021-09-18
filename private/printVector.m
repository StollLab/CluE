function printVector(vec, name)

N = numel(vec);
r = mod(N,10);
n = (N-r)/10;
fprintf([name,' = {\n']);
for ii = 1:n
  fprintf('%d, %d, %d, %d, %d, %d, %d, %d, %d, %d',vec( 1+ 10*(ii-1):10*ii ) )
  if ii==n && r == 0
    fprintf('\n');
  else
    fprintf(',\n');
  end
end

if r>0
  remStr = [ repmat('%d, ',1,r-1), '%d'];
  fprintf(remStr,vec( 10*n+1:end ))
end
fprintf('}.\n\n\n')
end