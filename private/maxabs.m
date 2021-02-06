function a = maxabs(A,varargin)

A_ = abs(A);
a = max(A_,varargin{:});
b = max(A,varargin{:});
N = numel(a);
for ii = 1:N

  if a(ii) ~= b(ii)
    a(ii) = -a(ii);
  end
  
end

end