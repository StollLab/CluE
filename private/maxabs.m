function a = maxabs(A,varargin)

A_ = abs(A);
a = max(A_,varargin{:});
b = max(A,varargin{:});
n = a~=b;
a(n) = -a(n);


end