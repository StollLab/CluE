function a = minabs(A,varargin)
a = 1./maxabs(1./A,varargin{:});
end