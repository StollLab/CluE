function logPrint(fileID,varargin)
try
    fprintf(fileID, varargin{:});
    fprintf(varargin{:});
catch 
  disp('logPrint() failed.');
end
end