function logPrint(fileID,varargin)
    fprintf(fileID, varargin{:});
    fprintf(varargin{:});
end