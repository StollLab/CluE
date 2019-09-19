function v= stretchedExp(Sys,varargin) 
  v = exp(-abs(Sys.twotau/Sys.TM).^Sys.k);
end
