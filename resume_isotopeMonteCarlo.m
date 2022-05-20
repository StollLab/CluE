function [signals,TM,statistics,twotau] = resume_isotopeMonteCarlo(...
    savefile, ...
    newThreshold,...
    newOptions...
    )

% Load data.
load(savefile);

% Get list of fields from new options.
fields = fieldnames(newOptions);

% Overwrite target options.
for iopt = 1:numel(fields)
  options.(fields{iopt}) = newOptions.(fields{iopt});
end

% Run IMC.
[signals,TM,statistics,twotau] = isotopeMonteCarlo(...                  
    System,Method,Data,...                                                       
    savefile,N,dN,newThreshold,options);
end
