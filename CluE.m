function [SignalMean, experiment_time, out] = CluE(Options,varargin)

if nargin==3
  System = Options;
  Options = struct('System',System,'Method',varargin{1},'Data',varargin{2});
  Options.mode = 'core';

elseif nargin==4
  System = Options;
  Options = struct('System',System,'Method',varargin{1},'Data',varargin{2},...
    'options',varargin{3});
  Options.mode = 'converge';

elseif nargin==8
  System = Options;
  Options = struct('System',System,'Method',varargin{1},'Data',varargin{2},...
    'savefile',varargin{3},...
    'N',varargin{4},...
    'dN',varargin{5},...
    'threshold',varargin{6},...
    'options',varargin{7});
  Options.mode = 'isotope Monte Carlo';
  
end

if ~isfield(Options,'mode')
  Options.mode = 'core';
end

if strcmp(Options.mode,'core')

  [SignalMean, experiment_time, ...
    out.TM_powder,out.Order_n_SignalMean,out.Nuclei,out.statistics] ...
    = CluE_core(Options.System,Options.Method,Options.Data);

elseif strcmp( Options.mode , 'isotope Monte Carlo')
[SignalMean,out.TM,out.statistics,experiment_time] = isotopeMonteCarlo(...
    Options.System,Options.Method,Options.Data,...
    Options.savefile,Options.N,Options.dN,Options.threshold,Options.options);

elseif strcmp( Options.mode , 'converge')
[out.parameters,out.System,out.Method,out.Data]  = converge_parameters(...
    Options.System,...
    Options.Method,...
    Options.Data, ...
    Options.options);

else
  error('The specified Options.mode, "%s", is unrecognized.\n',Options.mode);
end

end