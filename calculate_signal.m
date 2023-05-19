% calculateSignal Calculates the echo signal
%
%  [Signal,AuxSignals,Signals] = calculateSignal(System,Method,Nuclei,Clusters)
%
% Input:
%   System   system structure
%   Method   method structure
%   Nuclei   nuclei structure
%   Clusters information about clusters
%
% Output:
%   Signal
%   AuxSignals
%   Signals

% Clusters = Clusters(cluster index , 1:size ,order)
% Clusters(cluster index , size > order ,order) = 0.

function [total_signal,auxiliary_signals,order_n_signals,batch_name] ... 
       = calculate_signal(System,Method,Nuclei,clusters,OutputData)

if Method.use_calculate_signal_ckpt
  if Method.use_calculate_signal_cluster_groups
    [total_signal,auxiliary_signals,order_n_signals,batch_name] ...
      = calculate_signal_cluster_groups(System,Method,Nuclei,clusters,...
      OutputData);
  else
    [total_signal,auxiliary_signals,order_n_signals,batch_name] ...
      = calculate_signal_ckpt(System,Method,Nuclei,clusters,OutputData);
  end
else
  [total_signal,auxiliary_signals,order_n_signals] ...
    = calculate_signal_default(System,Method,Nuclei,clusters);
  batch_name = [];
end

end




