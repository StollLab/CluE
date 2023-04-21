function [total_signal,auxiliary_signals,order_n_signals,batch_name] ...
  = calculate_signal_cluster_groups(System,Method,Nuclei,clusters,OutputData)


cluster_groups = form_cluster_groups(System,Method,Nuclei,clusters);
total_signal = 1;


methylTunnelingSplitting = Nuclei.methylTunnelingSplitting;
for igroup = 1:numel(cluster_groups)

  group_of_clusters = cluster_groups{igroup};

  if isfile(group_name)
    [group_signal,auxiliary_signals,order_n_signals,batch_name] ...
      = calculate_signal_ckpt(System,Method,Nuclei,group_of_clusters,...
      OutputData);

    movefile(batch_name,group_name)
  else
    group_signal = readmatrix(group_name);
  end


  if strcmp(Method.cluster_grouping,'methyl') && Method.Methyl.turn_groups_off
    Nuclei.methylTunnelingSplitting = 0*Nuclei.methylTunnelingSplitting;

    no_tunneling_group_name = ['no_tunneling_',group_name];

    [group_signal,~,~,batch_name] ...
      = calculate_signal_ckpt(System,Method,Nuclei,group_of_clusters,...
      OutputData);

    movefile(batch_name,no_tunneling_group_name);

    Nuclei.methylTunnelingSplitting = methylTunnelingSplitting;
  end

  total_signal = total_signal.*group_signal;

end
end

%-------------------------------------------------------------------------------
function cluster_groups = form_cluster_groups(System,Method,Nuclei,clusters)

switch Method.cluster_grouping
  case 'methyl'
    cluster_groups = form_methyl_cluster_groups(...
      System,Method,Nuclei,clusters);
    %case 'particle'
  otherwise
    error('Unrecognized cluster grouping "%s".',Method.cluster_grouping);
end
end
%-------------------------------------------------------------------------------
function cluster_groups = form_methyl_cluster_groups(...
  System,Method,Nuclei,clusters)



if ~Method.useMethylPseudoParticles || Method.order >= 6
  error('Grouping by methyl group assumes that each cluster contain only one methyl.');
end

max_size = numel(clusters);
group_ids = cell(max_size,1);
n_groups = 0;

% Identify cluster groups.
for cluster_size = 1:max_size
  n_clusters = size(clusters,1);
  group_ids{cluster_size} = zeros(n_clusters,1);
  group_dictionary = dictionary(string([]),[]);
  for icluster = 1:n_clusters

    cluster = clusters{cluster_size}(icluster,:);
    methyl_ids = unique(Nuclei.MethylID(cluster));
    key = sprintf('%d,',methyl_ids);

    if group_dictionary.isKey(key)
      id = group_dictionary(key);
    else
      n_groups  = n_groups + 1;
      id = n_groups;
      group_dictionary(key) = id;
    end

    group_ids{cluster_size}(icluster) = id;
  end
  assert(~any(group_ids{cluster_size}==0));
end


% Group clusters.
cluster_groups  = cell(n_groups,1);
for igroup = 1:n_groups
  cluster_groups{igroup} = cell(max_size,1);
  for cluster_size = 1:max_size
    sele = group_ids{cluster_size} == igroup;
    cluster_groups{igroup}{cluster_size} = clusters{cluster_size}(sele,:);
  end
end

%TODO save clusters

end