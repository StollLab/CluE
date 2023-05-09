function [total_signal,auxiliary_signals,order_n_signals,batch_name] ...
  = calculate_signal_cluster_groups(System,Method,Nuclei,clusters,OutputData)

methyl_file_name = [OutputData,'_methyls.csv'];
analyze_methyls(Nuclei,methyl_file_name)

cluster_groups = form_methyl_cluster_groups(Method,Nuclei,clusters);
total_signal = 1;


for igroup = 1:numel(cluster_groups)

  group_of_clusters = cluster_groups{igroup};

  group_name = [OutputData,'_group_',int2str(igroup),'.csv'];

  if ~isfile(group_name) || ~Method.partialSave
    [group_signal,~,~,batch_name] ...
      = calculate_signal_ckpt(System,Method,Nuclei,group_of_clusters,...
      OutputData);

    save_group_signal(group_name,group_signal,igroup,group_of_clusters);
    if isfile(batch_name)
      delete(batch_name)
    end
  else
    group_signal = readmatrix(group_name);
  end


  total_signal = total_signal.*group_signal;

end
auxiliary_signals =[];
order_n_signals = [];
batch_name = [];
end
%-------------------------------------------------------------------------------
function save_group_signal(group_name,grp_sig,grp_idx,group_of_clusters)
  if ~isempty(group_of_clusters{1}) ...
      || numel(group_of_clusters)<3 || numel(group_of_clusters{3}) ~= 3
    grp_var = ['group_',int2str(grp_idx)];

  else
    cluster = group_of_clusters{3};
    grp_var = ['methyl_', sprintf('%d_',cluster)];
  end
  if size(grp_sig,1) == 1
    grp_sig = grp_sig.';
  end
  T = array2table(grp_sig);
  T.Properties.VariableNames(1) = {grp_var};
  writetable(T,group_name);
end

%-------------------------------------------------------------------------------
function cluster_groups = form_methyl_cluster_groups(Method,Nuclei,clusters)

if ~Method.useMethylPseudoParticles || Method.order >= 6
  error('Grouping by methyl group assumes that each cluster contain only one methyl.');
end

max_size = numel(clusters);
group_ids = cell(max_size,1);
n_groups = 0;

group_dictionary = dictionary(string([]),[]);

% Identify cluster groups.
for cluster_size = 1:max_size
  n_clusters = size(clusters{cluster_size},1);
  group_ids{cluster_size} = zeros(n_clusters,1);
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

for cluster_size = 1:max_size
  n_clusters = 0;
  for igroup = 1:n_groups
    n_clusters = n_clusters + size(cluster_groups{igroup}{cluster_size},1);
  end
  assert(n_clusters==size(clusters{cluster_size},1));
end
%TODO save clusters
% save_clusters(filename,clusters)

end