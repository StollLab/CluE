function [total_signal,auxiliary_signals,order_n_signals,batch_name] ...
  = calculate_signal_cluster_groups(System,Method,Nuclei,clusters,OutputData)


[total_signal,auxiliary_signals,order_n_signals,batch_name] ...
    = calculate_signal_ckpt(System,Method,Nuclei,clusters,OutputData);

assert(Method.ckptAuxiliarySignals);

methyl_file_name = [OutputData,'_methyls.csv'];
analyze_methyls(Nuclei,methyl_file_name)

[groups_ids,n_groups] = form_methyl_cluster_groups(Method,Nuclei,clusters);

partition_signals = ones(numel(total_signal),n_groups);
for cluster_size = 1:Method.order
  for ii = 1:numel(groups_ids{cluster_size})
    assert(abs(auxiliary_signals{cluster_size}(1,ii)-1)<1e-12);
    id = groups_ids{cluster_size}(ii);
    partition_signals(:,id) = partition_signals(:,id)...
      .*auxiliary_signals{cluster_size}(:,ii);
  end
end

total_signal_grps = prod(partition_signals,2).';
assert( max(max(abs(total_signal_grps-total_signal))) <1e-12 );
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
function [group_ids,n_groups] ...
    = form_methyl_cluster_groups(Method,Nuclei,clusters)

if ~Method.useMethylPseudoParticles || Method.order >= 6
  error(['Grouping by methyl group assumes that each cluster contains ',...
      'at most only one methyl.']);
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

end