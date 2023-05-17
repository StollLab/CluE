function [total_signal,auxiliary_signals,order_n_signals,batch_name] ...
  = calculate_signal_cluster_groups(System,Method,Nuclei,clusters,OutputData)


[total_signal,auxiliary_signals,order_n_signals,batch_name] ...
    = calculate_signal_ckpt(System,Method,Nuclei,clusters,OutputData);

assert(Method.ckptAuxiliarySignals);

methyl_file_name = [OutputData,'_methyls.csv'];
analyze_methyls(Nuclei,methyl_file_name)

[groups_ids,n_groups,group_names] ...
  = form_methyl_cluster_groups(Method,Nuclei,clusters);

partition_signals = ones(numel(total_signal),n_groups);
for cluster_size = 1:Method.order
  for ii = 1:numel(groups_ids{cluster_size})
    assert(abs(auxiliary_signals{cluster_size}(1,ii)-1)<1e-12);
    id = groups_ids{cluster_size}(ii);
    partition_signals(:,id) = partition_signals(:,id)...
      .*auxiliary_signals{cluster_size}(:,ii);
  end
end

save_name = [OutputData,'_methyl_partitions.csv'];
save_partition_signals(save_name,partition_signals,group_names);

total_signal_grps = prod(partition_signals,2).';
assert( max(max(abs(total_signal_grps-total_signal))) <1e-12 );
end
%-------------------------------------------------------------------------------
function save_partition_signals(save_name,partition_signals,group_names)
  T = array2table(partition_signals);
  T.Properties.VariableNames = group_names;
  writetable(T,save_name);
end

%-------------------------------------------------------------------------------
function [group_ids,n_groups,group_names] ...
    = form_methyl_cluster_groups(Method,Nuclei,clusters)

if ~Method.useMethylPseudoParticles || Method.order >= 6
  error(['Grouping by methyl group assumes that each cluster contains ',...
      'at most only one methyl.']);
end

max_size = numel(clusters);
group_ids = cell(max_size,1);
group_names = cell(max_size,1); 
n_groups = 0;

group_dictionary = dictionary(string([]),[]);

% Identify cluster groups.
for cluster_size = 1:max_size
  n_clusters = size(clusters{cluster_size},1);
  group_ids{cluster_size} = zeros(n_clusters,1);
  for icluster = 1:n_clusters

    cluster = clusters{cluster_size}(icluster,:);
    methyl_ids = unique(Nuclei.MethylID(cluster));
    methyl_ids = methyl_ids(methyl_ids > 0);
    key = sprintf('%d,',methyl_ids);

    if group_dictionary.isKey(key)
      id = group_dictionary(key);
    else
      n_groups  = n_groups + 1;
      id = n_groups;
      group_names{id} = get_partition_name(methyl_ids,Nuclei);
      group_dictionary(key) = id;
    end

    group_ids{cluster_size}(icluster) = id;
  end
  assert(~any(group_ids{cluster_size}==0));
end

end
%-------------------------------------------------------------------------------
function name = get_partition_name(methyl_ids,Nuclei)

  if isempty(methyl_ids) || all(methyl_ids == 0)
    name = 'bath';
    return;
  end

  
  name = '';
  for id = methyl_ids
     if id==0, continue; end

     hydrogens = find(Nuclei.MethylID == id);
     assert(numel(hydrogens)==3);
     methyl_str = sprintf('methyl_%i_%i_%i',...
       hydrogens(1),hydrogens(2),hydrogens(3));
     if isempty(name)
       name = methyl_str;
     else
       name = [name,'_',methyl_str];
     end
  end
end