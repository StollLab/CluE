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


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function [total_signal,auxiliary_signals,order_n_signals,batch_name] ...
  = calculate_signal_ckpt(System,Method,Nuclei,clusters,OutputData)

% Unpack spin operators.
Op = Nuclei.SpinOperators;
SpinXiXjOps = Nuclei.SpinXiXjOperators;

if any(System.Theory(:,10))
  [~,Nuclei.ZeemanSpinStates] = setRandomZeemanState(Nuclei);
end

dimensionality = 1 + strcmp(System.experiment,'CPMG-2D');
num_time_points = sum(System.nPoints)^dimensionality;


cluster_map = build_cluster_map(clusters);

order_n_signals = ones(num_time_points,Method.order);

[auxiliary_signals,loaded] = load_auxiliary_signals(...
  Nuclei.numberClusters,num_time_points,Method,OutputData);

batch_name = [];
for clusterSize = 1:Method.order

  if loaded(clusterSize)
    order_n_signals(:,clusterSize) = prod(auxiliary_signals{clusterSize},2);
  else
    [order_n_signals(:,clusterSize),auxiliary_signals,batch_name] ...
      = process_clusters(auxiliary_signals,...
      cluster_map,clusterSize,num_time_points,...
      Op,SpinXiXjOps,clusters,Nuclei,System,Method,OutputData);

    if Method.ckptAuxiliarySignals
      save_auxiliary_signals(auxiliary_signals{clusterSize},...
        clusterSize,clusters,OutputData);
    end
  end

  if clusterSize > 1
    order_n_signals(:,clusterSize) ...
      = order_n_signals(:,clusterSize).*order_n_signals(:,clusterSize-1);
  end
end


order_n_signals = order_n_signals.';
total_signal = order_n_signals(Method.order,:);


end
%-------------------------------------------------------------------------------
function [auxiliary_signals,successfully_loaded] ...
  =  load_auxiliary_signals(...
  numberClusters,num_time_points,Method,OutputData)

successfully_loaded = false(Method.order,1);
auxiliary_signals = cell(Method.order,1);
for clusterSize = 1:Method.order

  aux_file = get_auxiliary_file(clusterSize,OutputData);

  if Method.ckptAuxiliarySignals && isfile(aux_file)
    auxiliary_signals{clusterSize} = readmatrix(aux_file);    
    successfully_loaded(clusterSize) = true;
  else
    numClusters = numberClusters(clusterSize);
    auxiliary_signals{clusterSize} = zeros(num_time_points,numClusters);
  end

end
end
%------------------------------------------------------------------------------- 
function save_auxiliary_signals(...
  auxiliary_signals,cluster_size,clusters,OutputData)
  aux_file = get_auxiliary_file(cluster_size,OutputData);
  % writematrix(auxiliary_signals,aux_file);

  
  T = array2table(auxiliary_signals);

  num_clusters = size(clusters{cluster_size},1);
  for icluster = 1: num_clusters
    cluster = clusters{cluster_size}(icluster,:);
    var_name = ['clu_',sprintf('%d_',cluster)];
    T.Properties.VariableNames(icluster) = {var_name};
  end
  writetable(T,aux_file);
end
%------------------------------------------------------------------------------- 
function aux_file = get_auxiliary_file(clusterSize,OutputData)
aux_file = ['aux_',OutputData,'_clustersize_',...
  int2str(clusterSize), '.csv'];
end
%------------------------------------------------------------------------------- 
function [cum_prod_aux_sig, auxiliary_signals,batch_name] = process_clusters(...
  auxiliary_signals,cluster_map,clusterSize,...
  num_time_points,Op,SpinXiXjOps,clusters,Nuclei,System,Method,...
  OutputData)


batch_size = Method.batch_size;

numClusters = size(clusters{clusterSize},1);
auxiliary_signals{clusterSize} = zeros(num_time_points,numClusters);

is_highest_order = clusterSize == Method.order;
batch_name = [];
if is_highest_order
  batch_name = ['temp_',OutputData,'_batch_',int2str(clusterSize),...
    '_',int2str(batch_size), '.csv'];

  [cum_prod_aux_sig,end_cluster] = load_batch_ckpt(batch_name,num_time_points);
else
  cum_prod_aux_sig = ones(num_time_points,1);
  end_cluster = 0;
end

while end_cluster < numClusters
  start_cluster = end_cluster + 1;
  end_cluster = min(end_cluster + batch_size,numClusters);

  auxiliary_signals{clusterSize}(:,start_cluster:end_cluster) ...
    = calculate_cluster_signals(start_cluster,end_cluster,...
    clusterSize,num_time_points,Op,SpinXiXjOps,clusters,Nuclei,System,Method);

  auxiliary_signals = update_auxiliary_signals(auxiliary_signals,...
    start_cluster,end_cluster,...
    clusterSize, clusterSize,cluster_map,clusters);

  cum_prod_aux_sig = cum_prod_aux_sig...
    .*prod(auxiliary_signals{clusterSize}(:,start_cluster:end_cluster),2);

  if is_highest_order
    save_batch_ckpt(batch_name,cum_prod_aux_sig,end_cluster);

    if Method.save_partial_cce_product
      partial_name = ['partial_',OutputData,'_nClu_',int2str(end_cluster),...
        '.csv'];
      T = array2table(cum_prod_aux_sig);
      part_var = ['partial_sig_',int2str(end_cluster)];
      T.Properties.VariableNames(1) = {part_var};
      writetable(T,partial_name);

    end
  end
end


end
%-------------------------------------------------------------------------------
function [cum_prod_aux_sig,end_cluster] = load_batch_ckpt(batch_name,...
  num_time_points)

if ~isfile(batch_name)
  cum_prod_aux_sig = ones(num_time_points,1);
  end_cluster = 0;
  return;
end

A = readtable(batch_name);
end_cluster = str2double(A.Properties.VariableNames{1}(7:end));
cum_prod_aux_sig = table2array(A);

end
%-------------------------------------------------------------------------------
function save_batch_ckpt(batch_name,cum_prod_aux_sig,end_cluster)
  T = array2table(cum_prod_aux_sig);
  batch_var = ['batch_',int2str(end_cluster)];
  T.Properties.VariableNames(1) = {batch_var};
  writetable(T,batch_name);
end
%-------------------------------------------------------------------------------
function signal = calculate_cluster_signal(thisCluster,num_time_points,...
  Op,SpinXiXjOps,Nuclei,System,Method)

thisClusterSize = numel(thisCluster);
% Select the appropriate spin operator.
spinMultiplicity = Nuclei.StateMultiplicity(thisCluster(1));
SpinOp = Op{spinMultiplicity}{thisClusterSize};
if spinMultiplicity==3
  SpinXiXjOp = SpinXiXjOps{thisClusterSize};
else
  SpinXiXjOp = [];
end

signal = zeros(num_time_points,1);

for iave = 1:System.nStates(numel(thisClusterSize))

  [Ha,Hb,Ha_TR,Hb_TR] = get_Hamiltonians(thisCluster,SpinOp,SpinXiXjOp,...
    Nuclei,System,Method);

  % Calculate cluster coherence.
  new_sig = propagate(Hb,Ha,Hb_TR,Ha_TR, ...
    thisCluster, System);

  signal = signal  +  new_sig;

end
end
%-------------------------------------------------------------------------------
function cluster_signals = calculate_cluster_signals(...
  start_cluster,end_cluster,clusterSize,num_time_points,...
  Op,SpinXiXjOps,clusters,Nuclei,System,Method)

cluster_signals = zeros(num_time_points,end_cluster + 1 - start_cluster);

for iCluster = start_cluster:end_cluster
  idx = iCluster + 1 - start_cluster;

  thisCluster = extract_cluster(iCluster,clusterSize,...
    clusters,Nuclei,System);

  if do_skip_cluster(thisCluster,Nuclei,System)
    cluster_signals(:,idx) = 1.0;
    continue;
  end

  cluster_signals(:,idx) = calculate_cluster_signal(thisCluster,...
    num_time_points,Op,SpinXiXjOps,Nuclei,System,Method);

end



end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Generate Density Matrix
function densityMatrix = getDensityMatrix(states,multiplicities,Cluster)

clustersize = length(Cluster);

prod_state = states(clustersize);

offset_factor = multiplicities(clustersize);
for iSpin = clustersize-1:-1:1
  prod_state = prod_state + offset_factor*(states(iSpin) - 1);
  offset_factor = offset_factor*multiplicities(iSpin);
end
densityMatrix = zeros(offset_factor);
densityMatrix(prod_state,prod_state) = 1;

end
%-------------------------------------------------------------------------------
function  [Ha,Hb,Ha_TR,Hb_TR] = get_Hamiltonians(...
  thisCluster,SpinOp,SpinXiXjOp,Nuclei,System,Method)

B1y =  System.RF.B1y;
B1x =  System.RF.B1x;
nuRF = System.RF.nuRF;

% Methyl Groups
isMethylHydron = strcmp(Nuclei.Type,'CH3_1H');

theory = System.Theory(Method.order,:);


% Get interaction tensors.
[tensors,~] = pairwisetensors(Nuclei.Nuclear_g, ...
  Nuclei.Coordinates,thisCluster,Nuclei.Atensor,...
  System.magneticField,System.ge,System.gMatrix(3,3),...
  System.muB,System.muN,System.mu0,System.hbar,theory,...
  B1x,B1y,nuRF);


qtensors = Nuclei.Qtensor(:,:,thisCluster);


% Build Hamiltonians.
[Ha,Hb] = ...
  assembleHamiltonian(Nuclei.StateMultiplicity(thisCluster),...
  tensors,SpinOp,qtensors,SpinXiXjOp, theory,...
  isMethylHydron(thisCluster),...
  System.Methyl.method,...
  Nuclei.methylTunnelingSplitting(thisCluster,thisCluster),...
  Nuclei.MethylID(thisCluster) ...
  );


if strcmp(System.experiment,'Hahn-TR')

  B1x2 =  System.RF.B1x(2);
  B1y2 =  System.RF.B1y(2);
  nuRF2 = System.RF.nuRF(2);

  [tensors_TR,~] = pairwisetensors(Nuclei.Nuclear_g, ...
    Nuclei.Coordinates,thisCluster,Nuclei.Atensors,...
    System.magneticField,System.ge,System.gMatrix(3,3),...
    System.muB,System.muN,System.mu0,System.hbar,theory,...
    B1x2,B1y2,nuRF2);

  [Ha_TR,Hb_TR] = ...
    assembleHamiltonian(Nuclei.StateMultiplicity(thisCluster),...
    tensors_TR,SpinOp,qtensors,SpinXiXjOp, theory,...
    isMethylHydron(thisCluster), ...
    System.Methyl.method, ...
    Nuclei.methylTunnelingSplitting(thisCluster,thisCluster),...
    Nuclei.MethylID(thisCluster));
else
  Ha_TR = [];
  Hb_TR = [];
end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Propagate Function
function Signal = propagate(Ham_beta,Ham_alpha,...
  Ham_beta_TR,Ham_alpha_TR ,...
  thisCluster,System)

% Extract experimental parameters.
dt = System.dt(1);
dt2 = System.dt(2);
N1 = System.nPoints(1);
% N2 = System.nPoints(2);
Npoints = sum(System.nPoints);
nPulses = System.nPulses;


% Ensure sub-Hamiltonians are Hermitian.
Ham_beta = (Ham_beta+Ham_beta')/2;
Ham_alpha = (Ham_alpha+Ham_alpha')/2;


% Get density matrix.
if System.useThermalEnsemble
  betaT = 2*pi*System.hbar/System.kT; % 1/Hz.
  DensityMatrix = propagator_eig((Ham_alpha+Ham_beta)/2,-1i*betaT,0);
else
  DensityMatrix = getDensityMatrix(ZeemanStates(thisCluster),...
    Nuclei.StateMultiplicity(thisCluster),thisCluster);
end

% Vectorize density matrix.
vecDensityMatrixT = reshape(DensityMatrix.',1,[])/trace(DensityMatrix);

t_Uhrig = 0;

% Generate propagators for time steps dt and dt2
[dU_beta,dU_beta2] = propagator_eig(Ham_beta,dt,dt2);
[dU_alpha,dU_alpha2] = propagator_eig(Ham_alpha,dt,dt2);

doTR = ~isempty(Ham_beta_TR);
if doTR
  % Ensure subHamiltonians are Hermitian
  Ham_beta_TR = (Ham_beta_TR+Ham_beta_TR')/2;
  Ham_alpha_TR = (Ham_alpha_TR+Ham_alpha_TR')/2;
  % Generate propagators for time steps dt and dt2
  [dU_beta_TR, dU_beta2_TR] = propagator_eig(Ham_beta_TR,dt,dt2);
  [dU_alpha_TR, dU_alpha2_TR] = propagator_eig(Ham_alpha_TR,dt, dt2);
end

% Initialize propagators
nStates = length(Ham_beta);
U_beta = eye(nStates);
U_alpha = eye(nStates);

% Initialize signal
if strcmp(System.experiment,'CPMG-2D')
  v = ones(Npoints,Npoints);
else
  v = ones(Npoints,1);
end
if strcmp(System.experiment,'CPMG-const')
  total_time = 2*System.nPulses*System.Time(end);
end
% Loop over time points
for iTime = 1:Npoints
  switch System.experiment

    case 'FID'
      U = U_beta'*U_alpha;
      v(iTime) = vecDensityMatrixT*U(:);

    case 'Hahn'
      %       U = U_beta'*U_alpha'*U_beta*U_alpha;
      U = (U_alpha*U_beta)'*U_beta*U_alpha;
      v(iTime) = vecDensityMatrixT*U(:);

    case 'CPMG'

      U = U_alpha'*U_beta'  * ...
        (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha*U_beta;
      v(iTime) = vecDensityMatrixT*U(:);

    case 'CPMG-const'

      % THIS NEEDS TO BE UPDATED TO USE dt2 WHEN iTime > N1.
      [U_beta, U_beta_2] = propagator_eig(...
        Ham_beta,(iTime-1)*dt, total_time/4-(iTime-1)*dt);
      [U_alpha, U_alpha_2] = propagator_eig(...
        Ham_alpha,(iTime-1)*dt, total_time/4-(iTime-1)*dt);

      U = U_alpha_2'*U_beta_2'  *  ...
        (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;

      v(iTime) = vecDensityMatrixT*U(:);

    case 'CPMG-2D'

      U_beta_2 = eye(nStates);
      U_alpha_2 = eye(nStates);

      % Loop over second experimental dimension.
      for jTime = 1:iTime

        % Generate time dependent detection operator.
        U = U_alpha_2'*U_beta_2'  *  ...
          (U_beta'*U_alpha' * U_beta*U_alpha)  * U_alpha_2*U_beta_2;

        v(iTime,jTime) = vecDensityMatrixT*U(:);
        v(jTime,iTime) = v(iTime,jTime)';

        % Increment propagator.
        if jTime < N1
          U_beta_2 = dU_beta*U_beta_2;
          U_alpha_2 = dU_alpha*U_alpha_2;
        else
          U_beta_2 = dU_beta2*U_beta_2;
          U_alpha_2 = dU_alpha2*U_alpha_2;
        end

      end

    case 'Hahn-TR'
      U = U_beta'*U_beta_TR'*U_alpha'*U_alpha_TR' ...
        * U_beta_TR*U_beta*U_alpha_TR*U_alpha;

      v(iTime) = vecDensityMatrixT*U(:);

    case 'CP_N'
      U_aa = U_alpha*U_alpha;
      U_bb = U_beta*U_beta;

      exponent = fix((nPulses-1)/2);
      AABB = (U_aa*U_bb)^exponent;
      BBAA = (U_bb*U_aa)^exponent;
      if mod(nPulses,2)==0
        U = ( U_alpha*BBAA*U_bb*U_alpha )' ...
          *( U_beta*AABB*U_aa*U_beta  );
      else
        U = (  U_beta*AABB*U_alpha )' ...
          *( U_alpha*BBAA*U_beta  );

      end

      v(iTime) = vecDensityMatrixT*U(:);

    case 'Uhrig_N'
      % t_i = (2*tau)*sin^2[ (pi*i)/(2*(N+1)) ].

      if iTime< N1
        t_Uhrig = t_Uhrig + dt*2*nPulses;
      else
        t_Uhrig = t_Uhrig + dt2*2*nPulses;
      end

      v(iTime) = getUhrigN(t_Uhrig, Ham_beta,Ham_alpha,...
        vecDensityMatrixT, nPulses);

  end

  % Increment propagators
  if iTime < N1
    U_beta = dU_beta*U_beta;
    U_alpha = dU_alpha*U_alpha;
  else
    U_beta = dU_beta2*U_beta;
    U_alpha = dU_alpha2*U_alpha;
  end


  if doTR
    if iTime < N1
      U_beta_TR = dU_beta_TR*U_beta_TR;
      U_alpha_TR = dU_alpha_TR*U_alpha_TR;
    else
      U_beta_TR = dU_beta2_TR*U_beta_TR;
      U_alpha_TR = dU_alpha2_TR*U_alpha_TR;
    end
  end

end

% Return signal as a vector
if strcmp(System.experiment,'CPMG-2D')
  Signal = reshape(v.',1,[]);
else
  Signal = v;
end

if any(abs(v) - 1 > 1e-9)
  error('Coherence error: coherence cannot be larger than 1.');
end

end
%-------------------------------------------------------------------------------
% Calculate propagator using diagonalization
function [U1,U2] = propagator_eig(Ham,dt1,dt2)
[EigenVectors,EigenValues] = eig(Ham);
Udiag1 = exp(-2i*pi*diag(EigenValues)*dt1);
U1 = EigenVectors*diag(Udiag1)*EigenVectors';

if nargin>2 && dt2>0
  Udiag2 = exp(-2i*pi*diag(EigenValues)*dt2);
  U2 = EigenVectors*diag(Udiag2)*EigenVectors';
else
  U2 = 1;
end
end
%-------------------------------------------------------------------------------
function v_iTime = getUhrigN(total_time, Hamiltonian_beta,Hamiltonian_alpha,...
  vecDensityMatrixT, nPulses)
% t_i = T*sin^2[ (pi*i)/(2*(N+1)) ].

% delays = T*1/2*(cos(pi*([indices,N+1]-1)/( (N+1) ) ) ...
%        - cos(pi*[indices,N+1]/( (N+1) ) ));

% Since the Uhrig pulse sequence is symmetric about the  half-way point,
% only half the delays need to explicitly propogated.


if mod(nPulses,2)==0
  % Since the Uhrig pulse sequence is symmetric about the  half-way point,
  % only half the delays need to explicitly propogated.
  indices = 0:ceil( (nPulses+1)/2) ;
  delta_t = total_time * diff(sin(indices*pi/( 2*(nPulses+1) ) ).^2);

  % Split the last delat in two parts if needed.
  delta_t(end) = delta_t(end)/2;

  U_ket0 = 1;
  U_ket1 = 1;
  U_bra0 = 1;
  U_bra1 = 1;
  useHbeta = true;
  for t = delta_t
    if useHbeta
      H_ket = Hamiltonian_beta;
      H_bra = Hamiltonian_alpha;
    else
      H_ket = Hamiltonian_alpha;
      H_bra = Hamiltonian_beta;
    end
    useHbeta = ~useHbeta;

    U_prop = propagator_eig(H_ket,t);
    U_ket0 = U_prop*U_ket0;
    U_ket1 = U_ket1*U_prop;

    U_prop = propagator_eig(H_bra,t);
    U_bra0 = U_prop*U_bra0;
    U_bra1 = U_bra1*U_prop;
  end

else

  %   indices = 0:ceil(nPulses+1);
  %   delta_t = total_time * diff(sin(indices*pi/( 2*(nPulses+1) ) ).^2);
  indices = 0:( (nPulses+1)/2) ;
  delta_t = total_time * diff(sin(indices*pi/( 2*(nPulses+1) ) ).^2);

  U_ket0 = 1;
  U_ket1 = 1;
  U_bra0 = 1;
  U_bra1 = 1;

  useHbeta = true;
  for t = delta_t
    if useHbeta
      H_ket = Hamiltonian_beta;
      H_bra = Hamiltonian_alpha;
    else
      H_ket = Hamiltonian_alpha;
      H_bra = Hamiltonian_beta;
    end
    useHbeta = ~useHbeta;

    U_prop = propagator_eig(H_ket,t);
    U_ket0 = U_prop*U_ket0;
    U_bra1 = U_bra1*U_prop;

    U_prop = propagator_eig(H_bra,t);
    U_bra0 = U_prop*U_bra0;
    U_ket1 = U_ket1*U_prop;
  end

end
U_ = (U_bra1*U_bra0)' * (U_ket1*U_ket0);

v_iTime = vecDensityMatrixT*U_(:);

end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>






%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%-------------------------------------------------------------------------------
function subcluster_data = build_subcluster_data(cluster)

cluster_size = numel(cluster);

include_indices = include_vertices_in_which_subclusters(cluster_size);

buff_size = (cluster_size + 2)*2^(cluster_size-1) - cluster_size - 2;
subcluster_data = zeros(1,buff_size);
idx = 1;

for ii = 1:size(include_indices,1)
  subcluster = cluster(include_indices(ii,:));
  subcluster_size = numel(subcluster);

  subcluster_data(idx) = subcluster_size;
  idx = idx + 1;

  subcluster_data(idx:idx+subcluster_size-1) = subcluster;
  idx = idx + subcluster_size;

end


% Check resutls.
expected_size = (cluster_size + 2)*2^(cluster_size-1) - cluster_size - 2;
assert(numel(subcluster_data)==expected_size);
if cluster_size==3
  a = cluster(1); b = cluster(2); c = cluster(3);
  expected_data = [2,b,c, 2,a,c, 1,c, 2,a,b, 1,b, 1,a];
  assert(all(subcluster_data==expected_data))
end

end
%-------------------------------------------------------------------------------
function include_indices = include_vertices_in_which_subclusters(cluster_size)
n_subclusters = 2^cluster_size;

include_indices = false(n_subclusters-2,cluster_size);
increments = 2.^((1:cluster_size)-1);
keep = true(1,cluster_size);
for ii = 2:n_subclusters-1
  sele = mod(ii-1,increments)==0;
  keep(sele) = ~keep(sele);
  include_indices(ii-1,keep) = true;
end
end
%-------------------------------------------------------------------------------
% This function is not thread safe.  Do not use it inside parfor loops.
function auxiliary_signals = update_auxiliary_signals(auxiliary_signals,...
  start_cluster,end_cluster,...
  start_cluster_size,end_cluster_size,cluster_map, clusters)

assert(1 <= start_cluster_size);
assert(start_cluster_size <= end_cluster_size);
assert(end_cluster_size <= numel(clusters));


for cluster_size = start_cluster_size:end_cluster_size
  num_clusters = size(clusters{cluster_size},1);
  start_idx = 1;
  end_idx = num_clusters;
  if cluster_size ==end_cluster_size
    start_idx = start_cluster;
    end_idx = end_cluster;
  end
  for icluster = start_idx:end_idx

    % get_subclusters(cluster) returns all subcluster sizes and subclusters:
    % For the example cluster [a,b,c], build_subcluster_data([a,b,c])
    % returns [2,b,c, 2,a,c, 1,c, 2,a,b, 1,b, 1,a];
    % note that each subcluster is preceded by its size.
    cluster = clusters{cluster_size}(icluster,:);
    subcluster_data = build_subcluster_data(cluster);

    idx = 1;
    while idx < numel(subcluster_data)
      subcluster_size = subcluster_data(idx);
      assert(subcluster_size<cluster_size);
      idx = idx + 1;
      subcluster = zeros(1,subcluster_size);
      for ivertex = 1:subcluster_size
        subcluster(ivertex) = subcluster_data(idx);
        idx = idx + 1;
      end

      %key = num2str(subcluster);
      key = sprintf('%d,',subcluster);
      %key = {subcluster}; % TODO: Faster, but incorect, fix.
      if ~cluster_map.isKey(key)
        %         assert(sum(all(subcluster==clusters{subcluster_size},2))==0);
        continue;
      end



      %       assert(sum(all(subcluster==clusters{subcluster_size},2))==1);

      subcluster_index = cluster_map(key);

      auxiliary_signals{cluster_size}(:,icluster)...
        = auxiliary_signals{cluster_size}(:,icluster)...
        ./auxiliary_signals{subcluster_size}(:,subcluster_index);
    end

  end
end
end
%-------------------------------------------------------------------------------
function cluster_map = build_cluster_map(clusters)

cluster_map = dictionary(string([]),[]);
% cluster_map = dictionary({},[]);  % TODO: Faster, but incorect, fix.

for isize = 1:numel(clusters)
  for ii = 1:size(clusters{isize},1)
    cluster = clusters{isize}(ii,:);
    key = sprintf('%d,',cluster);
    cluster_map(key) = ii;
    %     cluster_map({cluster}) = ii;
  end
end

end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function thisCluster= extract_cluster(cluster_idx,clusterSize,...
  clusters,Nuclei,System)


isMethylCarbon = strcmp(Nuclei.Type,'CH3');

Cluster = clusters{clusterSize}(cluster_idx,:);

if System.Methyl.method ~= 2 && any(isMethylHydron(Cluster) )
  error(['Error in calculateSignal(): ',...
    'A methyl hydron is misplaced.']);
end
if System.Methyl.method == 2 && any(isMethylCarbon(Cluster) )
  error(['Error in calculateSignal(): ',...
    'A methyl carbon is misplaced.']);
end
if(  length(unique(Cluster)) <  length( Cluster )  )
  error(['Error in calculateSignal(): ',...
    'Cluster constains duplicate indices.']);
end

% Check if iCluster is valid.
if Cluster(1,1) == 0
  thisCluster = [];
  return;
end

% Determin the number of spins in the cluster
thisClusterSize = clusterSize + 2*sum(isMethylCarbon(Cluster));

% Initialize cluster
thisCluster = zeros(1,thisClusterSize);

% Initialize counters.
thisIndex = 0;
methyl_number = 0;

% Loop over the cluster to assign methyl rotors.
for ii = 1:clusterSize

  if isMethylCarbon(Cluster(ii))
    % Count the number of metyls.
    methyl_number = methyl_number + 1;

    % Unpack spins from methyl pseudo-particle.
    if Method.useCentralSpinSystem

      hydrons = ...
        nonzeros(Nuclei.associatedParticleMatrix(:,Cluster(ii)+1 ))' -1;
    else
      hydrons = Cluster(ii) + [1,2,3];
    end
    thisCluster(thisIndex+1:thisIndex+numel(hydrons)) = hydrons;
    thisIndex = thisIndex + numel(hydrons);
  else
    thisIndex = thisIndex + 1;
    thisCluster(thisIndex) = Cluster(ii);
  end

end

if(  length(unique(thisCluster)) <  length( thisCluster )  )

  error(['Error in calculateSignal(): ',...
    'This cluster constains duplicate indices.']);
end
if System.Methyl.method ~= 2 && mod(sum(isMethylHydron(thisCluster) ),3) ~= 0
  error(['Error in calculateSignal(): ',...
    'This cluster contains a partial methyl group.']);
end
if any(isMethylCarbon(thisCluster) )
  error(['Error in calculateSignal(): ',...
    'A methyl carbon is misplaced.']);
end
end
%-------------------------------------------------------------------------------
function skipCluster = do_skip_cluster(thisCluster,Nuclei,System)
% Decide if the cluster should be skipped.

if isempty(thisCluster)
  skipCluster = true;
  return;
end

skipCluster =  (System.spinHalfOnly ...
  && Nuclei.StateMultiplicity(thisCluster(1)) > 2);

skipCluster = skipCluster || (~all(Nuclei.StateMultiplicity(thisCluster) ...
  == Nuclei.StateMultiplicity(thisCluster(1))));
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

