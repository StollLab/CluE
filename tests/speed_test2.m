% This function runs an example time propagation multiple times, and return
% the mean time.

function [t_mean, delta_t] = speed_test2(d)
% Define beta and alpha Hamilonians.
% d = 4;

% Run trials.
N_trials = 10;
t_mean = 0;
delta_t = zeros(1,N_trials);
for itrials = 1:N_trials
  tic;
  
  Hamiltonian_beta  = rand(d) + 1i*rand(d);
  
  Hamiltonian_beta = (Hamiltonian_beta+Hamiltonian_beta')/2;
  
  Hamiltonian_alpha = rand(d) + 1i*rand(d);
  
  Hamiltonian_alpha = (Hamiltonian_alpha+Hamiltonian_alpha')/2;
  
  
  % Define density matrix in the high temperature limit.
  DensityMatrix = eye(d)/d;
  vecDensityMatrixT = reshape(DensityMatrix.',1,[])/trace(DensityMatrix);
  
  % Define time parameters.
  dt = 0.04*1e-6; % us.
  timepoints = 2^10;
  
  
  
  % Initialize propagators.
  U_beta = eye(d);
  U_alpha = eye(d);
  dU_beta = propagator_eig(Hamiltonian_beta,dt);
  dU_alpha = propagator_eig(Hamiltonian_alpha,dt);
  
  % Initialize signal.
  v= ones(1 ,timepoints);
  
  % Loop over time points.
  for iTime = 1:timepoints
    
    U_ = U_beta'*U_alpha'*U_beta*U_alpha;
    %     v(iTime) = vecDensityMatrixT*U_(:);
    v(iTime) =trace(DensityMatrix*U_);
    U_beta = dU_beta*U_beta;
    U_alpha = dU_alpha*U_alpha;
  end
  v = v./v(1);
  
  delta_t(itrials) = toc;
  t_mean = t_mean + delta_t(itrials);
end

t_mean = t_mean/N_trials;
fprintf('t_mean = %d μs.\n',t_mean*1e6);
fprintf('t_mean = %d ms.\n',1000*t_mean);
fprintf('t_mean = %d s.\n',t_mean);
end
% ========================================================================
% Calculate propagator using diagonalization
% ========================================================================
function v = runPropagation(Hamiltonian_alpha,Hamiltonian_beta,dt,timepoints,vecDensityMatrixT)
% Initialize propagators.
U_beta = eye(4);
U_alpha = eye(4);
dU_beta = propagator_eig(Hamiltonian_beta,dt);
dU_alpha = propagator_eig(Hamiltonian_alpha,dt);

% Initialize signal.
v= ones(1 ,timepoints);

% Loop over time points.
for iTime = 1:timepoints
  
  U_ = U_beta'*U_alpha'*U_beta*U_alpha;
  v(iTime) = vecDensityMatrixT*U_(:);
  
  U_beta = dU_beta*U_beta;
  U_alpha = dU_alpha*U_alpha;
end
v = v./v(1);
end
function U = propagator_eig(Ham,t)
%Ham = (Ham+Ham')/2; % "hermitianize" Hamiltonian
% hbar = 1.054571800e-34;
[EigenVectors, EigenValues] = eig(Ham);
Udiag = exp(-2i*pi*diag(EigenValues)*t);
U = EigenVectors*diag(Udiag)*EigenVectors';
end
