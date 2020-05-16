% This function runs an example time propagation multiple times, and return
% the mean time.

function [t_mean, v,twotau] = speed_test()
% Define beta and alpha Hamilonians.

Hb = 1.0e+10 * [ ...
  -1.6868         0         0         0; ...
  0   -1.6816         0         0; ...
  0         0   -1.6816         0; ...
  0         0         0   -1.6765];

Hamiltonian_beta  = Hb + 1.0e+03 * [ ...
  0.0000 + 0.0000i  -1.0295 - 0.0391i  -1.0295 - 0.0391i  -0.5563 - 0.0423i; ...
  -1.0295 + 0.0391i   0.0000 + 0.0000i   1.0825 + 0.0000i   1.0295 + 0.0391i; ...
  -1.0295 + 0.0391i   1.0825 + 0.0000i   0.0000 + 0.0000i   1.0295 + 0.0391i; ...
  -0.5563 + 0.0423i   1.0295 - 0.0391i   1.0295 - 0.0391i   0.0000 + 0.0000i];

Hamiltonian_beta = (Hamiltonian_beta+Hamiltonian_beta')/2;

Ha = 1.0e+10 * [ ...
  1.6765         0         0         0; ...
  0    1.6816         0         0; ...
  0         0    1.6816         0; ...
  0         0         0    1.6868];

Hamiltonian_alpha = Ha + 1.0e+03 * [...
  0.0000 + 0.0000i  -1.0295 - 0.0391i  -1.0295 - 0.0391i  -0.5563 - 0.0423i; ...
  -1.0295 + 0.0391i   0.0000 + 0.0000i   1.0825 + 0.0000i   1.0295 + 0.0391i; ...
  -1.0295 + 0.0391i   1.0825 + 0.0000i   0.0000 + 0.0000i   1.0295 + 0.0391i; ...
  -0.5563 + 0.0423i   1.0295 - 0.0391i   1.0295 - 0.0391i   0.0000 + 0.0000i];

Hamiltonian_alpha = (Hamiltonian_alpha+Hamiltonian_alpha')/2;

[Xb, Eb] = eig(Hamiltonian_beta);
Xb_inv = inv(Xb);
disp('Hb')
disp(Hamiltonian_beta)
disp('Eb')
disp(Eb)
disp('Xb')
disp(Xb)
disp('Hb -Xb*Eb*Xb_inv')
disp(Hb -Xb*Eb*Xb_inv)


% Define density matrix in the high temperature limit. 
DensityMatrix = eye(4)/4;
vecDensityMatrixT = reshape(DensityMatrix.',1,[])/trace(DensityMatrix);

% Define time parameters.
dt = 0.04*1e-6; % us.
timepoints = 2^7;

twotau = 2*dt*[0:timepoints-1];


% Run trials.
N_trials = 1;
t_mean = 0;
tic;
for itrial = 1:N_trials
  
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
%     v(iTime) = vecDensityMatrixT*U_(:);
    v(iTime) =trace(DensityMatrix*U_);
    U_beta = dU_beta*U_beta;
    U_alpha = dU_alpha*U_alpha;
  end
  v = v./v(1);
end

t_mean = toc;
t_mean = t_mean/N_trials;
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
