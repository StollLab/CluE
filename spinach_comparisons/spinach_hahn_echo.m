function fid = spinach_hahn_echo(System,Data,path2spinach)

% Get system.
[System, Tensors, Nuclei,Clusters] = setUpSystem(System,Data);
N = Nuclei.number;

oldpath = path;
fh = fopen('spinach_paths.txt');
allLines = textscan(fh,'%s','whitespace','','delimiter','\n');
fclose(fh);
allLines = allLines{1};

nLines = numel(allLines);
path_ = oldpath;
for iline = 1:nLines
  line_ = allLines{iline};
  path([path2spinach, line_],path_);
  path_ = path;
end

try
  tic
  %  sys     - spin system and instrument specification
  %            structure, see the spin system specification
  %            section of the online manual
  
  % Magnet field
  sys.magnet=System.magneticField;
  sys.isotopes = cell(1,N+1);
  sys.isotopes{1} = 'E';
  %  inter   - interaction specification structure, see
  %            see the spin system specification section
  %            of the online manual
  % inter.coordinates = cell(N,1);
  inter.zeeman.matrix = cell(1,N+1);
  inter.coupling.matrix = cell(N+1,N+1);
  inter.zeeman.matrix{1} = diag(System.g);
%   inter.coordinates{1} = [0,0,0];
  for ispin = 1:N
%       inter.coordinates{ispin+1}=Nuclei.Coordinates(ispin,:)*1e10;
    sys.isotopes{ispin + 1} = Nuclei.Type{ispin};
    
    if System.nuclear_quadrupole
      inter.coupling.matrix{ispin+1,ispin+1} = Nuclei.Qtensor(:,:,ispin);
    end
    for jspin = ispin+1:N+1
      inter.coupling.matrix{ispin,jspin} = Tensors{ispin,jspin};
    end
    
  end
  
  
  
  % Basis set
  bas.formalism='sphten-liouv';
  bas.approximation = 'none';
  bas.connectivity='full_tensors';
  
  % Generate spin system.
  spin_system=create(sys,inter);
  spin_system=basis(spin_system,bas);
  spin_system=assume(spin_system,'esr');
  
  % parameters --------------------------------------------------------------
  % L      - the Liouvillian to be used during evolution
  %
  %      rho    - the initial state vector or a horizontal stack thereof
  %
  %      output - a string giving the type of evolution that is required
  %
  %                'final' - returns the final state vector or a horizontal
  %                          stack thereof.
  %
  %                'trajectory' - returns the stack of state vectors giving
  %                               the trajectory of the system starting from
  %                               rho with the user-specified number of steps
  %                               and step length.
  %
  %                'total'   - returns the integral of the observable trace
  %                            from the simulation start to infinity. This
  %                            option requires the presence of relaxation.
  %
  %                'refocus' - evolves the first vector for zero steps,
  %                            second vector for one step, third vector for
  %                            two steps, etc., consistent with the second
  %                            stage of evolution in the indirect dimension
  %                            after a refocusing pulse.
  %
  %                'observable' - returns the time dynamics of an observable
  %                               as a vector (if starting from a single ini-
  %                               tial state) or a matrix (if starting from a
  %                               stack of initial states).
  %
  %                'multichannel' - returns the time dynamics of several
  %                                 observables as rows of a matrix. Note
  %                                 that destination state screening may be
  %                                 less efficient when there are multiple
  %                                 destinations to screen against.
  %
  %      coil   - the detection state, used when 'observable' is specified as
  %               the output option. If 'multichannel' is selected, the coil
  %               should contain multiple columns corresponding to individual
  %               observable vectors.
  %
  %      destination - (optional) the state to be used for destination state
  %                    screening.                                  acquisition begins (optional)
  
  % Hamiltonian -------------------------------------------------------------
  %    H     - Hamiltonian matrix, received from context function
  [I,Q]=hamiltonian(spin_system);
   
  params.spins = {'E'};
  % Apply the offsets
  params.decouple={};
  params.rframes={};
  params.verbose=1;
  params.offset=zeros(size(params.spins));
  I=frqoffset(spin_system,I,params);
  
  % Build the Hamiltonian and tidy up
  params.orientation = [0,0,0];
  H=I+orientation(Q,params.orientation);
  H=(H+H')/2; %clear('I','Q');
  L = H;
  
  params.npoints = System.timepoints;
  params.timestep = System.dt;
  params.rho0=state(spin_system,'Lz','E');
%   params.rho0 = equilibrium(spin_system,H,Q,[0,0,0]);
  params.coil = state(spin_system,'L+','E');
%   params.pulse_op = operator(spin_system,'L+','E','left');
  
  params.screen = state(spin_system,'L-','E');
  params.pulse_op=(operator(spin_system,'L+','E')-operator(spin_system,'L-','E'))/2i;
  params.zerofill=4096;
  
 % %{
  % params.rho0 = equilibrium(spin_system,H,Q,[0,0,0]);
  
  % First pulse
  rho = step(spin_system,params.pulse_op,params.rho0,pi/2);
  
  % Free evolution
  rho_stack = evolution(spin_system,L,[],rho,params.timestep,...
    (params.npoints-1),'trajectory');
  
  % Second pulse
  rho_stack=step(spin_system,params.pulse_op,rho_stack,pi);
  
  % Free evolution
  rho_stack=evolution(spin_system,L,[],rho_stack,params.timestep,...
    (params.npoints-1),'refocus',params.coil);
  
  % Detect
  fid=transpose(full(params.coil'*rho_stack));
  %}
  
% Make pulse operators
%{
Hp=operator(spin_system,'L+',params.spins{1});
% Hy=kron(speye(params.npts),(Hp-Hp')/2i);
Hy=params.pulse_op;
% Hard 90-degree pulse
rho=step(spin_system,Hy,params.rho0,pi/2);

% Evolution under the X gradient
rho=evolution(spin_system,L,[],rho,params.timestep,params.npoints,'final');

% Hard 180-degree pulse
rho=step(spin_system,Hy,rho,pi);

% Detection under the X gradient
fid=evolution(spin_system,L,params.coil,rho,...
  params.timestep, 2*params.npoints,'observable');
%}

% CPMG echo train with detection. Syntax:
%
%          fid=cpmg(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.nloops   - number of CPMG loops
%
%    parameters.timestep - time step
%
%    parameters.npoints  - number of steps per half-echo
%
% Outputs:
%
%    fid - free induction decay throughout the sequence
%
% i.kuprov@soton.ac.uk
%
% <http://spindynamics.org/wiki/index.php?title=cpmg.m>
% params.nloops = 1;  
% params.spc_dim = 1;
% 
% fid=cpmg(spin_system,params,H,0*H,0*H);

  
  
  if fid(1) ~= 0
%     fid =fid./fid(1);
  end
  path(oldpath);
  toc
catch
  path(oldpath);
end
end