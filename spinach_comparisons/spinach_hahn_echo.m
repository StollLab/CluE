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
 
  inter.zeeman.matrix = cell(1,N+1);
  inter.coupling.matrix = cell(N+1,N+1);
  inter.zeeman.matrix{1} = diag(System.g);

  for ispin = 1:N

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
  bas.approximation = 'IK-0';
  bas.level = 2;
  bas.connectivity='full_tensors';
%   sys.tols.prox_cutoff=5.0;
  % Generate spin system.
  spin_system=create(sys,inter);
  spin_system=basis(spin_system,bas);
  spin_system=assume(spin_system,'esr');
  

  % Hamiltonian -------------------------------------------------------------
  %    H     - Hamiltonian matrix, received from context function
  [I,Q]=hamiltonian(spin_system);
   
  params.spins = {'E'};
  
  % Apply the offsets
  params.decouple={};
  params.rframes={};
  params.verbose=0;
  params.offset=zeros(size(params.spins));
  I=frqoffset(spin_system,I,params);
  
  % Build the Hamiltonian and tidy up
  params.orientation = [0,0,0];
  H=I+orientation(Q,params.orientation);
  H=(H+H')/2; 
  L = H;
  
  params.npoints = System.timepoints;
  params.timestep = System.dt;
  params.rho0=state(spin_system,'Lz','E');
  params.coil = state(spin_system,'L+','E');
  
  params.screen = state(spin_system,'L-','E');
  params.pulse_op=(operator(spin_system,'L+','E')-operator(spin_system,'L-','E'))/2i;
  params.zerofill=4096;
  
  
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

  
  path(oldpath);
  
  toc
catch
  
  path(oldpath);

end
end