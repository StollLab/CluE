% This function creates a a set of matrices (each member correspond to a particular pairwise coupling) for each nuclei for each electronic state.
% The full Hamiltonian (not used) would be the direct product over all nuclei and
% electronic states of the sum over all couplings.  
% A cluster Hamiltonian, which is used, is the direct product over all
% nuclei in the cluster of the sum over all couplings in the cluster.

function [Hamiltonian,zeroIndex] = pairwiseHamiltonian(System,Nuclei,Cluster)

%   methyl_IDs = zeros(size(Cluster));
  Cluster_ = Cluster;
  for ispin = Cluster_
    if strcmp('CH3', Nuclei.Type{ispin})
      %       methyl_IDs = [methyl_IDs,  Nuclei.MethylID(ispin)];
      %       Cluster = [Cluster, Nuclei.Auxiliary_ID(ispin,:)]
      Cluster = [Cluster, [ispin+1,ispin+2,ispin+3]];
      Cluster(Cluster==ispin)=[];
      %       Nuclei.Coordinates(Nuclei.Auxiliary_ID(ispin,:),:) = Nuclei.Auxiliary_Coordinates{ispin};
    end
  end
  clear Cluster_;
  Cluster = sort(unique(Cluster));
  
%   if ~any(methyl_IDs>0)
%     getMethylHamiltonian(System,Nuclei,Cluster,methyl_IDs)
%   end
  zeroIndex = min(Cluster) - 1;
  Indices = fliplr(Cluster);
  N = size(Cluster,2);
  Hamiltonian = cell(N+1,N+1); % nspins by nspins
  
  % Electron Zeeman
  Hamiltonian{1,1} =  constructElectronZeeman(System);
  for i_index_nucleus = Indices
    inucleus = i_index_nucleus - zeroIndex;
    % Nuclear Zeeman
    Hamiltonian{1+inucleus,1+inucleus} =   constructNuclearZeeman(System, Nuclei,i_index_nucleus);
 
    % Hyperfine
    Hamiltonian{1,1+inucleus}= constructHyperfine(System, Nuclei, i_index_nucleus);
   
    % Nucleus-Nucleus Coupling
    for j_index_nucleus = Indices
      jnucleus = j_index_nucleus - zeroIndex;
      if jnucleus <= inucleus
        continue;
      end
 
      % Dipole Coupling
      Hamiltonian{1+inucleus,1+jnucleus} = constructNuclearDipoleCoupling(System, Nuclei,i_index_nucleus,j_index_nucleus);      
    end
  end
   
end

function electronZeeman = constructElectronZeeman(System) 
  ge = System.gMatrix(3,3);
  electronZeeman = ge*System.muB*System.magneticField; % J.
  electronZeeman = electronZeeman/(2*pi*System.hbar); % Hz. 
end

function NuclearZeeman = constructNuclearZeeman(System, Nuclei, i_index_nucleus)
   
    NuclearZeeman = Nuclei.Nuclear_g(i_index_nucleus)*System.magneticField*System.muN*eye(3); % J.
    NuclearZeeman = NuclearZeeman/(2*pi*System.hbar); % Hz
  
end

function Hyperfine = constructHyperfine(System, Nuclei, i_index_nucleus)

  km = System.mu0/(4*pi);
  muB = System.muB;
  muN = System.muN;
  if strcmp(Nuclei.Type{i_index_nucleus},'e')
    muN = -muB;
  end
  ge = System.Electron.g;
  gni = Nuclei.Nuclear_g(i_index_nucleus);
  r = Nuclei.Coordinates(i_index_nucleus,:);
  if size(r,2)==3
    r=r';
  end
  n = r/norm(r);
  nnt = n*n'; 
  r3 = norm(r)^3;
  Hdd = km*ge*muB*gni*muN/r3*(eye(3)-3*nnt);
  Hyperfine = -Hdd/(2*pi*System.hbar); % Hz.

  return;

end

function Hdd = constructNuclearDipoleCoupling(System, Nuclei,i_index_nucleus,j_index_nucleus)
  km = System.mu0/(4*pi);
  muNi = System.muN;
  if strcmp(Nuclei.Type{i_index_nucleus},'e')
    muNi = -muB;
  end
  muNj = System.muN;
  if strcmp(Nuclei.Type{j_index_nucleus},'e')
    muNj = -muB;
  end

  gni = Nuclei.Nuclear_g(i_index_nucleus);
  gnj = Nuclei.Nuclear_g(j_index_nucleus);
  r = Nuclei.Coordinates(i_index_nucleus,:)-Nuclei.Coordinates(j_index_nucleus,:);
  if size(r,2)==3
    r=r';
  end
  n = r/norm(r);
  nnt = n*n'; 
  r3 = norm(r)^3;
  Hdd = -km*gni*gnj*muNi*muNj/r3*(eye(3)-3*nnt);
  Hdd = Hdd/(2*pi*System.hbar); % Hz.
  
end

function cos_theta = cosFieldAngle(ParticleCoordinates_1,ParticleCoordinates_2)
   R = ParticleCoordinates_2 - ParticleCoordinates_1;
   cos_theta = R*[0;0;1]/norm(R); 
end


function  H  = getMethylHamiltonian(System,Nuclei,Cluster,methyl_IDs)
nonMethyls = Cluster(methyl_IDs == 0);
Methyls = Cluster(methyl_IDs > 0);

if length(Methyls) > length(unique(Methyls))
  return;
end
N = length(nonMethyls) + 3*length(Methyls);
Cluster_ =1:N;


Methyl_Spin.Type = {'1H','1H','1H'};
for imethyl = Methyls
 Methyl_Spin.Nuclear_g = Nuclei.Nuclear_g(imethyl)*[1,1,1];
 Methyl_Spin.Coordinates = Nuclei.Auxiliary_Coordinates{imethyl};
 [H_in,zeroIndex] = pairwiseHamiltonian(System,Methyl_Spin,Cluster);
 H_out = assembleHamiltonian(H_in,Cluster,System, eState,Nuclei,zeroIndex,clusterSize);
 H_out = Nuclei.Projection(imethyl)*Nuclei.Transform{imethyl}*H_out*Nuclei.Projection(imethyl);
 
end


end
