function [Alpha, Beta, gridWeight,GridInfo] = getOrientationGrid(System)

gridSize = System.gridSize;
GridInfo = [];
if gridSize==1 && ~strcmp(System.averaging,'custom')
  System.averaging = 'none';
end

if strcmp(System.averaging,'powder')
  
  Grid = lebedev_grid(gridSize);
  % Use the inversion symmetry of the spin-Hamiltonian to skip all the
  
  removeGridPoints = 1:gridSize;
  removeGridPoints(Grid.z < 0) = -1;
  removeGridPoints( ((Grid.z==0) & (Grid.x<0)) ) = -1;
  removeGridPoints( ((Grid.z==0) & (Grid.x==0)) & (Grid.y<0) ) =-1;
  keep = removeGridPoints>0;
  
  Grid.w = Grid.w(keep);
  Grid.w = Grid.w/sum(Grid.w);
  
  Grid.x = Grid.x(keep);
  Grid.y = Grid.y(keep);
  Grid.z = Grid.z(keep);
  
  % Convert xyz coordinates to alpha/beta angles
  Alpha(numel(Grid.z)) = 0;
  Beta(numel(Grid.z)) = 0;
  for iOri = 1:numel(Grid.z)
    
    % Determine the beta Euler angle.
    Beta(iOri) = acos(Grid.z(iOri));
    
    % Determine the alpha Euler angle.
    if (Grid.x(iOri)^2+Grid.y(iOri)^2)>0 % check for the pole singularities.
      
      % Determine which quadrant of the xy-plane the grid point is in.
      quadZ = sign(Grid.y(iOri));
      if abs(quadZ) < 0.1
        quadZ = 1;
      end
      Alpha(iOri) = quadZ*acos(Grid.x(iOri)/sqrt(Grid.x(iOri)^2 + Grid.y(iOri)^2));
    else % set alpha to zero for beta = 0 or pi.
      Alpha(iOri) = 0;
    end
    
  end
  gridSize = length(Alpha);
  gridWeight = Grid.w;
  
  GridInfo.gridSize   = gridSize;
  GridInfo.gridWeight = gridWeight;
  GridInfo.Alpha      = Alpha;
  GridInfo.Beta       = Beta;
  
elseif strcmp(System.averaging,'none')
  
  % Use only the PDB file orientation.
  Alpha = 0;
  Beta = 0;
  gridWeight = 1;
  
elseif strcmp(System.averaging,'xy')
  
  % Average over rotations about the B0 direction.
  Alpha = linspace(0,pi,gridSize+1);
  Alpha(end) = [];
  Beta = ones(1,gridSize)*pi/2;
  gridWeight = ones(gridSize,1)/gridSize;
  
elseif strcmp(System.averaging,'custom')

  if isstruct(System.Grid)
    Alpha = System.Grid.Alpha;
    Beta = System.Grid.Beta;
    gridWeight = System.Grid.gridWeight;  
  end

elseif strcmp(System.averaging,'Nitroxide_Wband_Weights')

  GridSizes = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, ...
         266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,...
         3074, 3470, 3890, 4334, 4802, 5294, 5810];
  gridIndex = find(GridSizes==gridSize);
  
  load([Data.path2CluE, 'grids/Lebedev_weighted_Nitroxide_Wband_Weights.mat'],'Grids');
  Alpha = Grids{gridIndex}.Alpha;
  Beta = Grids{gridIndex}.Beta;
  gridWeight = Grids{gridIndex}.Weight;

  gridSize = length(Alpha);
  
  GridInfo.gridSize   = gridSize;
  GridInfo.gridWeight = gridWeight;
  GridInfo.Alpha      = Alpha;
  GridInfo.Beta       = Beta;
  
elseif strcmp(System.averaging,'random')
  
  % Generate random Euler angles.
  Alpha = rand(1,gridSize)*2*pi;
  Beta = acos(2*rand(1,gridSize)-1);
  gridWeight = ones(gridSize,1)/gridSize;

end



end