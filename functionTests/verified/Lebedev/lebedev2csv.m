

function lebedev2csv() 
availableSizes = [...
     6,   14,   26,   38,   50,   74,   86,  110,  146,  170, ...
   194,  230,  266,  302,  350,  434,  590,  770,  974, 1202, ...
  1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, ...
  5294, 5810];


for isize=availableSizes
  
  lg = lebedev_grid(isize);
  lg.w = lg.w/sum(lg.w);

  filename = ['Lebedev_Grid_', num2str(isize), '.csv'];
  T = array2table([lg.x,lg.y,lg.z,lg.w]);
  T.Properties.VariableNames(1:4) = {'x','y','z','w'};
  writetable(T,filename);
end

exit();
end

