function minMax = getTotalMinMax( dataName, auto )
% get the min and max values for each data set according to
% the last time step of their real data points (including contour points)
% automatic value based on distance of cellDist = 25
if auto == 1
  if strcmp( dataName, '120830_raw' )
    minMax = [ 86.7711157299174 529.769701967227 -101.938395705689 55.9763331792212 ];
  elseif strcmp( dataName, '121204_raw_2014' )
    minMax = [ 112.469987230284 570.320702034532 -270.050023831757 -106.278711320592 ];
  elseif strcmp( dataName, '121211_raw' )
    minMax = [ -13.6716366404918 584.843509713204 -104.231971482375 117.236014022419 ];
  elseif strcmp( dataName, '130508_raw' )
    minMax = [ 24.6657222062463 648.628510237614 197.809309709361 334.809605754869 ];
  elseif strcmp( dataName, '130607_raw' )
    minMax = [ -173.735362370687 405.659180401568 135.735707809344 376.906008812735 ];
  end
else
  if strcmp( dataName, '120830_raw' )
    minMax = [ 99 521 -101 66 ];
  elseif strcmp( dataName, '121204_raw_2014' )
    minMax = [ 102 580 -270 -95 ];
  elseif strcmp( dataName, '121211_raw' )
    minMax = [ -11 601 -91 131 ];
  elseif strcmp( dataName, '130508_raw' )
    minMax = [ 20 405 195 355 ];
  elseif strcmp( dataName, '130607_raw' )
    minMax = [ 102 580 -270 -95 ];
  end
end