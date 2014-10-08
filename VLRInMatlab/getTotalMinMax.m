function minMax = getTotalMinMax( dataName )
% get the min and max values for each data set according to
% the last time step of their real data points (including contour points)
if strcmp( dataName, '120830_raw' )
    minMax = [ ];
  elseif strcmp( dataName, '121204_raw_2014' )
    minMax = [ ];
  elseif strcmp( dataName, '121211_raw' )
    minMax = [ -11 601 -91 131 ];
  elseif strcmp( dataName, '130508_raw' )
    minMax = [ ];
  elseif strcmp( dataName, '130607_raw' )
    minMax = [ ];
  elseif strcmp( dataName, '131203_raw' )
    minMax = [ ];
  end