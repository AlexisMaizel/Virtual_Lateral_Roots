function rotMat = getManualRotationMatrix( dataStr )
if strcmp( dataStr, '120830_raw' )
  angle = 0;
elseif strcmp( dataStr, '121204_raw_2014' )
  angle = 0;
elseif strcmp( dataStr, '121211_raw' )
  angle = 0;
elseif strcmp( dataStr, '130508_raw' )
  angle = -2;
elseif strcmp( dataStr, '130607_raw' )
  angle = 2;
elseif strcmp( dataStr, '131203_raw' )
  angle = 0;
end
rad = degtorad(angle);
rotMat = createRotationOz(rad);