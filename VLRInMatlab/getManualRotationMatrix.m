function rotMat = getManualRotationMatrix( dataStr, renderMasterFile )
if renderMasterFile == 1
  if strcmp( dataStr, '120830_raw' )
    angle = 0;
  elseif strcmp( dataStr, '121204_raw_2014' )
    angle = 0;
  elseif strcmp( dataStr, '121211_raw' )
    angle = 0;
  elseif strcmp( dataStr, '130508_raw' )
    angle = -3;
  elseif strcmp( dataStr, '130607_raw' )
    angle = 4;
  elseif strcmp( dataStr, '131203_raw' )
    angle = 0;
  end
else
  if strcmp( dataStr, '120830_raw' )
    angle = -4;
  elseif strcmp( dataStr, '121204_raw_2014' )
    angle = 2;
  elseif strcmp( dataStr, '121211_raw' )
    angle = 0;
  elseif strcmp( dataStr, '130508_raw' )
    angle = -5;
  elseif strcmp( dataStr, '130607_raw' )
    angle = 4;
  elseif strcmp( dataStr, '131203_raw' )
    angle = 0;
  end
end
rad = degTorad(angle);
rotMat = createRotationOz(rad);