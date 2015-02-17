function transMat = getDataSetTranslationMatrix( dataStr, registerBase )
if registerBase == 1
  if strcmp( dataStr, '120830_raw' )
    translation = [ 0 0 0 ];
  elseif strcmp( dataStr, '121204_raw_2014' )
    translation = [ 15 0 0 ];
  elseif strcmp( dataStr, '121211_raw' )
    translation = [ 9 0 0 ];
  elseif strcmp( dataStr, '130508_raw' )
    translation = [ 20 5 0 ];
  elseif strcmp( dataStr, '130607_raw' )
    translation = [ 8 9 0 ];
  end
else
  if strcmp( dataStr, '120830_raw' )
    translation = [ 0 0 0 ];
  elseif strcmp( dataStr, '121204_raw_2014' )
    translation = [ 10 -15 0 ];
  elseif strcmp( dataStr, '121211_raw' )
    translation = [ 12 -5 0 ];
  elseif strcmp( dataStr, '130508_raw' )
    translation = [ 9 -2 0 ];
  elseif strcmp( dataStr, '130607_raw' )
    translation = [ 5 -20 0 ];
  end
end
transMat = createTranslation3d( translation(1,1), translation(1,2), translation(1,3) );
