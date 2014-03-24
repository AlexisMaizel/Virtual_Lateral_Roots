function transMat = getManualTranslationMatrix( dataStr )
if strcmp( dataStr, '120830_raw' )
  translation = [ -29 9 0 ];
elseif strcmp( dataStr, '121204_raw_2014' )
  translation = [ 0 0 0 ];
elseif strcmp( dataStr, '121211_raw' )
  translation = [ -5 5 0 ];
elseif strcmp( dataStr, '130508_raw' )
  translation = [ 0 12 0 ];
elseif strcmp( dataStr, '130607_raw' )
  translation = [ -30 3 0 ];
elseif strcmp( dataStr, '131203_raw' )
  translation = [ 0 0 0 ];
end
transMat = createTranslation3d( translation(1,1), translation(1,2), translation(1,3) );

