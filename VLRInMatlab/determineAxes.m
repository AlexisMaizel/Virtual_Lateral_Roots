function semiAxes = determineAxes( X, Y, Z, center, dir )
  semiLength = 100000;
  maxLength = 0;
  maxSemiAxis = [];
  minSemiAxis = [];
  
  dim = size( X, 1 );
  for d=1:dim
    for c=1:dim
      pos = [ X(c,d) Y(c,d) Z(c,d) ];
      % compute length between current point and center
      length = norm( pos - center );
      
      if length > maxLength
        maxLength = length;
        maxSemiAxis = pos;
      end      
    end
  end
  
  % get cross product to determine the smaller semi axis
  C = cross( dir, maxSemiAxis - center );
  offset = 2;
  C = C * offset;
  
  for d=1:dim
    for c=1:dim
      % compute length between current point and center
      point = [ X(c,d) Y(c,d) Z(c,d) ];
      
      if distancePoints3d( point, center + C ) <= semiLength
        semiLength = distancePoints3d( point, center + C );
        minSemiAxis = [ X(c,d) Y(c,d) Z(c,d) ];
      end
    end
  end
  
  semiAxes = [ minSemiAxis maxSemiAxis ];