function [Q, curS] = updateBezierSurface( tileGridDir, curS,...
  totalMinAxes, totalMaxAxes, resGrid, rows, columns,...
  magnitudeScaling, dimCP, steps, interpolatedHeighGrowth )
for i=1:dimCP
  for j=1:dimCP
    % only interpolate inner cp
    if i == 1 || i == dimCP || j == 1 || j == dimCP
      % do nothing for the boundary points
    else
      % position of cp
      cp = [ curS( i, j, 1 ) curS( i, j, 2 ) ];
      translation = [ 0 0 ];
      % get tile index of cp in grid
      tileIndex = getTileIndex( cp, [totalMinAxes(1) totalMinAxes(2)],...
        [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns );
      % ignore empty tiles
      if size( tileGridDir{tileIndex}, 1 ) ~= 0
        [ averageDirection ] =...
          determineAverageDirection( tileGridDir{tileIndex} );
        translation = [ averageDirection(1,1)*magnitudeScaling...
          averageDirection(1,2)*magnitudeScaling ];
      end
      
      if interpolatedHeighGrowth == 1
        min = 1;
        max = dimCP;
        % translate the inner control point based on the growing
        % height which is determined by interpolation between the uppermost
        % and the lowermost cp
        factor = (i-min)/(max-min);
        curS( i, j, 2 ) =...
          (1-factor) * curS( min, j, 2 ) + factor * curS( max, j, 2 );
      end
      
      % then apply transformation to cp
      curS( i, j, 1 ) = curS( i, j, 1 ) + translation(1,1);
      curS( i, j, 2 ) = curS( i, j, 2 ) + translation(1,2);
    end
  end
end
% interpolation of surface
Q = bezierSurfaceInterp( curS, steps );