function [Q, curS] = updateBezierSurface( tileGridDir, curS,...
  totalMinAxes, totalMaxAxes, resGrid, rows, columns,...
  magnitudeScaling, dimCP, steps, interpolatedHeightGrowth )
oldVersion = 1;
attenuation = 2;
affectNeighborhood = 0;
if oldVersion == 1
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
          if averageDirection(1,2) > 0.
            translation = [ averageDirection(1,1)...
            averageDirection(1,2)*magnitudeScaling ];
          else
            translation = [ averageDirection(1,1)...
            averageDirection(1,2) ];
          end
          %translation = [ averageDirection(1,1)*magnitudeScaling...
            %averageDirection(1,2)*magnitudeScaling ];
        end
        
        if interpolatedHeightGrowth == 1
          minV = 1;
          maxV = dimCP;
          % translate the inner control point based on the growing
          % height which is determined by interpolation between the uppermost
          % and the lowermost cp
          factor = (i-minV)/(maxV-minV);
          curS( i, j, 2 ) =...
            (1-factor) * curS( minV, j, 2 ) + factor * curS( maxV, j, 2 );
        end
        
        % then apply transformation to cp
        curS( i, j, 1 ) = curS( i, j, 1 ) + translation(1,1);
        curS( i, j, 2 ) = curS( i, j, 2 ) + translation(1,2);
      end
    end
  end
  % interpolation of surface
  Q = bezierSurfaceInterp( curS, steps );
else
  for i=1:dimCP
    for j=1:dimCP
      % get position of cp
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
      
      % now translate each cp radially by decreasing factor with center of
      % the current cp
      if affectNeighborhood == 1
        for ii=1:dimCP
          for jj=1:dimCP
            distGrid = max( abs( ii - i ), abs( jj - j ) );
            curTrans = translation * 1./power( attenuation, distGrid );
            % then apply transformation to cp
            for c=1:2
              curS( ii, jj, c ) = curS( ii, jj, c ) + curTrans(1,c);
            end
          end
        end
      else
        curS( i, j, 1 ) = curS( i, j, 1 ) + translation(1);
        curS( i, j, 2 ) = curS( i, j, 2 ) + translation(2);
      end
    end
  end
  % interpolation of surface
  Q = bezierSurfaceInterp( curS, steps );
end