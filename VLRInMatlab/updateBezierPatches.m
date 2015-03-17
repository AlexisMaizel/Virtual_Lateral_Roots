function [Q, curS] = updateBezierPatches( tileGridDir, curS,...
  totalMinAxes, totalMaxAxes, resGrid, rows, columns,...
  magnitudeScaling, numPatches, interpolatedHeighGrowth )
%number of interpolated values between end control points
ni = 10;
%uniform parameterization
u = linspace(0,1,ni);
v = u;

if numPatches == 1
  for i=1:4
    for j=1:4
      % only interpolate inner control points (cp)
      if i == 1 || i == 4 || j == 1 || j == 4
        % do nothing for the boundary points
      else
        % position of cp
        cp = [ curS( i, j, 1, 1 ) curS( i, j, 2, 1 ) ];
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
        
        % translate the inner control point based on the growing
        % height which is determined by interpolation between the uppermost
        % and the lowermost cp
        if interpolatedHeighGrowth == 1
          min = 1;
          max = 4;
          factor = (i-min)/(max-min);
          curS( i, j, 2, 1 ) =...
            (1-factor) * curS( 1, j, 2, 1 ) + factor * curS( 4, j, 2, 1 );
        end
        % then apply transformation to cp
        curS( i, j, 1, 1 ) = curS( i, j, 1, 1 ) + translation(1,1);
        curS( i, j, 2, 1 ) = curS( i, j, 2, 1 ) + translation(1,2);
      end
    end
  end
  % interpolation of patch
  Q(:,:,:,1) = bezierpatchinterp(curS(:,:,:,1), u, v);
elseif numPatches == 4
  % loop over patches
  for p=1:numPatches
    % loop over cp
    for i=1:4
      for j=1:4
        % only interpolate inner cp
        if (i == 1 && (p == 1 || p == 2)) ||...
            (i == 4 && (p == 3 || p == 4)) ||...
            (j == 1 && (p == 1 || p == 3)) ||...
            (j == 4 && (p == 2 || p == 4))
          % do nothing for the boundary points
        else
          % position of cp
          cp = [ curS( i, j, 1, p ) curS( i, j, 2, p ) ];
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
            max = 4;
            % translate the inner control point based on the growing
            % height which is determined by interpolation between the uppermost
            % and the lowermost cp
            if p == 1
              factor = (i-min)/(max-min);
              % divide by two because there are 2 patches
              factor = factor/2.;
              curS( i, j, 2, p ) =...
                (1-factor) * curS( 1, j, 2, 1 ) + factor * curS( 4, j, 2, 3 );
            elseif p == 2
              factor = (i-min)/(max-min);
              % divide by two because there are 2 patches
              factor = factor/2.;
              curS( i, j, 2, p ) =...
                (1-factor) * curS( 1, j, 2, 2 ) + factor * curS( 4, j, 2, 4 );
            elseif p == 3
              factor = (i+2)/(max-min);
              % divide by two because there are 2 patches
              factor = factor/2.;
              curS( i, j, 2, p ) =...
                (1-factor) * curS( 1, j, 2, 1 ) + factor * curS( 4, j, 2, 3 );
            else
              factor = (i+2)/(max-min);
              % divide by two because there are 2 patches
              factor = factor/2.;
              curS( i, j, 2, p ) =...
                (1-factor) * curS( 1, j, 2, 2 ) + factor * curS( 4, j, 2, 4 );
            end
          end
          
          % then apply transformation to cp
          curS( i, j, 1, p ) = curS( i, j, 1, p ) + translation(1,1);
          curS( i, j, 2, p ) = curS( i, j, 2, p ) + translation(1,2);
        end
      end
    end
    % interpolation of patch
    Q(:,:,:,p) = bezierpatchinterp(curS(:,:,:,p), u, v);
  end
else
  disp( strcat( 'Number of patches is not supported yet! ', num2str(numPatches) ) );
end