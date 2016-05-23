function setSubPlots( f, numPlots, chosenData, spatialNormalization, width, height )
% subplot settings
clf(f)
if chosenData < 6
  for p=1:numPlots
    subplot( numPlots/2, numPlots/2, p );
    hold on
    if spatialNormalization == 0
      xMinMax = [ 0 width ];
      if p == numPlots
        yMinMax = [ 0 400 ];
      else
        yMinMax = [ 0 height ];
      end
    else
      xMinMax = [ -width/2 width/2 ];
      yMinMax = [ -height/2 height/2 ];
    end
    % axis([xmin xmax ymin ymax zmin zmax cmin cmax])
    axis( [ xMinMax(1) xMinMax(2) yMinMax(1) yMinMax(2) -10000 10000 0 1 ] );
    axis on
    daspect( [ 1 1 1 ] );
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    camproj( 'orthographic' );
  end
elseif chosenData == 6
  for p=1:numPlots
    subplot( 1, numPlots, p );
    xMinMax = [ 0 width ];
    yMinMax = [ 0 height ];
    % axis([xmin xmax ymin ymax zmin zmax cmin cmax])
    axis( [ xMinMax(1) xMinMax(2) yMinMax(1) yMinMax(2) -10000 10000 0 1 ] );
    axis on
    daspect( [ 1 1 1 ] );
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    camproj( 'orthographic' );
  end
else
  for p=1:numPlots
    subplot( 1, numPlots, p );
    xMinMax = [ 0 width ];
    yMinMax = [ 0 height ];
    % axis([xmin xmax ymin ymax zmin zmax cmin cmax])
    axis( [ xMinMax(1) xMinMax(2) yMinMax(1) yMinMax(2) -10000 10000 0 1 ] );
    axis on
    daspect( [ 1 1 1 ] );
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    camproj( 'orthographic' );
  end
end