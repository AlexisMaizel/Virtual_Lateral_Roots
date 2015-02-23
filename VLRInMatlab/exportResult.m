function exportResult( totalMinAxes, totalMaxAxes, bezierOffset, numData, pureDataStr,...
  curTC, curTN, curCellsC, curCellsN, allCellsC, allCellsN,...
  curI, deltaI, viewStr, imageDir, f, numPlots, renderTermType )
hold off;
set( f,'nextplot','replacechildren' );
for pl=1:numPlots
  subplot( numPlots, 1, pl );
  viewOffset = 100;
  % axis([xmin xmax ymin ymax zmin zmax cmin cmax])
  axis( [ totalMinAxes(1)-bezierOffset totalMaxAxes(1)+bezierOffset...
    totalMinAxes(2)-bezierOffset totalMaxAxes(2)+bezierOffset...
    totalMinAxes(3)-viewOffset totalMaxAxes(3)+viewOffset 0 1 ] );
  axis on
  daspect( [ 1 1 1 ] );
  
  % legend
  stringData = [];
  % first plot
  if pl == 1
    for d=1:numData
      stringData = [ stringData ; strcat( pureDataStr( 1, d ),...
        ' - T', num2str(curTC(d)),...
        ' - C', num2str(curCellsC(d)), '/', num2str(allCellsC(d)) ) ];
    end
    
    str = cellstr( stringData );
    leg = legend( str );
    set(leg, 'Location', 'NorthWestOutside');
    linecolors = { [ 1 0 1 ], [ 1 1 0 ], [ 1 0 0 ], [ 0 1 0 ], [ 0 0 1 ], [ 0 1 1 ] };
    legendlinestyles( leg, {}, {}, linecolors );
    
    grid off;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title( strcat( {'Normalized step '}, num2str(curI) ) );
  % second plot
  elseif pl == 2
    for d=1:numData
      stringData = [ stringData ; strcat( pureDataStr( 1, d ),...
        ' - T', num2str(curTN(d)),...
        ' - C', num2str(curCellsN(d)), '/', num2str(allCellsN(d)) ) ];
    end
    
    str = cellstr( stringData );
    leg = legend( str );
    set(leg, 'Location', 'NorthWestOutside');
    linecolors = { [ 1 0 1 ], [ 1 1 0 ], [ 1 0 0 ], [ 0 1 0 ], [ 0 0 1 ], [ 0 1 1 ] };
    legendlinestyles( leg, {}, {}, linecolors );
    
    grid off;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title( strcat( {'Normalized step '}, num2str(curI+deltaI) ) );
  % third plot
  else
    stringData = 'Deformation orientation';
    
    str = cellstr( stringData );
    leg = legend( str );
    set(leg, 'Location', 'NorthWestOutside');
    linecolors = { [ 0 0 0 ] };
    legendlinestyles( leg, {}, {}, linecolors );
    
    grid off;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title( strcat( {'Deformation between normalized step '}, num2str(curI),...
      {' and '} , num2str(curI+deltaI) ) );
  end
end

if curI < 10
  digit = strcat( viewStr, '_00' );
elseif curI < 100
  digit = strcat( viewStr, '_0' );
else
  digit = strcat( viewStr, '_' );
end

% B term only
if renderTermType == 1
  imageDir = strcat( imageDir, 'BTerm_' );
% T term only
elseif renderTermType == 2
  imageDir = strcat( imageDir, 'TTerm_' );
% Both terms
else
  imageDir = strcat( imageDir, 'AllTerms_' );
end

% output with number of cells
filePath = strcat( imageDir, digit, num2str(curI), '.png' );
saveas( gcf, char(filePath) );