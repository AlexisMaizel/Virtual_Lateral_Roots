function exportResult( totalMinAxes, totalMaxAxes, bezierOffset, numData, pureDataStr,...
  allCurT, allCurCells, allCells, curI, numNormCells,...
  viewStr, imageDir, f )
hold off;
set( f,'nextplot','replacechildren' );
viewOffset = 100;
% axis([xmin xmax ymin ymax zmin zmax cmin cmax])
axis( [ totalMinAxes(1)-bezierOffset totalMaxAxes(1)+bezierOffset...
  totalMinAxes(2)-bezierOffset totalMaxAxes(2)+bezierOffset...
  totalMinAxes(3)-viewOffset totalMaxAxes(3)+viewOffset 0 1 ] );
axis on
daspect( [ 1 1 1 ] );

% legend
stringData = [];
for d=1:numData
  stringData = [ stringData ; strcat( pureDataStr( 1, d ),...
    ' - T', num2str(allCurT(d, 1)),...
    ' - C', num2str(allCurCells(d, 1)), '/', num2str(allCells(d, 1)) ) ];
end
%stringData = [ stringData ; strcat( 'Average - C', num2str( numAverageCells ) ) ];
str = cellstr( stringData );
%leg = legend( '120830', '121204', '121211', '130508', '130607', 'Average' );
leg = legend( str );
set(leg, 'Location', 'NorthWestOutside');
linecolors = { [ 1 0 1 ], [ 1 1 0 ], [ 1 0 0 ], [ 0 1 0 ], [ 0 0 1 ], [ 0 1 1 ] };
legendlinestyles( leg, {}, {}, linecolors );

grid off;
xlabel('X');
ylabel('Y');
zlabel('Z');
title( strcat( 'Normalized Step ', num2str(curI), '-Cells', num2str(numNormCells) ) );

if curI < 10
  digit = strcat( viewStr, '_00' );
elseif curI < 100
  digit = strcat( viewStr, '_0' );
else
  digit = strcat( viewStr, '_' );
end

% output with number of cells
filePath = strcat( imageDir, digit, num2str(curI), '-Cells', num2str(numNormCells), '.png' );
saveas( gcf, char(filePath) );