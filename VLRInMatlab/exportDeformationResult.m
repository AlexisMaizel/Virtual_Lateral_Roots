function exportDeformationResult( totalMinAxes, totalMaxAxes, bezierOffset, numData, pureDataStr,...
  curTC, curTN, curCellsC, curCellsN, allCellsC, allCellsN,...
  curI, deltaI, viewStr, imageDir, f, numPlots, renderTermType,...
  contributions, positiveEigenvalueVector, startData, endData )
hold off;
set( f,'nextplot','replacechildren' );
fontSize = 7;
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
    set(leg, 'Location', 'NorthEast');
    set(leg, 'FontSize', fontSize );
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
    set(leg, 'Location', 'NorthEast');
    set(leg, 'FontSize', fontSize );
    linecolors = { [ 1 0 1 ], [ 1 1 0 ], [ 1 0 0 ], [ 0 1 0 ], [ 0 0 1 ], [ 0 1 1 ] };
    legendlinestyles( leg, {}, {}, linecolors );
    
    grid off;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title( strcat( {'Normalized step '}, num2str(curI+deltaI) ) );
  % third plot
  elseif pl == 3
    grid off;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title( strcat( {'Average deformation between normalized step '}, num2str(curI),...
      {' and '} , num2str(curI+deltaI) ) );
  % forth plot
  else
    BTerm = 0.;
    TTerm = 0.;
    PosEV = 0;
    NegEV = 0;
    dimTotal = 0;
    for d=startData:endData
      dim = size( contributions{d}, 1 );
      for t=1:dim
        BTerm = BTerm + contributions{d}( t, 1 );
        TTerm = TTerm + contributions{d}( t, 2 );
      end
      dimTotal = dimTotal + dim;
      
      dimP = size( positiveEigenvalueVector{d}, 1 );
      for c=1:dimP
        for p=1:size( positiveEigenvalueVector{d}, 2 )
          % positive ev
          if positiveEigenvalueVector{d}( c, p ) == 1
            PosEV = PosEV + 1;
          else
            NegEV = NegEV + 1;
          end
        end
      end
    end
    
    if dimTotal > 0
      BTerm = BTerm/dimTotal;
      TTerm = TTerm/dimTotal;
    end
    
    stringData = [ strcat( {'B Term: '}, num2str(BTerm) ) ;...
      strcat( {'T Term: '}, num2str(TTerm) ) ;...
      strcat( {'Positive eigenvalues: '}, num2str(PosEV) ) ;...
      strcat( {'Negative eigenvalues: '}, num2str(NegEV) ) ];
    
    str = cellstr( stringData );
    leg = legend( str );
    set(leg, 'Location', 'NorthEast');
    set(leg, 'FontSize', fontSize );
    linecolors = { [ 1 1 1 ], [ 1 1 1 ], [ 0 0 1 ], [ 1 0 0 ] };
    legendlinestyles( leg, {}, {}, linecolors );
      
    grid off;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title( strcat( {'Real Deformations between normalized step '}, num2str(curI),...
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
filePath = strcat( imageDir, digit, num2str(curI) );
export_fig( gcf, char(filePath), '-m2', '-png' );
%saveas( gcf, char(filePath) );