function L = drawGridArrow( startPos, endPos, tileIndex, min, max, res,...
  rows, columns, color )
topleft = [ min(1) max(2) 0 ];
topright = [ max(1) max(2) 0 ];
bottomleft = [ min(1) min(2) 0 ];
tilesizeX = distancePoints3d(topright, topleft)/res;
tilesizeY = distancePoints3d(topright, bottomleft)/res;
rowIndex = 0;
columnIndex = 0;
for r=1:rows
  for c=1:columns    
    index = c + (r-1)*columns;
    if index == tileIndex
      rowIndex = r;
      columnIndex = c;
      break;
    end
  end
end

% determine current tile bbox
cSize = (columnIndex-1)*tilesizeX;
rSize = (rowIndex-1)*tilesizeY;
tileLeftBottom = [ min(1)+cSize min(2)+rSize ];
tileRightTop = [ min(1)+cSize+tilesizeX min(2)+rSize+tilesizeY ];

% center of current tile
centerTilePos = [ (tileRightTop(1) + tileLeftBottom(1))/2. (tileRightTop(2) + tileLeftBottom(2))/2. ];

averageDirection = endPos-startPos;
averageDirection = normalize( averageDirection );

averageDirection = averageDirection.*tilesizeX/2.;
lineX = [ centerTilePos(1)-averageDirection(1), centerTilePos(1)+averageDirection(1) ];
lineY = [ centerTilePos(2)-averageDirection(2), centerTilePos(2)+averageDirection(2) ];

hold on;
averageDirection = normalizeVector3d( averageDirection );
L = quiver3( lineX(1)-averageDirection(1), lineY(1)-averageDirection(2), 50, averageDirection(1), averageDirection(2), 0,...
  tilesizeX, 'LineWidth', 0.8, 'Color', color, 'MaxHeadSize', 1.0 );
