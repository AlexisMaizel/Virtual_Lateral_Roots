function rec = drawContributionRectangle( tileIndex, min, max, res,...
  rows, columns, color, termType )
topleft = [ min(1) max(2) 0 ];
topright = [ max(1) max(2) 0 ];
bottomleft = [ min(1) min(2) 0 ];
tilesizeX = distancePoints3d(topright, topleft)/res;
% for rectangle tiles:
%tilesizeY = distancePoints3d(topright, bottomright)/res;
% for square tiles:
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

% width and height of rectangle
width = tilesizeX/2.;
height = tilesizeY/2.;

% top left rectangle
if termType == 0
  x = tileLeftBottom(1);
  y = tileLeftBottom(2) + height;
% top right rectangle
elseif termType == 1
  x = tileLeftBottom(1) + width;
  y = tileLeftBottom(2) + height;
% bottom left rectangle
elseif termType == 2
  x = tileLeftBottom(1);
  y = tileLeftBottom(2);
% bottom right rectangle
elseif termType == 3
  x = tileLeftBottom(1) + width;
  y = tileLeftBottom(2);
else
  % do nothing
end

hold on;
%averageDirection = normalizeVector3d( averageDirection );
rec = rectangle( 'Position', [ x, y, width, height ], 'FaceColor', color, 'EdgeColor', 'none' );

