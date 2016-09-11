function tile = drawColoredTile( tileIndex, min, max, res,...
  rows, columns, color )
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

% width and height of rectangle
width = tilesizeX;
height = tilesizeY;
x = tileLeftBottom(1);
y = tileLeftBottom(2);

hold on;
tile = rectangle( 'Position', [ x, y, width, height ], 'FaceColor', color, 'EdgeColor', 'none' );

