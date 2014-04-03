function SD = drawStandardDeviationArea( sd, averageSlope, tileIndex, min, max, res,...
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
tileRightTop = [ min(1)+cSize+tilesizeX min(2)+rSize+tilesizeY ];

% center of current tile
centerTilePos = [ (tileRightTop(1) + tileLeftBottom(1))/2. (tileRightTop(2) + tileLeftBottom(2))/2. ];

% create line with desired slope
p1 = [ 0 0 ];
p2 = [ 1 averageSlope ];
averageDirection = [ p2(1)-p1(1) p2(2)-p1(2) ];
averageDirection = normalize( averageDirection );
averageDirection = [ averageDirection(1) averageDirection(2) 0 1 ];

% perform rotation depending on slope
rad = degtorad(sd);
rotMat = createRotationOz(rad);
sd1 = rotMat * averageDirection';
sd1 = sd1.*tilesizeX/2.;
rad = degtorad(-sd);
rotMat = createRotationOz(rad);
sd2 = rotMat * averageDirection';
sd2 = sd2.*tilesizeX/2.;

X = [ centerTilePos(1)+sd1(1) centerTilePos(1)+sd2(1) centerTilePos(1)...
  centerTilePos(1)-sd1(1) centerTilePos(1)-sd2(1) centerTilePos(1) ];
Y = [ centerTilePos(2)+sd1(2) centerTilePos(2)+sd2(2) centerTilePos(2)...
  centerTilePos(2)-sd1(2) centerTilePos(2)-sd2(2) centerTilePos(2) ];

hold on;
SD = fill( X, Y, color );
