function index = getTileIndex( pos, min, max, res, rows, columns )
topleft = [ min(1) max(2) 0 ];
topright = [ max(1) max(2) 0 ];
bottomleft = [ min(1) min(2) 0 ];
tilesizeX = distancePoints3d(topright, topleft)/res;
% for rectangle tiles:
%tilesizeY = distancePoints3d(topright, bottomright)/res;
% for square tiles:
tilesizeY = distancePoints3d(topright, bottomleft)/res;
index = 0;
for r=1:rows
  for c=1:columns
    % determine current tile bbox
    cSize = (c-1)*tilesizeX;
    rSize = (r-1)*tilesizeY;
    tileLeftBottom = [ min(1)+cSize min(2)+rSize ];
    tileRightTop = [ min(1)+cSize+tilesizeX min(2)+rSize+tilesizeY ];
    
    % determine index in which pos is located
    if pos(1) >= tileLeftBottom(1) &&...
        pos(1) < tileRightTop(1) &&...
        pos(2) >= tileLeftBottom(2) &&...
        pos(2) < tileRightTop(2)
      index = c + (r-1)*columns;
      return;
    end
  end
end

% default value whihch
disp('Point was not found in grid!')
index = 1;