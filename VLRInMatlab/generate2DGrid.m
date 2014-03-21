function generate2DGrid( min, max, res )
topleft = [ min(1) max(2) 0 ];
topright = [ max(1) max(2) 0 ];
bottomleft = [ min(1) min(2) 0 ];
bottomright = [ max(1) min(2) 0 ];
tilesizeX = distancePoints3d(topright, topleft)/res;
tilesizeY = distancePoints3d(topright, bottomleft)/res;
for r=1:res+1
  % row lines
  indexR = (r-1)*tilesizeY;
  lineX = [ bottomleft(1) bottomright(1) ];
  lineY = [ bottomleft(2)+indexR bottomright(2)+indexR ];
  hold on;
  line( lineX, lineY, [ 50 50 ], 'Color', [ 0.9 0.9 0.9 ], 'LineWidth', 1.2 );
  % column lines
  indexC = (r-1)*tilesizeX;
  lineX = [ bottomleft(1)+indexC topleft(1)+indexC ];
  lineY = [ bottomleft(2) topleft(2) ];
  hold on;
  line( lineX, lineY, [ 50 50 ], 'Color', [ 0.9 0.9 0.9 ], 'LineWidth', 1.2 );
end