function [ xmax, ymax, zmax ] = determineCCThickness( coords, bb )
xmax = 0;
ymax = 0;
zmax = 0;
% in tiff images x and y and changed
xWidth = bb(1,5);
yWidth = bb(1,4);
zWidth = bb(1,6);
xbbStart = floor(bb(1,2));
ybbStart = floor(bb(1,1));
zbbStart = floor(bb(1,3));
yProj = zeros( xWidth, zWidth );
zProj = zeros( xWidth, yWidth );
for i=1:size(coords,1)
  yProj( coords(i, 2) - xbbStart, coords(i, 3) - zbbStart ) = 1;
  zProj( coords(i, 2) - xbbStart, coords(i, 1) - ybbStart ) = 1;
end

for x=1:xWidth
  rows = find( zProj( x, : ) == 1 );
  length = size( rows, 2 );
  if length > 0
    dist = abs( rows(1, length) - rows(1, 1) ) + 1;
    if dist > xmax
      xmax = dist;
    end
  end
end

for y=1:yWidth
  col = find( zProj( :, y ) == 1 );
  length = size( col, 1 );
  if length > 0
    dist = abs( col(length, 1) - col(1, 1) ) + 1;
    if dist > ymax
      ymax = dist;
    end
  end
end

for z=1:zWidth
  depth = find( yProj( :, z ) == 1 );
  length = size( depth, 1 );
  if length > 0
    dist = abs( depth(length, 1) - depth(1, 1) ) + 1;
    if dist > zmax
      zmax = dist;
    end
  end
end