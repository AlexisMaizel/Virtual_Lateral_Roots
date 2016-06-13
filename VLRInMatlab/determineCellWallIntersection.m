function result = determineCellWallIntersection( pos1, pos2, memImageStack, sampleStep, minVoxels )
numCellWallVoxels = 0;
startS = 0.1;
endS = 0.9;
nRadius = 5;

% TODO: check boundaries!!!

[ x1, y1, z1 ] = ndgrid( int16(pos1(1,1)-nRadius):int16(pos1(1,1)+nRadius),...
  int16(pos1(1,2)-nRadius):int16(pos1(1,2)+nRadius),...
  int16(pos1(1,3)-nRadius):int16(pos1(1,3)+nRadius) );
start1 = [ int16(pos1(1,1)-nRadius) int16(pos1(1,2)-nRadius) int16(pos1(1,3)-nRadius) ];
[ x2, y2, z2 ] = ndgrid( int16(pos2(1,1)-nRadius):int16(pos2(1,1)+nRadius),...
  int16(pos2(1,2)-nRadius):int16(pos2(1,2)+nRadius),...
  int16(pos2(1,3)-nRadius):int16(pos2(1,3)+nRadius) );
start2 = [ int16(pos2(1,1)-nRadius) int16(pos2(1,2)-nRadius) int16(pos2(1,3)-nRadius) ];
sphereVoxels1 = ( int16(x1) - int16(pos1(1,1)) ).^2 + ( int16(y1) - int16(pos1(1,2)) ).^2 + ( int16(z1) - int16(pos1(1,3)) ).^2 <= nRadius.^2;
sphereVoxels2 = ( int16(x2) - int16(pos2(1,1)) ).^2 + ( int16(y2) - int16(pos2(1,2)) ).^2 + ( int16(z2) - int16(pos2(1,3)) ).^2 <= nRadius.^2;

% index of all non-zero elements
ind1 = find(sphereVoxels1);
[ i11, i21, i31 ] = ind2sub( size(sphereVoxels1), ind1 );
ind2 = find(sphereVoxels2);
[ i12, i22, i32 ] = ind2sub( size(sphereVoxels2), ind2 );

si = size( i12, 1 );
for n=1:si
  p1 = [ int16(i11( n, 1 )) int16(i21( n, 1 )) int16(i31( n, 1 )) ] + start1;
  p2 = [ int16(i12( n, 1 )) int16(i22( n, 1 )) int16(i32( n, 1 )) ] + start2;
  for lambda = startS:sampleStep:endS
    %line = p(n,:) + lambda * ( q(n,:) - p(n,:) );
    line = p1 + lambda * ( p2 - p1 );
    curPos = [ int16(line(1,1)) int16(line(1,2)) int16(line(1,3)) ];
    % NOTE that reading the TIFF, x and y are switched
    intens = memImageStack( curPos(1,2), curPos(1,1), curPos(1,3) );
    if intens > 0
      numCellWallVoxels = numCellWallVoxels + 1;
    end
  end
end

%numCellWallVoxels
if numCellWallVoxels > minVoxels
  result = 1;
else
  result = 0;
end