function result = determineCellWallIntersection( pos1, pos2, memImageStack, searchRadius, sampleStep, minVoxels )
numCellWallVoxels = 0;
startS = 0.2;
endS = 0.8;
nRadius = searchRadius;
mVoxels = minVoxels;
height = size( memImageStack, 1 );
width = size( memImageStack, 2 );
slices = size( memImageStack, 3 );
% TODO: check boundaries
if pos1(1,1)-nRadius < 0 || pos1(1,2)-nRadius < 0 || pos1(1,3)-nRadius < 0 ||...
   pos1(1,1)+nRadius > height || pos1(1,2)+nRadius > width || pos1(1,3)+nRadius > slices ||...
   pos2(1,1)-nRadius < 0 || pos2(1,2)-nRadius < 0 || pos2(1,3)-nRadius < 0 ||...
   pos2(1,1)+nRadius > height || pos2(1,2)+nRadius > width || pos2(1,3)+nRadius > slices
 nRadius = 0;
 mVoxels = 5;
end 

% create cube and starting positions
[ xs, ys, zs ] = ndgrid( -nRadius:nRadius, -nRadius:nRadius, -nRadius:nRadius );
start1 = [ pos1(1,1)-nRadius pos1(1,2)-nRadius pos1(1,3)-nRadius ];
start2 = [ pos2(1,1)-nRadius pos2(1,2)-nRadius pos2(1,3)-nRadius ];
sphereVoxels = (xs).^2 + (ys).^2 + (zs).^2 <= nRadius.^2;

% index of all non-zero elements and we only have to add the start1 and
% start2 values to get all spherical values with center pos1 and pos2
ind = find(sphereVoxels);
[ xi, yi, zi ] = ind2sub( size(sphereVoxels), ind );

si = size( xi, 1 );
prevPos = start1;
for n=1:si
  startIndex = [ xi( n, 1 )-1 yi( n, 1 )-1 zi( n, 1 )-1 ];
  p1 = startIndex + start1;
  p2 = startIndex + start2;
  for lambda = startS:sampleStep:endS
    line = p1 + lambda * ( p2 - p1 );
    curPos = int16(line);
    % only check the membrane channel if we are located at a new position
    if curPos ~= prevPos
      % NOTE that reading the TIFF, x and y are switched
      intens = memImageStack( curPos(1,2), curPos(1,1), curPos(1,3) );
      if intens > 0
        numCellWallVoxels = numCellWallVoxels + 1;
      end
    end
    prevPos = curPos;
  end
end

% subdivide by distance and number of "lines" to normalize number of voxels
%distance = distancePoints3d( pos1, pos2 );
%numCellWallVoxels
%numCellWallVoxels = numCellWallVoxels/si;
%numCellWallVoxels
if numCellWallVoxels > 1
  result = 1;
else
  result = 0;
end