function conPoints = generateAutomaticContourPoints( curPos, bDist, eps, first, dataName )
% I always choose 16 contour points for which 7 are used
% at the top and bottom while 1 is used for left and right
numContourMarks = 16;

% initialize 16 contour points as the boundary of all
% real data points
offset = 200.;
minX = min( curPos(:,1) );
minY = min( curPos(:,2) );
maxX = max( curPos(:,1) );
maxY = max( curPos(:,2) );
minX = minX - offset;
maxX = maxX + offset;
minY = minY - offset;
maxY = maxY + offset;
cPoints = zeros( numContourMarks, 3 );
numRow = 7;
numColumn = 1;
xStep = (maxX-minX)/(numRow-1);
% top
for c=1:numRow
  cPoints(c, 1:2) = [ minX+(c-1)*xStep maxY ];
end
% bottom
for c=1:numRow
  cPoints(c+numRow+numColumn, 1:2) = [ maxX-(c-1)*xStep minY ];
end
% left
cPoints(2.*(numRow+numColumn), 1:2) = [ minX minY+(maxY-minY)/2. ];
% right
cPoints(numRow+numColumn, 1:2) = [ maxX minY+(maxY-minY)/2. ];

% find nearest real data point for each contour point
for c=1:size( cPoints, 1 )
  dist = 100000;
  index = 1;
  for p=1:size( curPos, 1 )
    curDist = distancePoints3d( cPoints(c,:), curPos(p,:) );
    if( curDist < dist )
      dist = curDist;
      index = p;
    end
  end
  % direction vector
  dir = curPos(index,:)-cPoints(c,:);
  % translate contour point such that the new distance is smaller
  % than the desired bDist
  factor = 0.5;
  step = 1;
  newPos = zeros(1,3);
  while dist > bDist 
    newPos = cPoints(c,:) + (1.-factor)*dir;
    dist = distancePoints3d( newPos, curPos(index,:) );
    factor = factor/2.;
    step = step + 1;
  end
  cPoints(c,:) = newPos;
end

% set the subdivision level for the initial 16 contour points
% which means that between each pair of points an additional
% contour point is inserted in the middle
% 1 -> original 16 points
% 2 -> 2*16 contour points
level = 1;
if level ~= 1
  conPoints = zeros( level*numContourMarks, 3 );
  for p=1:size(cPoints,1)
    startPos = cPoints( p, : );
    if p == size(cPoints,1)
      endPos = cPoints( 1, : );
    else
      endPos = cPoints( p+1, : );
    end
    for l=1:level
      factor = 1./l;
      newPos = factor * startPos + (1-factor) * endPos;
      conPoints( (p-1)*level+l, : ) = newPos;
    end
  end
else
  conPoints = cPoints;
end

% export of initial contour positions used in the model before the
% offset is applied such that the model contour points are always
% within the increased surface realized with the eps offset
if first == true
  fileName = strcat( '/tmp/conPoints-', dataName, '_auto.txt' );
  fileId = fopen( char(fileName), 'w' );
  fprintf( fileId, '%1d\n', size(conPoints,1) );
  for p=1:size(conPoints,1)
    fprintf( fileId, '%4f %4f\n', conPoints(p,1:2) );
  end
  fprintf( fileId, '\n' );
  fclose( fileId );
end

% apply the offset to the points according to their position
% on the boundary of the surface
% top left point
cPoints(1,1:2) = [ cPoints(1,1)-eps cPoints(1,2)+eps ];
for p=2:6
  cPoints(p,1:2) = [ cPoints(p,1) cPoints(p,2)+eps ];
end
% top right point
cPoints(7,1:2) = [ cPoints(7,1)+eps cPoints(7,2)+eps ];
cPoints(8,1:2) = [ cPoints(8,1)+eps cPoints(8,2) ];
% bottom right point
cPoints(9,1:2) = [ cPoints(9,1)+eps cPoints(9,2)-eps ];
for p=10:14
  cPoints(p,1:2) = [ cPoints(p,1) cPoints(p,2)-eps ];
end
% bottom left point
cPoints(15,1:2) = [ cPoints(15,1)-eps cPoints(15,2)-eps ];
cPoints(16,1:2) = [ cPoints(16,1)-eps cPoints(16,2) ];

% set the subdivision level again for the slightly increased surface
if level ~= 1
  for p=1:size(cPoints,1)
    startPos = cPoints( p, : );
    if p == size(cPoints,1)
      endPos = cPoints( 1, : );
    else
      endPos = cPoints( p+1, : );
    end
    for l=1:level
      factor = 1./l;
      newPos = factor * startPos + (1-factor) * endPos;
      conPoints( (p-1)*level+l, : ) = newPos;
    end
  end
else
  conPoints = cPoints;
end
