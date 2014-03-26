function [sdLeft, sdRight] = determineSDDirections( lineDirections )

numLines = size( lineDirections, 1 );

if numLines == 0
  sdLeft = zeros( 1, 3 );
  sdRight = zeros( 1, 3 );
  return;
end

averageDirection = zeros( 1, 3 );
for l=1:numLines
  averageDirection = averageDirection + lineDirections(l, :);
end

averageDirection = averageDirection./numLines;
averageDirection = normalizeVector3d( averageDirection );

angles = zeros( numLines, 1 );
for l=1:numLines
  dir = normalizeVector3d( lineDirections(l, :) );
  angles(l, :) = acos( dot(averageDirection, dir) );
  angles(l, :) = angles(l, :)*180/pi;
end

% compute the sd of all angles
sd = std(angles);
avDir = [ averageDirection 1 ];
sdLeft = createRotationOz(-sd) * avDir';
sdRight = createRotationOz(sd) * avDir';
sdLeft = [ sdLeft(1) sdLeft(2) sdLeft(3) ];
sdRight = [ sdRight(1) sdRight(2) sdRight(3) ];
  
  