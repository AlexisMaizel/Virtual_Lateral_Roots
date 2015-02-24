function [ averageSlope, averageDirection ] = determineAverageSlope( lineDirections )
numLines = size( lineDirections, 1 );
numPositive = 0;
averagePSlope = 0;
numNegative = 0;
averageNSlope = 0;
averageStart = zeros(1,2);
averageEnd = zeros(1,2);
for l=1:numLines
  % first compute the slope
  startPos = lineDirections(l, 1:2);
  endPos = lineDirections(l, 3:4);
  slope = (endPos(2)-startPos(2))/(endPos(1)-startPos(1));
  
  p1 = [ 0 0 ];
  p2 = [ 1 slope ];
  dir = [ p2(1)-p1(1) p2(2)-p1(2) ];
  averageStart = averageStart + startPos;
  averageEnd = averageEnd + endPos;
  dir = normalize( dir );
  angle = acos( dot([ 1 0 ], dir) );
  angle = angle*180/pi;
  
  if slope >= 0
    numPositive = numPositive + 1;
    averagePSlope = averagePSlope + angle;
  end
  
  if slope < 0
    numNegative = numNegative + 1;
    averageNSlope = averageNSlope + angle;
  end
end

% compute the average of the magnitude
averageDirection = averageEnd./numLines - averageStart./numLines;

if numPositive > 0
  averagePSlope = averagePSlope/numPositive;
  if numNegative == 0
    rad = degtorad(averagePSlope);
    rotMat = createRotationOz(rad);
    dir = [ 1 0 0 1 ];
    dir = rotMat * dir';
    averageSlope = dir(2)/dir(1);
    return;
  end
end
if numNegative > 0
  averageNSlope = averageNSlope/numNegative;
  if numPositive == 0
    rad = degtorad(-averageNSlope);
    rotMat = createRotationOz(rad);
    dir = [ -1 0 0 1 ];
    dir = rotMat * dir';
    averageSlope = dir(2)/dir(1);
    return;
  end
end

rad = degtorad(averagePSlope);
rotMat = createRotationOz(rad);
dir1 = [ 1 0 0 1 ];
dir1 = rotMat * dir1';

rad = degtorad(-averageNSlope);
rotMat = createRotationOz(rad);
dir2 = [ 1 0 0 1 ];
dir2 = rotMat * dir2';

angle = acos( dot([dir1(1) dir1(2)], [dir2(1) dir2(2)]) );
angle = angle*180/pi;

% if the angle is smaller than 90 degrees then the average should be in
% this area else it is located in the mid of the other angle
if angle < 90
  rad = degtorad(-angle/2.);
  rotMat = createRotationOz(rad);
  dir = [ dir1(1) dir1(2) 0 1 ];
  dir = rotMat * dir';
  averageSlope = dir(2)/dir(1);
else
  angle = 180 - angle;
  rad = degtorad(-angle/2.);
  rotMat = createRotationOz(rad);
  dir = [ dir2(1) dir2(2) 0 1 ];
  dir = rotMat * dir';
  averageSlope = dir(2)/dir(1);
end
