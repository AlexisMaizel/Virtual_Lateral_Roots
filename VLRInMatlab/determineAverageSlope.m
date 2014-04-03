function averageSlope = determineAverageSlope( lineDirections )
numLines = size( lineDirections, 1 );
averageSlope = 0;
slopeLimit = 2;
numPositive = 0;
averagePSlope = 0;
numNegative = 0;
averageNSlope = 0;
for l=1:numLines
  % first compute the slope
  startPos = lineDirections(l, 1:2);
  endPos = lineDirections(l, 3:4);
  slope = (endPos(2)-startPos(2))/(endPos(1)-startPos(1));
  
  if slope > slopeLimit && numLines > 1
    slope = slopeLimit;
  elseif slope < -slopeLimit && numLines > 1
    slope = -slopeLimit;
  end
  
  if slope >= 0
    numPositive = numPositive + 1;
    averagePSlope = averagePSlope + slope;
  end
  
  if slope < 0
    numNegative = numNegative + 1;
    averageNSlope = averageNSlope + slope;
  end
  
  %averageSlope = averageSlope + slope;
end

%averageSlope = averageSlope/numLines;

if numPositive > 0
  averagePSlope = averagePSlope/numPositive;
  if numNegative == 0
    averageSlope = averagePSlope;
    return;
  end
end
if numNegative > 0
  averageNSlope = averageNSlope/numNegative;
  if numPositive == 0
    averageSlope = averageNSlope;
    return;
  end
end

% determine smallest angle between positive and negative slope sets
% create line with desired slope
p1 = [ 0 0 ];
q1 = [ 1 averagePSlope ];
dir1 = [ q1(1)-p1(1) q1(2)-p1(2) ];
dir1 = normalize( dir1 );
p2 = [ 0 0 ];
q2 = [ 1 averageNSlope ];
dir2 = [ q2(1)-p2(1) q2(2)-p2(2) ];
dir2 = normalize( dir2 );

angle = acos( dot(dir1, dir2) );
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
  rad = degtorad(-angle/2.);
  rotMat = createRotationOz(rad);
  dir = [ dir2(1) dir2(2) 0 1 ];
  dir = rotMat * dir';
  averageSlope = dir(2)/dir(1);
end
