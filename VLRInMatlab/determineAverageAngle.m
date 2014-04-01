function averageAngle = determineAverageDirection( lineDirections )

numLines = size( lineDirections, 1 );

% if numLines == 0
%   averageDirection = zeros( 1, 4 );
%   return;
% end

averageAngle = 0;
for l=1:numLines
  % first compute the slope
  startPos = lineDirections(l, 1:2);
  endPos = lineDirections(l, 3:4);
  slope = (endPos(2)-startPos(2))/(endPos(1)-startPos(1));
  negSlope = 0;
  if slope < 0
    negSlope = 1;
  end
  
  % check both direction vectors
  firstDir = endPos - startPos;
  secondDir = startPos - endPos;
  normLine = norm(firstDir);
  
  % and choose this one that has the smallest angle compared with the
  % x-axis
  angle1 = acos( (dot(firstDir, [ 1 0 ]))/normLine );
  angle1 = angle1*180/pi;
  angle2 = acos( (dot(secondDir, [ 1 0 ]))/normLine );
  angle2 = angle2*180/pi;
  if angle1 < angle2
    angle = angle1;
  else
    angle = angle2;
  end
  
  % if the slope is negative then subtract the angle from 180 degrees
  if negSlope == 1
    angle = 180 - angle;
  end

  averageAngle = averageAngle + angle;
end

averageAngle = averageAngle/numLines;