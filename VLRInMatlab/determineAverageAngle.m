function averageAngle = determineAverageAngle( lineDirections )

numLines = size( lineDirections, 1 );

if numLines == 0
  averageAngle = -1;
  return;
end

averageAngle = 0;
for l=1:numLines
  % first compute the slope
  startPos = lineDirections(l, 1:2);
  endPos = lineDirections(l, 3:4);
  slope = (endPos(2)-startPos(2))/(endPos(1)-startPos(1));
  
  % determine direction from left to right and compute the smallest angle
  dir = endPos - startPos;
  normLine = norm(dir);
  angle = acos( (dot(dir, [ 1 0 ]))/normLine );
  angle = angle*180/pi;
  
  % if the slope is negative then subtract the angle from 180 degrees
%   if slope < 0
%     if angle >= 45
%       angle = 180 - angle;
%       %angle = angle - 180;
%     else
%       angle = -angle;
%       %angle = angle - 180;
%       %angle = 180 - angle;
%     end
%   end

  if slope < 0
      angle = 180 - angle;
  end
  
  if slope == 0
    zero
  end
  
  averageAngle = averageAngle + angle;
end

averageAngle = averageAngle/numLines;