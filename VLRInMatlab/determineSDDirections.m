function sd = determineSDDirections( lineDirections )

numLines = size( lineDirections, 1 );

if numLines == 0
  sd = 0;
  return;
end

angles = zeros( numLines, 1 );
for l=1:numLines
  slope = (lineDirections(l, 4)-lineDirections(l, 2))/(lineDirections(l, 3)-lineDirections(l, 1));
  p1 = [ 0 0 ];
  p2 = [ 1 slope ];
  dir = [ p2(1)-p1(1) p2(2)-p1(2) ];
  dir = normalize( dir );
  angles(l, :) = acos( dot([ 1 0 ], dir) );
  angles(l, :) = angles(l, :)*180/pi;
end

% compute the sd of all angles
sd = std(angles);
  