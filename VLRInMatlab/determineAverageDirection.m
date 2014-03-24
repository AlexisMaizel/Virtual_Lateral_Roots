function averageDirection = determineAverageDirection( lineDirections )

numLines = size( lineDirections, 1 );

if numLines == 0
  averageDirection = [ 0 0 0 ];
  return;
end

averageDirection = zeros( 1, 3 );
for l=1:numLines
  averageDirection = averageDirection + lineDirections(l, :);
end

averageDirection = averageDirection./numLines;

averageDirection = normalizeVector3d( averageDirection );