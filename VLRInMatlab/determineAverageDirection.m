function [ averageDirection ] = determineAverageDirection( lineDirections )
numLines = size( lineDirections, 1 );
averageStart = zeros(1,2);
averageEnd = zeros(1,2);
for l=1:numLines
  startPos = lineDirections(l, 1:2);
  endPos = lineDirections(l, 3:4);
  averageStart = averageStart + startPos;
  averageEnd = averageEnd + endPos;
end

% compute the average of the magnitude
averageDirection = averageEnd./numLines - averageStart./numLines;