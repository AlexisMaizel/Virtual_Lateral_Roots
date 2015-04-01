function [ averageDirection ] = determineAverageDirection( lineDirections )
numLines = size( lineDirections, 1 );
averageStart = zeros(1,2);
averageEnd = zeros(1,2);
for l=1:numLines
  averageStart = averageStart + lineDirections(l, 1:2);
  averageEnd = averageEnd + lineDirections(l, 3:4);
end

% compute the average of the magnitude
averageDirection = averageEnd./numLines - averageStart./numLines;