function [ averageBTerm, averageTTerm ] = determineAverageTerm( contributions )
numLines = size( contributions, 1 );

averageBTerm = 0.;
averageTTerm = 0.;

for i=1:numLines
  averageBTerm = averageBTerm + contributions( i, 1 );
  averageTTerm = averageTTerm + contributions( i, 2 );
end

averageBTerm = averageBTerm / numLines;
averageTTerm = averageTTerm / numLines;