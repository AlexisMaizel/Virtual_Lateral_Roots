function [dist] = determinePointSetSimilarity( X, Y )

[ IDX1, D1 ] = knnsearch(X, Y, 'k', 1,'distance','euclidean');
[ IDX2, D2 ] = knnsearch(Y, X, 'k', 1,'distance','euclidean');

% compute cross nearest neighbor distance
sumX = 0;
sumY = 0;
for i=1:size(D1,1)
  sumX = sumX + D1(i,1);
end
for i=1:size(D2,1)
  sumY = sumY + D2(i,1);
end
dist = ( sumX + sumY )/( size(X,1) + size(Y,1) );