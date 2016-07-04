function [S, cellCenters, maxIntensities] = mergeNearestCCs( imageStack, cellRadius, randomIntensities )
% determine connected components in image
% can be 6, 18 or 26
connectivity = 26;
CC = bwconncomp( imageStack( :, :, : ), connectivity );
% determine properties of ccs
S = regionprops( CC, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );
% merge nearest cc
if size(S, 1) > 1
  X = zeros( size(S, 1), 3 );
  for i=1:size(S, 1)
    X( i, : ) = S(i,:).Centroid;
  end
  %Y = pdist(X);
  Z = linkage( X, 'centroid' );
  T = cluster( Z, 'cutoff', cellRadius, 'criterion', 'distance' );
  numC = max( T );
  cellCenters = zeros( numC, 3 );
  maxIntensities = zeros( numC, 1 );
  for i=1:size(T, 1)
    if randomIntensities == 1
      a = 800;
      b = 1300;
      r = (b-a).*rand(1,1) + a;
      maxInt = int64(r);
    else
      maxInt = max( imageStack( S(i, :).PixelIdxList ) );
    end
    cellCenters( T(i,1), : ) = cellCenters( T(i,1), : ) + S(i, :).Centroid;
    maxIntensities( T(i,1), 1 ) = maxIntensities( T(i,1), 1 ) + maxInt;
  end
  for c=1:numC
    cellCenters( c, : ) = cellCenters( c, : ) / nnz(T==c);
    maxIntensities( c, 1 ) = maxIntensities( c, 1 ) / nnz(T==c);
  end
end