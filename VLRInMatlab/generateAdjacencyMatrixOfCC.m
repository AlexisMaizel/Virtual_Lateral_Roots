function ccAdjacencyMap = generateAdjacencyMatrixOfCC( CC, S, connectivity )
ccAdjacencyMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );

% init adjacency matrix for each cc with zeros
for c=1:CC.NumObjects
  list = S(c,:).PixelList;
  numPixels = size( list, 1 );
  bb = S(c,:).BoundingBox;
  
  xmin = min( list(:, 1) );
  ymin = min( list(:, 2) );
  zmin = min( list(:, 3) );
  
  mat = zeros( bb(1, 4), bb(1, 5), bb(1, 6) );
  for p=1:numPixels
    pos = list( p, : );
    pos = pos - [ xmin-1 ymin-1 zmin-1 ];
    mat( pos(1,1), pos(1,2), pos(1,3) ) = p;
  end
  
  % create adjacency matrix
  adjMat = zeros( numPixels, numPixels );
  for p=1:numPixels
    pos = list( p, : );
    pos = pos - [ xmin-1 ymin-1 zmin-1 ];
    [ connList ] = checkConnectivity( pos, mat, [ bb(1, 4) bb(1, 5) bb(1, 6) ], connectivity );
    for l=1:size( connList, 2 )
      if p ~= connList(l)
        adjMat( p, connList(l) ) = 1;
        adjMat( connList(l), p ) = 1;
      end
    end
  end
  ccAdjacencyMap(c) = adjMat;
end
