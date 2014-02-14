function nVec = getNeighborsOfObjectId( objectId, objectLinks )
% store the neighbor list for the cell with objectId
nVec = [];
dim = size( objectLinks, 1 );
for i=1:dim
  if objectLinks( i, 1 ) == objectId
    nVec( length(nVec)+1 ) = objectLinks( i, 2 );
  elseif objectLinks( i, 2 ) == objectId
    nVec( length(nVec)+1 ) = objectLinks( i, 1 );
  end
end