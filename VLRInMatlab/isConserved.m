function exists = isConserved( objectId, objectLinks )
dim = size( objectLinks, 1 );
exists = 0;
for i=1:dim
  if ( objectId == objectLinks( i, 1 ) ) || ( objectId == objectLinks( i, 2 ) )
    exists = 1;
    break;
  end
end