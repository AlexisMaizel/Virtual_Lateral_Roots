function objectIdLinks = generateLinksBetweenObjectsId( edgeList, objectIdVector )

  dim = size( edgeList( :, 1 ) );
  objectIdLinks = [ dim:2 ];
  
  for i=1:dim
    objectIdLinks( i, 1 ) = objectIdVector( edgeList( i, 1 ), 1 );
    objectIdLinks( i, 2 ) = objectIdVector( edgeList( i, 2 ), 1 );
  end