function objectID = getPrecursorID( startID, precursors, previousObjectIdList )

% return startID if the list of precursors is empty
if strcmp( precursors, '{}' )
  objectID = startID;
  return;
end

% remove the first and last bracket and split the string into the
% object IDs
S = strsplit( precursors(2:end-1), ', ' );
dim = size( S, 2 );

% traverse the object list and find the earliest precursor
while dim > 0
  for o=1:size( previousObjectIdList, 2 )
    if S( 1, dim ) == previousObjectIdList( 1, o )
      objectID = S( 1, dim );
      return;
    end
  end
  dim = dim-1;
end

% never should be here
objectID = startID;
