function objectID = getPrecursorID( startID, precursors, previousObjectIdList )

% return startID if the list of precursors is empty
if strcmp( precursors, '{}' )
  objectID = startID;
  return;
end

% remove the first and last bracket and split the string into the
% object IDs
str = char(precursors);
str = str(2:end-1);

S = strsplit( char(str), ', ' );
dim = size( S, 2 );

% traverse the object list and find the earliest precursor
while dim > 0
  for o=1:size( previousObjectIdList, 1 )
    if str2double(char( S( 1, dim ))) == previousObjectIdList( o, 1 )
      objectID = previousObjectIdList( o, 1 );
      return;
    end
  end
  dim = dim-1;
end

% never should be here
objectID = startID;
