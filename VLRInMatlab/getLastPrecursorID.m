function objectID = getLastPrecursorID( precursors )
% return startID if the list of precursors is empty
if strcmp( precursors, '{}' )
  objectID = -1;
  return;
end

% remove the first and last bracket and split the string into the
% object IDs
str = char(precursors);
str = str(2:end-1);

S = strsplit( char(str), ', ' );
dim = size( S, 2 );

% store the last precursor and return the object id
objectID = str2double(char( S( 1, dim )));
