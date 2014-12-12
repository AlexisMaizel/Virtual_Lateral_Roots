function [ IDList, numEntries ] = getPrecursorIDList( precursors )
IDList = [];
% return startID if the list of precursors is empty
if strcmp( precursors, '{}' )
  IDList = -1;
  numEntries = 0;
  return;
end

% remove the first and last bracket and split the string into the
% object IDs
str = char(precursors);
str = str(2:end-1);

S = strsplit( char(str), ', ' );
numEntries = size( S, 2 );

% traverse the object list and find the earliest precursor
for o=1:numEntries
  IDList = [ IDList ; str2double(char( S( 1, o ))) ];
end
