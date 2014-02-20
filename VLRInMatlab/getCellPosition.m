function pos = getCellPosition( objectId, tri, cellIds )
dim = size( cellIds, 1 );
for c=1:dim
  if objectId == cellIds( c, 1 )
    pos = [ tri.Points( c, 1 ) tri.Points( c, 2 ) tri.Points( c, 3 ) ];
    break;
  end
end