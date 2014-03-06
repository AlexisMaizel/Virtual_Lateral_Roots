function pos = getCellPosition( objectId, tri, cellIds, triangulationType, matPos )
dim = size( cellIds, 1 );
for c=1:dim
  if objectId == cellIds( c, 1 )
    if triangulationType == 1
      pos = [ tri.Points( c, 1 ) tri.Points( c, 2 ) tri.Points( c, 3 ) ];
    else
      pos = [ matPos( c, 1 ) matPos( c, 2 ) matPos( c, 3 ) ];
    end
    break;
  end
end