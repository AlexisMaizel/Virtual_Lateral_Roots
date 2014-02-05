function lineageId = getLineageId( objectId, LEntry )

% if no precursor exists then the object id is the lineage id
if strcmp( LEntry, '{}' ) == 1
  lineageId = objectId;
  % else the lineage id is given by the first numerical entry of the
  % precursor string
else
  S = char( LEntry );
  C = strsplit( S, {',','{','}',' '} );
  lineageId = str2double( C(2) );
end