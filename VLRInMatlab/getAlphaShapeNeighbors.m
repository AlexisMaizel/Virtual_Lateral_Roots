function nVec = getAlphaShapeNeighbors( vertexID, tri )
% store the neighbor list for each vertex ID
nVec = [];
% for all cells in the current time step
% get all the vertex IDs of its neighbors based
% on the triangulation of the alpha shape
dim = size( tri, 1 );
for v=1:dim
  elements = [ 0 0 0 0 ];
  vertexIdIncluded = 0;
  for l=1:4
    elements( 1, l ) = tri( v, l );
    if vertexID == tri( v, l )
      vertexIdIncluded = 1;
    end
  end
  
  if vertexIdIncluded == 1
    for l=1:4
      if vertexID ~= elements( 1, l )
        nVec( length(nVec)+1 ) = elements( 1, l );
      end
    end
  end
end

% at last get only the unique vector of nVec
nVec = unique( nVec );