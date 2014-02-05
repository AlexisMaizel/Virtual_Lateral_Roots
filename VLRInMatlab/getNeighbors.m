function nVec = getNeighbors( vertexID, tri, numCells )
% store the neighbor list for each vertex ID
nVec = [];
% for all cells in the current time step
% get all the vertex IDs of its neighbors based
% on the delaunay triangulation
for v=1:numCells
  if v ~= vertexID
    E = [ vertexID v ];
    if isConnected( tri,E ) == 1
      nVec( length(nVec)+1 ) = v;
    end
  end
end