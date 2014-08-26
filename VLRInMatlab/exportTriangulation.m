function exportTriangulation( tri, timeStep, dataName, triangulationType )
fileName = strcat( '/tmp/triangulation-', dataName, '.txt' );
fileId = fopen( char(fileName), 'a' );
% first write the current time step
fprintf( fileId, '%1d\n', timeStep );

% then write points in their order of index
% before we have to transpose the matrix
if triangulationType == 1
  mat = tri.Points';
else
  mat = tri;
end

% number of points
fprintf( fileId, '%1d\n', size( mat, 2 ) );
fprintf( fileId, '%4f %4f\n', mat );

% and finally write the list of triangulations
mat = tri.ConnectivityList';
% number of triangles
fprintf( fileId, '%1d\n', size( mat, 2 ) );
if size( mat, 2 ) > 0
  fprintf( fileId, '%1d %1d %1d\n', mat );
end

fprintf( fileId, '\n' );
fclose( fileId );