function exportTriangulation( tri, pos, timeStep, dataName, triangulationType )
fileName = strcat( '/tmp/triangulation-', dataName, '.txt' );
fileId = fopen( char(fileName), 'a' );
% first write the current time step
fprintf( fileId, '%1d\n', timeStep );

% then write points in their order of index
% before we have to transpose the matrix
if triangulationType == 1
  mat = tri.Points';
else
  mat = pos';
end

% number of points
fprintf( fileId, '%1d\n', size( mat, 2 ) );
if triangulationType == 1
  fprintf( fileId, '%4f %4f\n', mat );
else
  fprintf( fileId, '%4f %4f\n', mat( 1:2, : ) );
end

% and finally write the list of triangulations
if triangulationType == 1
  mat = tri.ConnectivityList';
else
  mat = tri';
end

% number of triangles
fprintf( fileId, '%1d\n', size( mat, 2 ) );
if size( mat, 2 ) > 0
  fprintf( fileId, '%1d %1d %1d\n', mat );
end

fprintf( fileId, '\n' );
fclose( fileId );