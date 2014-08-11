function exportTriangulationWithNewPos( tri, timeStep, nextPos, dataName )
fileName = strcat( '../FinalVLRForMatlab/triangulation-', dataName, '.txt' );
fileId = fopen( char(fileName), 'a' );
% first write the current time step
fprintf( fileId, '%1d\n', timeStep );

nextPos = nextPos';
% number of new points
fprintf( fileId, '%1d\n', size( nextPos, 2 ) );
fprintf( fileId, '%4f %4f\n', nextPos( 1:2, : ) );

% and finally write the list of triangulations
mat = tri.ConnectivityList';
% number of triangles
fprintf( fileId, '%1d\n', size( mat, 2 ) );
if size( mat, 2 ) > 0
  fprintf( fileId, '%1d %1d %1d\n', mat );
end

fprintf( fileId, '\n' );
fclose( fileId );