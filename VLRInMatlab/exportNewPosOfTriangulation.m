function exportNewPosOfTriangulation( contourPoints, nextPos, dataName )
fileName = strcat( '../FinalVLRForMatlab/triangulation-', dataName, '.txt' );
fileId = fopen( char(fileName), 'a' );

nextPos = nextPos';
fprintf( fileId, '%4f %4f\n', nextPos( 1:2, : ) );

% then write the contour points
contourPoints = contourPoints';
fprintf( fileId, '%4f %4f\n', contourPoints( 1:2, : ) );

fprintf( fileId, '\n' );
fclose( fileId );