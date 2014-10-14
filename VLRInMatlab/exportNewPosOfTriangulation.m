function exportNewPosOfTriangulation( contourPoints, nextPos, dataName, autoContour )
if autoContour == 0
  fileName = strcat( '/tmp/triangulation-', dataName, '.txt' );
else
  fileName = strcat( '/tmp/triangulation-', dataName, '_auto.txt' );
end
fileId = fopen( char(fileName), 'a' );

nextPos = nextPos';
fprintf( fileId, '%4f %4f\n', nextPos( 1:2, : ) );

% then write the contour points
contourPoints = contourPoints';
fprintf( fileId, '%4f %4f\n', contourPoints( 1:2, : ) );

fprintf( fileId, '\n' );
fclose( fileId );