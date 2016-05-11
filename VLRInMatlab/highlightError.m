function highlightError( dataStr, startT, endT, paramT, paramTau )
fileName = strcat( 'ParameterOptimization\ErrorAnalysis', dataStr, '.txt' );
fileID = fopen( char(fileName), 'a' );
fprintf( fileID, '%s %d %d %.2f %.2f\n', 'An error occurred for', startT, endT, paramT, paramTau );
fclose( fileID );