function highlightError( dataStr, timeStep, paramT, paramTau )
fileName = strcat( 'ParameterOptimization\ErrorAnalysis', dataStr, '.txt' );
fileID = fopen( char(fileName), 'a' );
fprintf( fileID, '%s %d %.2f %.2f\n', 'An error occurred for', timeStep, paramT, paramTau );
fclose( fileID );