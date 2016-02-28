function exportCurrentMinError( dataStr, timeStep, minThres, minTau, minError, numCellsA, numCellsM )
fileName = strcat( 'ParameterOptimization\MinErrorAnalysis', dataStr, '.txt' );
fileID = fopen( char(fileName), 'a' );
fprintf( fileID, '%d %.2f %.2f %.2f %d %d\n', timeStep, minThres, minTau, minError, numCellsA, numCellsM );
fclose( fileID );