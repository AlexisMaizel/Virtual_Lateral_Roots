function exportCurrentMinError( dataStr, timeStep, minThres, minTau, minError, numCellsA, numCellsM )
fileName = strcat( 'ParameterOptimization\MinErrorAnalysis', dataStr, '.txt' );
fileID = fopen( char(fileName), 'a' );
date = datestr( datetime('now') );
fprintf( fileID, '%d %.2f %.2f %.2f %d %d %s\n', timeStep, minThres, minTau, minError, numCellsA, numCellsM, date );
fclose( fileID );