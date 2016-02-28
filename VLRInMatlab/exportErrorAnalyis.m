function exportErrorAnalyis( dataStr, timeStep, paramT, paramTau, error, numCellsA, numCellsM )
%mkdir( 'ParameterOptimization\' );
fileName = strcat( 'ParameterOptimization\ErrorAnalysis', dataStr, '.txt' );
fileID = fopen( char(fileName), 'a' );
fprintf( fileID, '%d %.2f %.2f %.2f %d %d\n', timeStep, paramT, paramTau, error, numCellsA, numCellsM );
fclose( fileID );