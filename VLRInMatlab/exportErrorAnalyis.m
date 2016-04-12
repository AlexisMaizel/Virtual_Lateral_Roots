function exportErrorAnalyis( dataStr, timeStep, paramT, paramTau, error, numCellsA, numCellsM )
%mkdir( 'ParameterOptimization\' );
fileName = strcat( 'ParameterOptimization\ErrorAnalysis', dataStr, '.txt' );
fileID = fopen( char(fileName), 'a' );
date = datestr( datetime('now') );
fprintf( fileID, '%d %.2f %.2f %.2f %d %d %s\n', timeStep, paramT, paramTau, error, numCellsA, numCellsM, date);
fclose( fileID );