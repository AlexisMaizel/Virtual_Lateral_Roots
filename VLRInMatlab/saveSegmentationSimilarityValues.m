function saveSegmentationSimilarityValues( dataStr, timeStep, curParam, endParam, dist, minParam, minDist, numCellsX, numCellsY )
%mkdir( 'ParameterOptimization\' );
cd('I:\SegmentationResults\ParameterOptimization')
fileName = strcat( 'SimilarityValues_', dataStr, '.txt' );
fileID = fopen( char(fileName), 'a' );
date = datestr( datetime('now') );
fprintf( fileID, 'T: %d P: %.1f D: %.1f NX: %d NY: %d Date: %s\n', timeStep, curParam, dist, numCellsX, numCellsY, date);
if curParam == endParam
  fprintf( fileID, 'MinP: %.1f MinD: %.1f\n\n', minParam, minDist);
else
  fprintf( fileID, 'MinP: %.1f MinD: %.1f\n', minParam, minDist);
end
fclose( fileID );
cd('C:\Jens\VLRRepository\VLRInMatlab')