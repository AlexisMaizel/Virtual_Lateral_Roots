function [S, cellCenters, maxIntensities] = identifyCellObjects( imageStack,...
  type, membraneFileName, minVoxelCount, voxelCountPerCell )
randomIntensities = 1;
cellRadius = 10;
if type == 0
  [S, cellCenters, maxIntensities] = mergeNearestCCs( imageStack, cellRadius, randomIntensities );
elseif type == 1
  [S, cellCenters, maxIntensities] = determineCellsBasedOnSizeAndMembrane(...
    imageStack, voxelCountPerCell, randomIntensities, minVoxelCount, membraneFileName );
elseif type == 2
  [S, cellCenters, maxIntensities] = determineLocalMaxima( imageStack, cellRadius, randomIntensities, membraneFileName );
end
