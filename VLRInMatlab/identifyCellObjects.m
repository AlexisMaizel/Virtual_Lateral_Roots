function [S, cellCenters, maxIntensities, width, height, slices, thresholdStack] =...
  identifyCellObjects( fileName, membraneFileName, minVoxelCount,...
  voxelCountPerCell )
imageStack = readTIFstack( char(fileName) );
% create 3D array of binary image data
height = size( imageStack, 1 );
width = size( imageStack, 2 );
slices = size( imageStack, 3 );

% type of method how to map a connected component (cc) onto a cell
% type == 0 -> nearest cc are merged based on cellRadius and represent a cell.
% type == 1 -> consider a user-chosen cell size to determine an initial
% number of cells for a cc which is then clustered and verified by
% considering cell walls in the membrane channel. Ignore ccs smaller than
% some threshold.
% type == 2 -> generate cc using regional maxima, merge ccs with distance
% smaller than a user-defined cellRadius, check for nearest neighbors and
% additionally check membrane channel.
type = 2;

randomIntensities = 1;
cellRadius = 10;
if type == 0
  [S, cellCenters, maxIntensities] = mergeNearestCCs( imageStack, cellRadius, randomIntensities );
elseif type == 1
  [S, cellCenters, maxIntensities] = determineCellsBasedOnSizeAndMembrane(...
    imageStack, voxelCountPerCell, randomIntensities, minVoxelCount, membraneFileName );
elseif type == 2
  [S, cellCenters, maxIntensities, thresholdStack] = determineLocalMaxima( imageStack, cellRadius, randomIntensities, membraneFileName );
end
