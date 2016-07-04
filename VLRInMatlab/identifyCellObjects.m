function [S, cellCenters, maxIntensities, width, height, slices] =...
  identifyCellObjects( fileName, membraneFileName, minVoxelCount,...
  voxelCountPerCell )
imageStack = readTIFstack( char(fileName) );
% create 3D array of binary image data
height = size( imageStack, 1 );
width = size( imageStack, 2 );
slices = size( imageStack, 3 );

% type of method how to map a connected component (cc) onto a cell
% type == 0 -> nearest cc are merged and these represent a cell
% type == 1 -> consider a user-chosen cell size to determine an initial
% number of cells a cc consists of which is then clustered and verified by
% considering cell walls in the membrane channel. Ignore ccs smaller than
% some threshold
% type == 2 -> consider a distance mask and its number of local maxima to
% check if a cc consists of one or more cells, also checking membrane
% channel. Ignore cc smaller than some threshold
type = 1;

randomIntensities = 1;
cellRadius = 50;
if type == 0
  [S, cellCenters, maxIntensities] = mergeNearestCCs( imageStack, cellRadius, randomIntensities );
elseif type == 1
  [S, cellCenters, maxIntensities] = determineCellsBasedOnSizeAndMembrane(...
    imageStack, voxelCountPerCell, randomIntensities, minVoxelCount, membraneFileName );
elseif type == 2
  [S, cellCenters, maxIntensities] = determineRegionalMaxima( imageStack,...
    cellRadius, randomIntensities, membraneFileName );
end
