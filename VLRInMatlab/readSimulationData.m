function [ cellData, connMat ] = readSimulationData( dataStr, simulationSteps, numTimeSteps )
% reading simulation data
path = strcat( 'I:\GraphTopologyAnalysis\SimulationData\AdjacencyMatrices', char(dataStr), '.csv' );
pathCP = strcat( 'I:\GraphTopologyAnalysis\SimulationData\AdjacencyMatricesCellProperties', char(dataStr), '.csv' );
fileID = fopen( char(path) );
fileIDCP = fopen( char(pathCP) );

% initialze a cell array for number of simulation steps
connMat = cell( simulationSteps, 1 );
cellData = cell( simulationSteps, 1 );

% read adjacency matrix information with format:
% first line of each simulation step: SimulationStep TimeStep NumCells \n
% then the (NumCells x NumCells) adjacency matrix
tline = fgetl( fileID );
for s=1:simulationSteps
  connMat{ s, 1 } = cell( numTimeSteps, 1 );
  for t=1:numTimeSteps
    % first read header of current time and simulation step
    h = str2num( tline );
    numCells = h( 1, 3 );
    % then read adjacency matrix
    adjMat = zeros( numCells, numCells );
    for i=1:numCells
      tline = fgetl( fileID );
      adjMat( i, : ) = str2num( tline );
    end
    connMat{ s, 1 }{ t, 1 } = adjMat;
    tline = fgetl( fileID );
  end
end
fclose(fileID);

% read cell Data information with format:
% first line of each simulation step: SimulationStep TimeStep NumCells \n
% uniqueIDPerTimeStep uniqueGlobalId parentID posx posy posz layerValue \n
tline = fgetl( fileIDCP );
for s=1:simulationSteps
  cellData{ s, 1 } = cell( numTimeSteps, 1 );
  for t=1:numTimeSteps
    % first read header of current time and simulation step
    h = str2num( tline );
    numCells = h( 1, 3 );
    % then read cell data information
    data = zeros( numCells, 7 );
    for i=1:numCells
      tline = fgetl( fileIDCP );
      data( i, : ) = str2num( tline );
    end
    cellData{ s, 1 }{ t, 1 } = data;
    tline = fgetl( fileIDCP );
  end
end
fclose(fileIDCP);