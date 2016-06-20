function [curS, curCellCenters, curMaxIntensities, timeS] =...
  applyCellSegmentation( t, tStart, tMax, tSteps, timeS, minVoxelCount,...
  voxelCountPerCell, anisotropyZ, inputPath, storeTIFF )
cellRadius = 15;
%   if endT ~= startT
% 		curTRatio = double(t-1)/double(abs( endT - startT ));
%   else
% 		curTRatio = 0;
%   end

%voxelCountPerCell = voxelCountPerCellStart + curTRatio*double(abs( voxelCountPerCellEnd - voxelCountPerCellStart ));
%msg = strcat( msg, {' '}, 'VoxelSize', {' '}, num2str(voxelCountPerCell) );
%disp( msg );

% apply cell segmentation based only on current time step
if tSteps == 0
  if t < 10
    digit = '00';
  elseif t < 100
    digit = '0';
  else
    digit = '';
  end
  fileName = strcat( inputPath, '\changed_t', digit, num2str(t), '.tif' );
  %fileName = strcat( inputPath, '\manThres_t', digit, num2str(t), '.tif' );
  membraneFileName = strcat( inputPath, '\changed_membrane_t', digit, num2str(t), '.tif' );
  newFileName = strcat( inputPath, '\preprocessed_t', digit, num2str(t), '.tif' );
  [curS, curCellCenters, curMaxIntensities, width, height, slices] =...
    findConnectedComponents( fileName, membraneFileName, minVoxelCount, voxelCountPerCell, anisotropyZ );
  if storeTIFF == 1
    generateCellShape( width, height, slices, curCellCenters, curMaxIntensities, cellRadius, newFileName );
  end
else
  % else apply a segmentation with respect to future segmentation results
  tEnd = t+tSteps;
  % TODO
%   if tEnd > tMax
%     tEnd = tMax;
%     tSteps = tEnd-tStart;
%   end
  
  if t == tStart
    % if it is the first time step then apply the segmentation for all coming tSteps
    numSteps = 1;
    for tt=tStart:tEnd
      if tt < 10
        digit = '00';
      elseif tt < 100
        digit = '0';
      else
        digit = '';
      end
      fileName = strcat( inputPath, '\changed_t', digit, num2str(tt), '.tif' );
      membraneFileName = strcat( inputPath, '\changed_membrane_t', digit, num2str(tt), '.tif' );
      [S, cellCenters, maxIntensities, width, height, slices] =...
        findConnectedComponents( fileName, membraneFileName,...
        minVoxelCount, voxelCountPerCell, anisotropyZ );
      timeS{ 1, numSteps } = S;
      timeS{ 2, numSteps } = cellCenters;
      numSteps = numSteps + 1;
    end
  else
    % only determine the connected components of the newly added time step
    tt = tEnd;
    if tt < 10
      digit = '00';
    elseif tt < 100
      digit = '0';
    else
      digit = '';
    end
    fileName = strcat( inputPath, '\changed_t', digit, num2str(tt), '.tif' );
    membraneFileName = strcat( inputPath, '\changed_membrane_t', digit, num2str(tt), '.tif' );
    [S, cellCenters, maxIntensities, width, height, slices] =...
      findConnectedComponents( fileName, membraneFileName,...
      minVoxelCount, voxelCountPerCell, anisotropyZ );
    timeS( 1:2, 1:tSteps ) = timeS( 1:2, 2:tSteps+1 );
    timeS{ 1, tSteps+1 } = S;
    timeS{ 2, tSteps+1 } = cellCenters;
  end
  
  % compare current segmentation results with future results and choose the
  % most adequate result for the current timestep
  r = 30;
  numCells = size( timeS{2,1}, 1 );
  changedStatus = cell( numCells, 1 );
  for c=1:numCells
    center =  timeS{2,1}( c, : )
    status = zeros( 1, tSteps );
    [IDX,D] = rangesearch( timeS{2,1}, timeS{2,1}( c, : ), r );
    numNeighbors = size(D{1,:},2)
    for s=1:tSteps
      [IDX,D] = rangesearch( timeS{2,s+1}, timeS{2,1}( c, : ), r );
      status( 1, s ) = size(D{1,:},2);
    end
    changedStatus{ c, 1 } = status;
    %status
    app = zeros( size( status ) );
    for i=1:length(status)
      app(i) = sum(numNeighbors==status(i));
    end
    app
  end
  
  curS = timeS{ 1, 1 };
  curCellCenters = timeS{ 2, 1 };
  % TODO
  curMaxIntensities = maxIntensities;
end

