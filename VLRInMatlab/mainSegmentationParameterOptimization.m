cd('C:\Jens\VLRRepository\VLRInMatlab')
addpath( genpath( '/geom3d/' ) );
addpath( genpath( '/export_fig/' ) );
addpath( genpath( '/@tree/' ) );

chosenData = 3;
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };
rawDataStr = { '120830' '121204' '121211' '130508' '130607' };
inputPath = strcat( 'I:\SegmentationResults\Preprocessing\', rawDataStr( 1, chosenData ), '\changed_t' );
startT = 1;
endT = 150;
spatialNormalization = 1;
minVoxelCount = 300;
startCellSize = 2500;
endCellSize = 3500;
stepCellSize = 100;

storeTIFF = 1;
cellRadius = 17;

% output format of values
format longG

% read raw data (manual segmentation and tracking)
[ cellData, dimData, centerPosPerTimeStep, numCellsPerTimeStep ] =...
  readRawData( dataStr( 1, chosenData ) );

% start measuring elapsed time
%tic

for t=startT:endT
  msg = strcat( 'Timestep', {' '}, num2str(t) );
  disp( msg );
        
  % image output options
  if t < 10
    fileName = strcat( inputPath, '00', num2str(t), '_c1.tif' );
    digit = strcat( rawDataStr( 1, chosenData ), 'TimeStepSeg', '_00' );
    newFileName = strcat( inputPath, '00', num2str(t), '_c1M.tif' );
  elseif t < 100
    fileName = strcat( inputPath, '0', num2str(t), '_c1.tif' );
    digit = strcat( rawDataStr( 1, chosenData ), 'TimeStepSeg', '_0' );
    newFileName = strcat( inputPath, '0', num2str(t), '_c1M.tif' );
  else
    fileName = strcat( inputPath, num2str(t), '_c1.tif' );
    digit = strcat( rawDataStr( 1, chosenData ), 'TimeStepSeg', '_' );
    newFileName = strcat( inputPath, num2str(t), '_c1M.tif' );
  end
  
  % store cell positions of ground truth data in X
  numCellsManual = numCellsPerTimeStep(t, 1);
  cellCounter = 1;
  X = zeros(numCellsManual, 2);
  for j=1:dimData
    if cellData{j, 5} == t
      % get position of current cell
      p = [ cellData{j, 2} cellData{j, 3} cellData{j, 4} ];
      if spatialNormalization == 1
        p = p - centerPosPerTimeStep( t, : );
      end
      X(cellCounter, :) = [ p(1) p(2) ];
      cellCounter = cellCounter + 1;
    end
  end
  
  minDist = realmax;
  minParam = startCellSize;
  % loop over parameter range
  for curParam=startCellSize:stepCellSize:endCellSize
    msg = strcat( 'Param', {' '}, num2str(curParam) );
    disp( msg );
    % get connected components
    [S, cc, maxInt, imageStack] = findConnectedComponents( fileName, minVoxelCount, curParam );
    
    % store cell centers determined automatically
    centerPosAuto = zeros(1, 3);
    numCellsAuto = 0;
    for i=1:size(cc,1)
      cen = cc(i, :);
      centerPosAuto = centerPosAuto + cen;
      numCellsAuto = numCellsAuto+1;
    end
    centerPosAuto = centerPosAuto./numCellsAuto;
    
    % store all auto positions
    cellCounter = 1;
    Y = zeros(numCellsAuto, 2);
    for i=1:size(cc,1)
      cen = cc(i, :);
      if spatialNormalization == 1
        cen = cen - centerPosAuto;
      end
      Y(cellCounter, :) = [ cen(1) cen(2) ];
      cellCounter = cellCounter + 1;
    end
    
    % apply similarity check between ground truth data and automatic
    % segmentation result
    [ dist ] = determinePointSetSimilarity( X, Y );
    
    if storeTIFF == 1 && dist < minDist
      generateCellShape( imageStack, cc, maxInt, cellRadius, newFileName );
    end
    
    % update min values
    if dist < minDist
      minDist = dist;
      minParam = curParam;
    end
    
    % store current param settings
    saveSegmentationSimilarityValues( rawDataStr( 1, chosenData ), t, curParam,...
      endCellSize, dist, minParam, minDist, numCellsManual, numCellsAuto );
  end
end

% print elapsed time
%toc
