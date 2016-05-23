setWorkingPathProperties()

chosenData = 7;
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };
rawDataStr = { '120830' '121204' '121211' '130508' '130607' '20160427' '20160428' };
inputPath = strcat( 'I:\SegmentationResults\Preprocessing\', rawDataStr( 1, chosenData ), '\changed_t' );
startT = 0;
endT = 34;
radEllip = 10;
lineWidth = 1;
if chosenData < 6
  numPlots = 6;
else
  numPlots = 5;
end
numColorsCells = 40;
cmapNumCells = colorcube(numColorsCells);
ELLIP = [];
ELLIPPATCH = [];
nucleiCounter = 1;
spatialNormalization = 0;

% parameters for finding connected components (cc)
% both parameters have to be chosen carefully depending on
% the choice of threshold used in the Fiji script for thresholding the
% image
% thres = 1000 -> minV = 200, voxelC = 2700
% thres = 900 -> minV = 400, voxelC = 3200
% thres = 750 -> minV = 200, voxelC = 3200
minVoxelCount = 800;
voxelCountPerCellStart = 8000;%4600;
voxelCountPerCellEnd = 8000;%4600;

storeTIFF = 1;
storePNGs = 1;
cellRadius = 15;

% output format of values
format shortG %longG %shortG

% path to image output
imageDir = strcat( 'I:/SegmentationResults/Matlab/Segmentation/', rawDataStr( 1, chosenData ), '/' );
mkdir( char(imageDir) );

% read raw data (manual segmentation and tracking)
if chosenData < 6
  [ cellData, dimData, centerPosPerTimeStep, numCellsPerTimeStep ] =...
    readRawData( dataStr( 1, chosenData ) );
end

%generateTreeStructureFromData( cellData, size(numCellsPerTimeStep,1), 1, 20 );

if storePNGs == 1
  f = figure( 'Name', 'Segmentation', 'Position', [ 50 50 1400 600 ] );
end

% start measuring elapsed time
tic

for t=startT:endT
  msg = strcat( 'Timestep', {' '}, num2str(t) );
        
  % image output options
  if t < 10
    fileName = strcat( inputPath, '00', num2str(t), '.tif' );
    digit = strcat( rawDataStr( 1, chosenData ), 'TimeStepSeg', '_00' );
    newFileName = strcat( inputPath, '00', num2str(t), '_M.tif' );
  elseif t < 100
    fileName = strcat( inputPath, '0', num2str(t), '.tif' );
    digit = strcat( rawDataStr( 1, chosenData ), 'TimeStepSeg', '_0' );
    newFileName = strcat( inputPath, '0', num2str(t), '_M.tif' );
  else
    fileName = strcat( inputPath, num2str(t), '.tif' );
    digit = strcat( rawDataStr( 1, chosenData ), 'TimeStepSeg', '_' );
    newFileName = strcat( inputPath, num2str(t), '_M.tif' );
  end
  
  % get connected components
  if endT ~= startT
		curTRatio = double(t-1)/double(abs( endT - startT ));
  else
		curTRatio = 0;
  end
  
	voxelCountPerCell = voxelCountPerCellStart + curTRatio*double(abs( voxelCountPerCellEnd - voxelCountPerCellStart ));
  msg = strcat( msg, {' '}, 'VoxelSize', {' '}, num2str(voxelCountPerCell) );
  %disp( msg );
  [S, cc, maxInt, imageStack] = findConnectedComponents( fileName, minVoxelCount, voxelCountPerCell );
  if storeTIFF == 1
    generateCellShape( imageStack, cc, maxInt, cellRadius, newFileName );
  end
  
  if storePNGs == 1
    height = size( imageStack, 1);
    width = size( imageStack, 2);
    setSubPlots( f, numPlots, chosenData, spatialNormalization, width, height );
    
    % determine center of data points
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
    X = zeros(numCellsAuto, 2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% automatic segmentation result before clustering %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chosenData < 6
      ax1 = subplot( numPlots/2, numPlots/2, 1 );
    else
      ax1 = subplot( 1, numPlots, 1 );
    end
    hold on
    
    numColorsArea = 0;
    for i=1:size(S,1)
      areaSize = S(i, :).Area;
      if areaSize > numColorsArea
        numColorsArea = areaSize;
      end
    end
    cmapArea = cool(numColorsArea);
    for i=1:size(S,1)
      areaSize = S(i, :).Area;
      cen = S(i, :).Centroid;
      if spatialNormalization == 1
        cen = cen - centerPosAuto;
      else
        cen = [ cen(1) -cen(2)+height cen(3) ];
      end
      
      color = cmapArea( mod( areaSize-1, numColorsArea )+1, : );
      [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
        drawEllipse3d( cen(1), cen(2), cen(3), radEllip, radEllip, 0, 0 );
      set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
      set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
      nucleiCounter = nucleiCounter+1;
    end
    
    tit = strcat( 'AutoSegBefore\_', 'T', {' '}, num2str(t), ',',...
      {' '}, '#Cells', {' '}, num2str(size(S,1)) );
    title( char(tit) );
    
%     colormap(ax1, cmapArea)
%     colorbar
%     caxis([0, numColorsArea])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% automatic segmentation result before clustering without small ccs %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chosenData < 6
      ax2 = subplot( numPlots/2, numPlots/2, 2 );
    else
      ax2 = subplot( 1, numPlots, 2 );
    end
    hold on
    
    remainingCells = 0;
    for i=1:size(S,1)
      areaSize = S(i, :).Area;
      
      if areaSize < minVoxelCount
        continue;
      else
        remainingCells = remainingCells + 1;
      end
      
      cen = S(i, :).Centroid;
      if spatialNormalization == 1
        cen = cen - centerPosAuto;
      else
        cen = [ cen(1) -cen(2)+height cen(3) ];
      end
      
      color = cmapArea( mod( areaSize-1, numColorsArea )+1, : );
      [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
        drawEllipse3d( cen(1), cen(2), cen(3), radEllip, radEllip, 0, 0 );
      set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
      set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
      nucleiCounter = nucleiCounter+1;
    end
    
    tit = strcat( 'AutoSegBeforeCropped\_', 'T', {' '}, num2str(t), ',',...
      {' '}, '#Cells', {' '}, num2str(remainingCells) );
    title( char(tit) );
    
%     colormap(ax2, cmapArea)
%     colorbar
%     caxis([0, numColorsArea])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% automatic segmentation result after clustering %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chosenData < 6
      ax3 = subplot( numPlots/2, numPlots/2, 3 );
    else
      ax3 = subplot( 1, numPlots, 3 );
    end
    hold on
    
    for i=1:size(cc,1)
      cen = cc(i, :);
      if spatialNormalization == 1
        cen = cen - centerPosAuto;
      else
        cen = [ cen(1) -cen(2)+height cen(3) ];
      end
      X(cellCounter, :) = [ cen(2) cen(1) ];
      
      color = cmapNumCells( mod( cellCounter, numColorsCells )+1, : );
      [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
        drawEllipse3d( cen(1), cen(2), cen(3), radEllip, radEllip, 0, 0 );
      set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
      set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
      nucleiCounter = nucleiCounter+1;
      
      cellCounter = cellCounter + 1;
    end
    
    tit = strcat( 'AutoSegAfter\_', 'T', {' '}, num2str(t), ',',...
      {' '}, '#Cells', {' '}, num2str(numCellsAuto) );
    title( char(tit) );
    
    %colormap(ax3, cmapNumCells)
    %colorbar
    %caxis([0, numColorsCells])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MIP image of raw data %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chosenData < 6
      ax4 = subplot( numPlots/2, numPlots/2, 4 );
    else
      ax4 = subplot( 1, numPlots, 4 );
    end
    
    h1 = showMIP( rawDataStr( 1, chosenData ), 0, t );
    tit = strcat( 'rawMIP\_', 'T', {' '}, num2str(t) );
    title( char(tit) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MIP image of thresholding %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chosenData < 6
      ax5 = subplot( numPlots/2, numPlots/2, 5 );
    else
      ax5 = subplot( 1, numPlots, 5 );
    end
    h2 = showMIP( rawDataStr( 1, chosenData ), 1, t );
    tit = strcat( 'thresholdMIP\_', 'T', {' '}, num2str(t) );
    title( char(tit) );
    
    if chosenData < 6
      %%% manual segmentation result %%%
      ax6 = subplot( numPlots/2, numPlots/2, 6 );
      hold on
      
      % store all manual positions
      numCellsManual = numCellsPerTimeStep(t, 1);
      cellCounter = 1;
      Y = zeros(numCellsManual, 2);
      
      for j=1:dimData
        if cellData{j, 5} == t
          % get position of current cell
          p = [ cellData{j, 2} cellData{j, 3} cellData{j, 4} ];
          if spatialNormalization == 1
            p = p - centerPosPerTimeStep( t, : );
          end
          Y(cellCounter, :) = [ p(1) p(2) ];
          
          color = cmapNumCells( mod( cellCounter, numColorsCells )+1, : );
          [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
            drawEllipse3d( p(1), p(2), p(3), radEllip, radEllip, 0, 0 );
          set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
          set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
          nucleiCounter = nucleiCounter+1;
          
          cellCounter = cellCounter + 1;
        end
      end
      
      tit = strcat( 'ManualSeg\_', 'T', {' '}, num2str(t), ',',...
        {' '}, '#Cells', {' '}, num2str(numCellsManual) );
      title( char(tit) );
      
      %colormap(ax6, cmapNumCells)
      %colorbar
      %caxis([0, numColorsCells])
    end
    
    filePath = strcat( imageDir, digit, num2str(t), '.png' );
    export_fig( gcf, char(filePath), '-m2', '-png' );
  end
end

% print elapsed time
toc
