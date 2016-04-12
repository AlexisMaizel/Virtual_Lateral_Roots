cd('C:\Jens\VLRRepository\VLRInMatlab')
addpath( genpath( '/geom3d/' ) );
addpath( genpath( '/export_fig/' ) );
addpath( genpath( '/@tree/' ) );

chosenData = 3;
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };
rawDataStr = { '120830' '121204' '121211' '130508' '130607' };
inputPath = strcat( 'I:\SegmentationResults\Preprocessing\', rawDataStr( 1, chosenData ), '\changed_t' );
startT = 150;
endT = 150;
radEllip = 10;
lineWidth = 1;
numPlots = 4;
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
minVoxelCount = 300;
voxelCountPerCell = 3300;

storeTIFF = 1;
storePNGs = 1;
cellRadius = 17;

% output format of values
format longG

% path to image output
imageDir = strcat( 'I:/SegmentationResults/Matlab/Segmentation/' );
mkdir( char(imageDir) );

% read raw data (manual segmentation and tracking)
[ cellData, dimData, centerPosPerTimeStep, numCellsPerTimeStep ] =...
  readRawData( dataStr( 1, chosenData ) );

%generateTreeStructureFromData( cellData, size(numCellsPerTimeStep,1), 1, 20 );

if storePNGs == 1
  f = figure( 'Name', 'Segmentation', 'Position', [ 50 50 1600 800 ] );
end

% start measuring elapsed time
tic

for t=startT:endT
  msg = strcat( 'Timestep', {' '}, num2str(t) );
  disp( msg );
  if storePNGs == 1
    % subplot settings
    clf(f)
    for p=1:numPlots
      subplot( numPlots/2, numPlots/2, p );
      hold on
      if spatialNormalization == 0
        xMinMax = [ -50 750 ];
        if p == numPlots
          yMinMax = [ 0 400 ];
        else
          yMinMax = [ -450 -50 ];
        end
      else
        xMinMax = [ -400 400 ];
        yMinMax = [ -200 200 ];
      end
      % axis([xmin xmax ymin ymax zmin zmax cmin cmax])
      axis( [ xMinMax(1) xMinMax(2) yMinMax(1) yMinMax(2) -10000 10000 0 1 ] );
      axis on
      daspect( [ 1 1 1 ] );
      xlabel('X');
      ylabel('Y');
      zlabel('Z');
      camproj( 'orthographic' );
    end
  end
        
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
  
  % get connected components
  [S, cc, maxInt, imageStack] = findConnectedComponents( fileName, minVoxelCount, voxelCountPerCell );
  if storeTIFF == 1
    generateCellShape( imageStack, cc, maxInt, cellRadius, newFileName );
  end
  
  if storePNGs == 1
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
    
    %%% automatic segmentation result before clustering %%%
    ax1 = subplot( numPlots/2, numPlots/2, 1 );
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
      end
      
      color = cmapArea( mod( areaSize-1, numColorsArea )+1, : );
      [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
        drawEllipse3d( cen(1), -cen(2), cen(3), radEllip, radEllip, 0, 0 );
      set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
      set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
      nucleiCounter = nucleiCounter+1;
    end
    
    tit = strcat( 'AutoSegBefore\_', 'T', {' '}, num2str(t), ',',...
      {' '}, '#Cells', {' '}, num2str(size(S,1)) );
    title( char(tit) );
    
    colormap(ax1, cmapArea)
    colorbar
    caxis([0, numColorsArea])
    
    %%% automatic segmentation result before clustering without small ccs %%%
    ax2 = subplot( numPlots/2, numPlots/2, 2 );
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
      end
      
      color = cmapArea( mod( areaSize-1, numColorsArea )+1, : );
      [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
        drawEllipse3d( cen(1), -cen(2), cen(3), radEllip, radEllip, 0, 0 );
      set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
      set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
      nucleiCounter = nucleiCounter+1;
    end
    
    tit = strcat( 'AutoSegBeforeCropped\_', 'T', {' '}, num2str(t), ',',...
      {' '}, '#Cells', {' '}, num2str(remainingCells) );
    title( char(tit) );
    
    colormap(ax2, cmapArea)
    colorbar
    caxis([0, numColorsArea])
    
    %%% automatic segmentation result after clustering %%%
    ax3 = subplot( numPlots/2, numPlots/2, 3 );
    hold on
    
    for i=1:size(cc,1)
      cen = cc(i, :);
      if spatialNormalization == 1
        cen = cen - centerPosAuto;
      end
      X(cellCounter, :) = [ cen(1) cen(2) ];
      
      color = cmapNumCells( mod( cellCounter, numColorsCells )+1, : );
      [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
        drawEllipse3d( cen(1), -cen(2), cen(3), radEllip, radEllip, 0, 0 );
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
    
    %%% manual segmentation result %%%
    ax4 = subplot( numPlots/2, numPlots/2, 4 );
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
    
    %colormap(ax4, cmapNumCells)
    %colorbar
    %caxis([0, numColorsCells])
    
    filePath = strcat( imageDir, digit, num2str(t), '.png' );
    export_fig( gcf, char(filePath), '-m2', '-png' );
  end
end

% print elapsed time
toc
