% Complete nuclei segmentation step from, thresholding, different methods for identifying cell objects
% as using the membrane channel, nearest neighbor search or clustering of regional maxima.
setWorkingPathProperties()

chosenData = 8;
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };
rawDataStr = { '120830' '121204' '121211' '130508' '130607' '20160427' '20160428' '20160426' };
startT = 10;
endT = 10;
radEllip = 5;%10
lineWidth = 1;
if chosenData < 6
  numPlots = 5;
else
  numPlots = 4;
end
numColorsCells = 40;
cmapNumCells = colorcube(numColorsCells);
ELLIP = [];
ELLIPPATCH = [];
nucleiCounter = 1;
spatialNormalization = 0;

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

% two parameters for finding connected components (cc) which are only
% considered when type == 1
% number of voxels a cc has to have else it will be ignored
minVoxelCount = 150;
% number of voxels for a cell in average
voxelCountPerCell = 10000;

% store TIFF of perfectly created cell shapes?
storeTIFF = 1;

% store png output
storePNGs = 1;

% store cell positions into csv file
storeCSVResult = 0;

% output format of values
format shortG

% input path
inputPath = strcat( 'I:\SegmentationResults\Preprocessing\', rawDataStr( 1, chosenData ) );
% path to image output
imageOutputPath = strcat( 'I:\SegmentationResults\Matlab\Segmentation\', rawDataStr( 1, chosenData ), '\' );
mkdir( char(imageOutputPath) );

% read raw data (manual segmentation and tracking) of Daniel's data
if chosenData < 6
  [ cellData, dimData, centerPosPerTimeStep, numCellsPerTimeStep ] =...
    readRawData( dataStr( 1, chosenData ) );
end

if storePNGs == 1
  if chosenData < 6
    w = 1800;
    h = 900;
  elseif chosenData == 7
    w = 1000;
    h = 750;
  else
    w = 1500;
    h = 750;
  end
  f = figure( 'Name', 'Segmentation', 'Position', [ 50 50 w h ] );
end

% start measuring elapsed time
tic
for t=startT:endT
  if t < 10
    digit = '00';
  elseif t < 100
    digit = '0';
  else
    digit = '';
  end
  
  %nucleiFileName = strcat( inputPath, '\changed_t', digit, num2str(t), '.tif' );
  %nucleiFileName = strcat( inputPath, '\manThres_t', digit, num2str(t), '.tif' );
  %nucleiFileName = strcat( 'I:\NewDatasets\Zeiss\20160427\red\cropped_spim_TL_slice', digit, num2str(t), '_Angle1.tif' );
  nucleiFileName = strcat( 'I:\NewDatasets\ilastikWorkshopData\20160427\nuclei\small_cropped_nuclei_T', digit, num2str(t), '_Angle1.tif' );
  %nucleiFileName = strcat( 'I:\NewDatasets\ilastikWorkshopData\20160427\nuclei\MAX_small_cropped_nuclei_T', digit, num2str(t), '_Angle1.tif' );
  
  %membraneFileName = strcat( inputPath, '\changed_membrane_t', digit, num2str(t), '.tif' );
  membraneFileName = strcat( 'I:\SegmentationResults\Matlab\Segmentation\20160427\Membrane\20160427_Membrane_T', digit, num2str(t), '.tif' );
  %membraneFileName = strcat( 'I:\NewDatasets\Zeiss\20160427\green\cropped_spim_TL_slice', digit, num2str(t), '_Angle1.jpg' );
  %membraneFileName = strcat( 'I:\NewDatasets\ilastikWorkshopData\20160427\membrane\small_cropped_membrane_T', digit, num2str(t), '_Angle1.tif' );
  
  % create 3D array of binary image data
  imageStack = readTIFstack( char(nucleiFileName) );
  height = size( imageStack, 1 );
  width = size( imageStack, 2 );
  slices = size( imageStack, 3 );
  
  %nImageStack = normalizeImage( imageStack, 16 );
  outputFileName = strcat( inputPath, '\preprocessed_t', digit, num2str(t), '.tif' );
  
  % generate gradient magnitude image
  %[ gradientStack ] = generate3DGradientMagnitudeImage( imageStack );
  %newFileName = strcat( inputPath, '\gradient_', digit, num2str(t), '.tif' );
  %writeTIFstack( gradientStack, char(newFileName), 2^31 );
  
  % apply a thresholding
  [ thresholdStack ] = applyThresholding( imageStack );
  
  % identify cell objects from the preprocessed threshold stack
  [S, cc, intensities] = identifyCellObjects( thresholdStack, type, membraneFileName, minVoxelCount, voxelCountPerCell );
  
  % store the extracted cells as newly perfect generated cells in tiff
  if storeTIFF == 1
    segmentedStack = generateCellShape( width, height, slices, cc, intensities, outputFileName );
  end
  
  if storePNGs == 1
    setSubPlots( f, numPlots, chosenData, spatialNormalization, width, height, 0 );
    
    % determine center of data points
    centerPosAuto = zeros(1, 3);
    numCellsAuto = 0;
    for i=1:size(cc,1)
      cen = cc(i, :);
      centerPosAuto = centerPosAuto + cen;
      numCellsAuto = numCellsAuto+1;
    end
    if numCellsAuto ~= 0
      centerPosAuto = centerPosAuto./numCellsAuto;
    end
    
    % store all auto positions
    cellCounter = 1;
    X = zeros(numCellsAuto, 2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% original regional maxima %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chosenData > 5
      ax1 = subplot( 1, numPlots, 3 );
    else
      ax1 = subplot( 2, (numPlots+1)/2, 1 );
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
    
    tit = strcat( 'LocalMaxima\_', 'T', {' '}, num2str(t), ',',...
      {' '}, '#CCs', {' '}, num2str(size(S,1)) );
    title( char(tit) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% automatic segmentation result %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chosenData > 5
      ax2 = subplot( 1, numPlots, 4 );
    else
      ax2 = subplot( 2, (numPlots+1)/2, 2 );
    end
    hold on
    
    for i=1:size(cc,1)
      cen = cc(i, :);
      if spatialNormalization == 1
        cen = cen - centerPosAuto;
      else
        cen = [ cen(1) -cen(2)+height cen(3) ];
      end
      X(cellCounter, :) = [ cen(1) cen(2) ];
      
      color = cmapNumCells( mod( cellCounter, numColorsCells )+1, : );
      [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
        drawEllipse3d( cen(1), cen(2), cen(3), radEllip, radEllip, 0, 0 );
      set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
      set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
      nucleiCounter = nucleiCounter+1;
      
      cellCounter = cellCounter + 1;
    end
    
    if storeCSVResult == 1
      csvPath = strcat( imageOutputPath, rawDataStr( 1, chosenData ), '_T', digit, num2str(t), '.dat' );
      csvwrite( char(csvPath), X );
    end
    
    clear title
    tit = strcat( 'AutoSeg\_', 'T', {' '}, num2str(t), ',',...
      {' '}, '#Cells', {' '}, num2str(numCellsAuto) );
    title( char(tit) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MIP image of raw data %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chosenData > 5
      ax3 = subplot( 1, numPlots, 1 );
    else
      ax3 = subplot( 2, (numPlots+1)/2, 3 );
    end
    
    h1 = showMIP( imageStack );
    clear title
    tit = strcat( 'rawMIP\_', 'T', {' '}, num2str(t) );
    title( char(tit) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MIP image of thresholding %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chosenData > 5
      ax4 = subplot( 1, numPlots, 2 );
    else
      ax4 = subplot( 2, (numPlots+1)/2, 4 );
    end
    h2 = showMIP( thresholdStack );
    clear title
    tit = strcat( 'thresholdMIP\_', 'T', {' '}, num2str(t) );
    title( char(tit) );
    
    if chosenData < 6
      %%% manual segmentation result %%%
      ax6 = subplot( 2, numPlots/2, 6 );
      hold on
      
      % store all manual positions
      numCellsManual = numCellsPerTimeStep(t, 1);
      cellCounter = 1;
      
      for j=1:dimData
        if cellData{j, 5} == t
          % get position of current cell
          p = [ cellData{j, 2} cellData{j, 3} cellData{j, 4} ];
          if spatialNormalization == 1
            p = p - centerPosPerTimeStep( t, : );
          end
          
          color = cmapNumCells( mod( cellCounter, numColorsCells )+1, : );
          [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
            drawEllipse3d( p(1), p(2), p(3), radEllip, radEllip, 0, 0 );
          set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
          set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
          nucleiCounter = nucleiCounter+1;
          
          cellCounter = cellCounter + 1;
        end
      end
      
      clear title
      tit = strcat( 'ManualSeg\_', 'T', {' '}, num2str(t), ',',...
        {' '}, '#Cells', {' '}, num2str(numCellsManual) );
      title( char(tit) );
    end
    filePath = strcat( imageOutputPath, rawDataStr( 1, chosenData ), '_T', digit, num2str(t), '.png' );
    export_fig( gcf, char(filePath), '-m2', '-png' );
  end
end
% print elapsed time
toc
