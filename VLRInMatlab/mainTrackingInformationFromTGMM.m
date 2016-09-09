% execute compileMex.m before to generate mex files
setWorkingPathProperties()

% output format of values
format longG

%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
requireNewSegmentation = 1;
requireNewTracking = 1;

% render specific properties
renderAutoLineages = 0;
renderManLineages = 0;
renderNuclei = 1;
showAdditionalImages = 1;

% parameters for TGMM
paramThreshold = 1000; % 600 for all preprocessed data sets
% the higher paramTau the fewer cells are merged to one cell
% 200 for 121211
% 20 for chosenData > 5
paramTau = 200;
minNucleiSize = 10;
maxNucleiSize = 3000;

% which kind of data should be chosen: raw (=0) or preproccessed (=1)
dataType = 0;
startT = 1;
endT = 50;
chosenData = 3;
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };
rawDataStr = { '120830' '121204' '121211' '130508' '130607' '20160427' '20160428' '20160426' };
% normalize data points based on their center
spatialNormalization = 0;
resultPath = strcat( 'I:\SegmentationResults\TGMM\', rawDataStr( 1, chosenData ), '\' );
radEllip = 10;% 5, 10
lineWidth = 1;
numColors = 40;
% amount of plots
if chosenData < 6
  numPlots = 4;
else
  numPlots = 3;
end
if showAdditionalImages == 0
  numPlots = numPlots - 2;
end

% colormap for ellipses
cmap = colorcube(numColors);
ELLIP = [];
ELLIPPATCH = [];
nucleiCounter = 1;
lineageAutoMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
lineageManualMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
lineageManualColorMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%

% path to image output
imageDir = strcat( 'I:/SegmentationResults/Matlab/Tracking/', rawDataStr( 1, chosenData ), '/' );
mkdir( char(imageDir) );

% start measuring elapsed time
tic

% read raw data (manual segmentation and tracking)
if chosenData < 6
  [ cellData, dimData, centerPosPerTimeStep, numCellsPerTimeStep ] =...
    readRawData( dataStr( 1, chosenData ) );
end

% execute TGMM segmentation and tracking
tRange = strcat( num2str(startT), {' '}, num2str(endT) );
dataSelect = strcat( 'TGMM_configFile', rawDataStr( 1, chosenData ), '.txt' );
% first create TGMMconfigfile
rawDataPath = exportTGMMConfigFile( rawDataStr( 1, chosenData ), paramThreshold, paramTau, minNucleiSize, maxNucleiSize, dataType );
cd('C:\Jens\TGMM_Supplementary_Software_1_0\build')
if requireNewSegmentation == 1
  cmdSegmentation = strcat( 'nucleiChSvWshedPBC\Release\ProcessStackBatchMulticore.exe AutoTGMMConfig\', dataSelect, {' '}, tRange );
  system( char(cmdSegmentation), '-echo');
  disp( 'Segmentation done.' );
end
if requireNewTracking == 1
  % execute tracking
  cmdTracking = strcat( 'Release\TGMM AutoTGMMConfig\', dataSelect, {' '}, tRange );
  system( char(cmdTracking), '-echo');
  disp( 'Tracking done.' );
end
cd('C:\Jens\VLRRepository\VLRInMatlab')

newestResult = getNewestFolder( char( resultPath ) );
totalPath = strcat( resultPath, newestResult, '\XML_finalResult_lht\GMEMfinalResult_frame');

% read TGMM data (automatic segmentation and tracking)
% svIdxCell:		cell array of length N, where N is the total number of objects tracked. svIdxCell{i} contains the indexes of the supervoxels belonging to the i-th object in trackingMatrix array. This index is necessary to rescue the segmentation from the .svb files output by TGMM software. The index starts in 0 following C convention.
% trackingMatrix:		numericall array of size Nx10, where N is the number of points tracked bt TGMM over time. Each of the columns contains the following information:
%
% 1.	Unique Id from the database to identify the point ( a large integer number)
% 2.	Cell type (represented by an integer). It is 0 if no cell type has been set for this object.
% 3.	x location of the nucleus centroid in world coordinates. Use the variable stackRes to convert from world coordinates to pixel unites.
% 4.	Same as 3 but for y location.
% 5.	Same as 3 but for z location.
% 6.	Estimated radius of the nucleus. It is 0 if this parameter was not estimated.
% 7.	Id of the cell in the previous time point. It is -1 if there is no linkage. Otherwise it has the unique id of the parent from column 1, so you can reconstruct the lineage.
% 8.	Time point of the nucleus.
% 9.	Confidence level in the tracking result. Value of 3 indicates high confidence that the object was correctly tracked. Value of 0 indicates low confidence.
% 10.	Skeleton id. All cells belonging to the same lineage have the same unique skeleton id.
[trackingMatrix, svIdxCell, errorOccurred] = parseMixtureGaussiansXml2trackingMatrixCATMAIDformat( char(totalPath), startT, endT );

if errorOccurred == 1
  highlightError( rawDataStr( 1, chosenData ), startT, endT, paramThreshold, paramTau );
end

% figure settings
if renderNuclei == 1
  f = figure( 'Name', 'Tracking', 'Position', [ 50 50 1600 600 ] );
end

maxAutoLineages = max( trackingMatrix( :, 10 ) );
autoNodes = cell( maxAutoLineages, 1 );
cellAutoNumLin = ones( maxAutoLineages, 1 );
autoNodeIdToNumCellMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );

if chosenData < 6
  maxManLineages = max( [ cellData{:, 6} ] );
  manNodes = cell( maxManLineages, 1 );
  cellManNumLin = ones( maxManLineages, 1 );
  manNodeIdToNumCellMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
end

imageStack = readTIFstack( char(rawDataPath) );
height = size( imageStack, 1);
width = size( imageStack, 2);

for t=startT:endT
  % subplot settings
  if renderNuclei == 1
    setSubPlots( f, numPlots, chosenData, spatialNormalization, width, height, 1 );
  end
  
  numCellsAuto = 0;
  
  % determine center of data points
  centerPosAuto = zeros(1, 3);
  for i=1:size(trackingMatrix,1)
    timeStep = trackingMatrix( i, 8 );
    if timeStep == t
      cen = trackingMatrix( i, 3:5 );
      centerPosAuto = centerPosAuto + cen;
      numCellsAuto = numCellsAuto+1;
    end
  end
  if numCellsAuto ~= 0
    centerPosAuto = centerPosAuto./numCellsAuto;
  end
  
  % store all auto positions
  cellCounter = 1;
  X = zeros(numCellsAuto, 2);

  % automatic segmentation result
  if renderNuclei == 1
    if chosenData < 6
      ax1 = subplot( numPlots/2, numPlots/2, 1 );
    else
      ax1 = subplot( 1, numPlots, 1 );
    end
    hold on
  end
  for i=1:size(trackingMatrix,1)
    timeStep = trackingMatrix( i, 8 );
    if timeStep == t
      lineage = trackingMatrix( i, 10 );
      if isKey( lineageAutoMap, lineage ) == 1
        lineageAutoMap( lineage ) = lineageAutoMap( lineage ) + 1;
      else
        lineageAutoMap( lineage ) = 1;
      end
      cen = trackingMatrix( i, 3:5 );
      if spatialNormalization == 1
        cen = cen - centerPosAuto;
      else
        cen = [ cen(1) -cen(2)+height cen(3) ];
      end
      
      % store the lineage information
      nodeID = trackingMatrix( i, 1 );
      prevNodeID = trackingMatrix( i, 7 );
      if prevNodeID == -1
        autoNodes{lineage}(cellAutoNumLin( lineage, 1 ) ) = 0;
        autoNodeIdToNumCellMap(nodeID) = cellAutoNumLin( lineage, 1 );
      else
        prevID = autoNodeIdToNumCellMap(prevNodeID);
        autoNodes{lineage}(cellAutoNumLin( lineage, 1 ) ) = prevID;
        autoNodeIdToNumCellMap(nodeID) = cellAutoNumLin( lineage, 1 );
      end
      cellAutoNumLin( lineage, 1 ) = cellAutoNumLin( lineage, 1 ) + 1;
      
      X(cellCounter, :) = [ cen(1) cen(2) ];
      cellCounter = cellCounter + 1;
      
      if renderNuclei == 1
        color = cmap( mod( lineage, numColors )+1, : );
        [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
          drawEllipse3d( cen(1), cen(2), cen(3), radEllip, radEllip, 0, 0 );
        set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
        set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
        nucleiCounter = nucleiCounter+1;
      end
    end
  end
  
  if renderNuclei == 1 && showAdditionalImages == 1
    tit = strcat( 'AutoTracking\_', 'T', {' '}, num2str(t), ',',...
      {' '}, '#Cells', {' '}, num2str(numCellsAuto), ',',...
      {' '}, '#Lineages', {' '}, num2str(lineageAutoMap.Count) );
    title( char(tit) );
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Segmentation result of raw data %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chosenData < 6
      ax2 = subplot( numPlots/2, numPlots/2, 2 );
    else
      ax2 = subplot( 1, numPlots, 2 );
    end
    hold on
    if t < 10
      digit = '00';
    elseif t < 100
      digit = '0';
    else
      digit = '';
    end
    csvPath = strcat( 'I:/SegmentationResults/Matlab/Segmentation/',...
      rawDataStr( 1, chosenData ), '/', rawDataStr( 1, chosenData ),...
      '_T', digit, num2str(t), '.dat' );
    clusRes = csvread( char(csvPath) );
    
    for c=1:size(clusRes, 1)
      color = cmap( 10, : );
      [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
        drawEllipse3d( clusRes(c, 1), clusRes(c, 2), 0, radEllip, radEllip, 0, 0 );
      set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
      set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
      nucleiCounter = nucleiCounter+1;
    end
    
    tit = strcat( 'segResult\_', 'T', {' '}, num2str(t), ', #Cells', num2str(size(clusRes, 1)) );
    title( char(tit) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MIP image of raw data %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if chosenData < 6
      ax3 = subplot( numPlots/2, numPlots/2, 3 );
    else
      ax3 = subplot( 1, numPlots, 3 );
    end
    hold on
    h1 = showMIP( rawDataStr( 1, chosenData ), 0, t );
    tit = strcat( 'rawMIP\_', 'T', {' '}, num2str(t) );
    title( char(tit) );
    
    if chosenData < 6
      % manual segmentation result
      ax4 = subplot( numPlots/2, numPlots/2, 4 );
      hold on
    end
  end

  if chosenData < 6
    % store all manual positions
    numCellsManual = numCellsPerTimeStep(t, 1);
    cellCounter = 1;
    Y = zeros(numCellsManual, 2);
    manLinCounter = 1;
    
    for j=1:dimData
      if cellData{j, 5} == t
        lin = cellData{j, 6};
        if isKey( lineageManualMap, lin ) == 1
          lineageManualMap( lin ) = lineageManualMap( lin ) + 1;
        else
          lineageManualMap( lin ) = 1;
          lineageManualColorMap( lin ) = manLinCounter;
          manLinCounter = manLinCounter + 1;
        end
        % get position of current cell
        p = [ cellData{j, 2} cellData{j, 3} cellData{j, 4} ];
        if spatialNormalization == 1
          p = p - centerPosPerTimeStep( t, : );
        end
        
        % store the lineage information
        if cellData{j, 9} == -1
          prevNodeID = -1;
        else
          prevNodeID = pairingFunction( cellData{j, 9}, t-1 );
        end
        nodeID = pairingFunction( cellData{j, 1}, t );
        
        if prevNodeID == -1
          manNodes{lin}(cellManNumLin( lin, 1 ) ) = 0;
        else
          prevLinID = manNodeIdToNumCellMap(prevNodeID);
          manNodes{lin}(cellManNumLin( lin, 1 ) ) = prevLinID;
        end
        
        manNodeIdToNumCellMap(nodeID) = cellManNumLin( lin, 1 );
        cellManNumLin( lin, 1 ) = cellManNumLin( lin, 1 ) + 1;
        
        Y(cellCounter, :) = [ p(1) p(2) ];
        cellCounter = cellCounter + 1;
        
        if renderNuclei == 1
          color = cmap( mod( lineageManualColorMap( lin ), numColors )+1, : );
          [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
            drawEllipse3d( p(1), p(2), p(3), radEllip, radEllip, 0, 0 );
          set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
          set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
          nucleiCounter = nucleiCounter+1;
        end
      end
    end
    if renderNuclei == 1
      tit = strcat( 'ManualSeg\_', 'T', {' '}, num2str(t), ',',...
        {' '}, '#Cells', {' '}, num2str(numCellsManual), ',',...
        {' '}, '#Lineages', {' '}, num2str(lineageManualMap.Count) );
      title( char(tit) );
    end
  end
  
  if renderNuclei == 1
    % image output options
    if t < 10
      digit = strcat( rawDataStr( 1, chosenData ), 'TimeStep', '_00' );
    elseif t < 100
      digit = strcat( rawDataStr( 1, chosenData ), 'TimeStep', '_0' );
    else
      digit = strcat( rawDataStr( 1, chosenData ), 'TimeStep', '_' );
    end
    
    filePath = strcat( imageDir, digit, num2str(t), '.png' );
    export_fig( gcf, char(filePath), '-m2', '-png' );
  end
end
%disp( strcat( 'Parameters:', num2str(paramThreshold), {' '}, num2str(paramTau) ) );

if renderAutoLineages == 1
  % tree plot figure settings
  f(2) = figure( 'Name', 'Tree Vis Auto', 'Position', [ 50 50 800 800 ] );
  for l=1:maxAutoLineages
    subplot( 1, maxAutoLineages, l );
    treeplot( autoNodes{l} );
  end
end

if renderManLineages == 1
  % tree plot figure settings
  f(3) = figure( 'Name', 'Tree Vis Man', 'Position', [ 50 50 800 800 ] );
  for l=1:maxManLineages
    subplot( 1, double(maxManLineages), double(l) );
    treeplot( manNodes{l} );
  end
end

% print elapsed time
toc