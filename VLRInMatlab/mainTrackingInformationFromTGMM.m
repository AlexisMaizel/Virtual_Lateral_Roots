% execute compileMex.m before to generate mex files

cd('C:\Jens\VLRRepository\VLRInMatlab')
GeomPath = strcat( pwd, '/geom3d' );
addpath( GeomPath );
ExportPath = strcat( pwd, '/export_fig' );
addpath( ExportPath );
TGMMPath = strcat( pwd, '/readTGMM_XMLoutput' );
addpath( TGMMPath );

% output format of values
format longG

% path to image output
imageDir = strcat( 'images/Segmentation/' );
mkdir( char(imageDir) );

%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
requireNewSegmentation = 1;
requireNewTracking = 1;
performParameterOptimization = 0;
paramThresStart = 500;
paramThresEnd = 500;
paramThresStep = 100;
paramTauStart = 150;
paramTauEnd = 150;
paramTauStep = 50;
cellDiffWeight = 100;
startT = 1;
endT = 10;
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };
rawDataStr = { '120830' '121204' '121211' '130508' '130607' };
chosenData = 3;
% normalize data points based on their center
spatialNormalization = 1;
resultPath = strcat( 'I:\SegmentationResults\TGMM\', rawDataStr( 1, chosenData ), '\' );
radEllip = 10;
lineWidth = 1;
numColors = 40;
%  amount of plots
numPlots = 2;
% colormap for ellipses
cmap = colorcube(numColors);
ELLIP = [];
ELLIPPATCH = [];
nucleiCounter = 1;
lineageAutoMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
lineageManualMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%

% start measuring elapsed time
tic

% read raw data (manual segmentation and tracking)
[ cellData, dimData, centerPosPerTimeStep, numCellsPerTimeStep ] =...
  readRawData( dataStr( 1, chosenData ) );

% instance storing the parameters for having the minimal errors
minErrorsPerTimeStep = ones( endT-startT+1, 1 ).*realmax('double');
minThresPerTimeStep = zeros( endT-startT+1, 1 );
minTausPerTimeStep = zeros( endT-startT+1, 1 );

oldThres = 0;
for thres=paramThresStart:paramThresStep:paramThresEnd
  for tau=paramTauStart:paramTauStep:paramTauEnd
    % execute TGMM segmentation and tracking
    cd('C:\Jens\TGMM_Supplementary_Software_1_0\build')
    tRange = strcat( num2str(startT), {' '}, num2str(endT) );
    dataSelect = strcat( 'TGMM_configFile', rawDataStr( 1, chosenData ), '.txt' );
    % first create TGMMconfigfile
    exportTGMMConfigFile( rawDataStr( 1, chosenData ), 2.0, thres, tau );
    if requireNewSegmentation == 1
      cmdSegmentation = strcat( 'nucleiChSvWshedPBC\Release\ProcessStackBatchMulticore.exe AutoTGMMConfig\', dataSelect, {' '}, tRange );
      % execute segementation only if the background threshold was changed
      if oldThres ~= thres
        system( char(cmdSegmentation), '-echo');
        disp( 'Segmentation done.' );
        oldThres = thres;
      end
    end
    if requireNewTracking == 1
      % execute tracking
      cmdTracking = strcat( 'Release\TGMM AutoTGMMConfig\', dataSelect, {' '}, tRange );
      system( char(cmdTracking), '-echo');
      disp( 'Tracking done.' );
      cd('C:\Jens\VLRRepository\VLRInMatlab')
    end
    
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
      highlightError( rawDataStr( 1, chosenData ), t, thres, tau );
      continue;
    end
    
    % figure settings
    if performParameterOptimization == 0
      f = figure( 'Name', 'Segmentation', 'Position', [ 50 50 800 800 ] );
    end
    
    for t=startT:endT
      % subplot settings
      if performParameterOptimization == 0
        clf(f)
        for p=1:2
          subplot( numPlots, 1, p );
          hold on
          if spatialNormalization == 0
            xMinMax = [ -50 750 ];
            if p == 1
              yMinMax = [ -450 -50 ];
            else
              yMinMax = [ 0 400 ];
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
      
      numCellsAuto = 0;
      numCellsManual = 0;
      
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
      centerPosAuto = centerPosAuto./numCellsAuto;
      
      % store all auto positions
      cellCounter = 1;
      X = zeros(numCellsAuto, 2);
      
      % automatic segmentation result
      if performParameterOptimization == 0
        subplot( numPlots, 1, 1 );
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
          end
          X(cellCounter, :) = [ cen(1) cen(2) ];
          cellCounter = cellCounter + 1;
          
          if performParameterOptimization == 0
            color = cmap( mod( lineage, numColors )+1, : );
            [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
              drawEllipse3d( cen(1), -cen(2), cen(3), radEllip, radEllip, 0, 0 );
            set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
            set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
            nucleiCounter = nucleiCounter+1;
          end
        end
      end
      
      if performParameterOptimization == 0
        tit = strcat( 'AutoSeg\_', 'T', {' '}, num2str(t), ',',...
          {' '}, '#Cells', {' '}, num2str(numCellsAuto), ',',...
          {' '}, '#Lineages', {' '}, num2str(lineageAutoMap.Count) );
        title( char(tit) );
        
        % manual segmentation result
        subplot( numPlots, 1, 2 );
        hold on
      end
      
      % store all manual positions
      numCellsManual = numCellsPerTimeStep(t, 1);
      cellCounter = 1;
      Y = zeros(numCellsManual, 2);
      
      for j=1:dimData
        if cellData{j, 5} == t
          lin = cellData{j, 6};
          if isKey( lineageManualMap, lin ) == 1
            lineageManualMap( lin ) = lineageManualMap( lin ) + 1;
          else
            lineageManualMap( lin ) = 1;
          end
          % get position of current cell
          p = [ cellData{j, 2} cellData{j, 3} cellData{j, 4} ];
          if spatialNormalization == 1
            p = p - centerPosPerTimeStep( t, : );
          end
          Y(cellCounter, :) = [ p(1) p(2) ];
          cellCounter = cellCounter + 1;
          
          if performParameterOptimization == 0
            color = cmap( mod( lin, numColors )+1, : );
            [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
              drawEllipse3d( p(1), p(2), p(3), radEllip, radEllip, 0, 0 );
            set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
            set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
            nucleiCounter = nucleiCounter+1;
          end
        end
      end
      
      if performParameterOptimization == 0
        tit = strcat( 'ManualSeg\_', 'T', {' '}, num2str(t), ',',...
          {' '}, '#Cells', {' '}, num2str(numCellsManual), ',',...
          {' '}, '#Lineages', {' '}, num2str(lineageManualMap.Count) );
        title( char(tit) );
        
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
      else
        % perform kNN search for determing similarity between auto and
        % manual segmentation
        [ IDX, D ] = knnsearch(X, Y, 'k',1,'distance','euclidean');
        % calculate mean squared error
        error = 0;
        for d=1:size(D,1)
          error = error + (D(d,1)*D(d,1));
        end
        % also add the difference in the number of cells with a weight
        cellDiff = abs(numCellsAuto-numCellsManual);
        error = error + cellDiffWeight*cellDiff*cellDiff;
        error = error/(size(D,1)+1);
        
        %disp( strcat( 'Error=', num2str(error), {' '}, 'at T=', num2str(t),...
        %  ' with', {' '}, num2str(numCellsAuto), '(A) and', {' '},...
        %  num2str(numCellsManual), '(M) cells.') );
        
        % determine current minimum settings for time step
        if error < minErrorsPerTimeStep( t, 1 )
          minErrorsPerTimeStep( t, 1 ) = error;
          minThresPerTimeStep( t, 1 ) = thres;
          minTausPerTimeStep( t, 1 ) = tau;
        end
        exportCurrentMinError( rawDataStr( 1, chosenData ), t,...
          minThresPerTimeStep( t, 1 ), minTausPerTimeStep( t, 1 ),...
          minErrorsPerTimeStep( t, 1 ), numCellsAuto, numCellsManual );
        exportErrorAnalyis( rawDataStr( 1, chosenData ), t, thres, tau,...
          error, numCellsAuto, numCellsManual );
      end
    end
    disp( strcat( 'Parameters:', num2str(thres), {' '}, num2str(tau) ) );
  end
end

% print elapsed time
toc
