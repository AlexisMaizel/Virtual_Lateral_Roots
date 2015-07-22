geomPath = strcat( pwd, '/geom3d' );
splinesPath = strcat( pwd, '/hobbysplines' );
bezierPath = strcat( pwd, '/BezierPatchSurface' );
exportLibPath = strcat( pwd, '/export_fig' );
exportSVGPath = strcat( pwd, '/plot2svg' );
addpath( geomPath );
addpath( splinesPath );
addpath( bezierPath );
addpath( exportLibPath );
addpath( exportSVGPath );

setenv('LC_ALL','C')

%%%%% setting of properties %%%%%%
% color range for different data sets
% { 'm' 'y' 'r' 'g' 'b' 'c' 'dr' };
colors = [ [ 1 0 1 ]; [ 1 1 0 ]; [ 1 0 0 ]; [ 0 1 0 ]; [ 0 0 1 ]; [ 0 1 1 ]; [ 0.5 0 0 ] ];
switchColor = 0; % switch color between white (0) or black (1)
% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
cView = 3;
% startIndex
startI = 19;
% endIndex
endI = 20;
% step size
deltaI = 1;
% min and max index
minI = 1;
maxI = 20;
renderSingleContours = 1;
renderNuclei = 0;
renderContour = 0;
renderMasterContour = 0;
contourEps = 20.;
% exclude nuclei outliers
excludeOutliers = 1;
% render only master file?
renderMasterFile = 1;
renderLineType = 2;
lineStr = { 'renderAllDeformationsInGrid'...
  'renderOnlyAveragedDeformationsInGrid' };
% line width of contour
lineWidth = 2;
% enable z overlapping
overlapping = 1;
% radius of ellipses
radEllip = 3;
% use triangulation based on delaunay or alpha shape
% 1 -> delaunay
% 2 -> alpha shape
triangulationType = 1;
% vector of data strings
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' '131203_raw' };
pureDataStr = { '120830' '121204' '121211' '130508' '130607' };
% vector of view strings
viewStr = { 'Top' 'Side' 'Radial' };
% data Index:
% 1 -> 120830 -> m
% 2 -> 121204 -> y
% 3 -> 121211 -> r
% 4 -> 130508 -> g
% 5 -> 130607 -> b
% 6 -> 131203
% start id of data
startData = 1;
% end id of data
endData = 5;
% num of data
numData = 5;
% resolution of grid
resGrid = 30;
% register data sets based on base instead of dome tip
registerBase = 0;
% apply the registration to all existing cell ranges and not only the one
% that all data sets share (cells in [18,143])
considerAllCells = 0;
% represent the deformation only based on the longest deformation or all
% principal components of the deformation
onlyLongestDeformation = 1;

% figure properties
f = figure( 'Name', 'Mesh Deformation', 'Position', [0 0 1600 1200] );
% activate orbit rotation by default
cameratoolbar( 'SetMode', 'orbit' );
% activate none coord system by default for not resetting the camera up
% vector when rotating
cameratoolbar( 'SetCoordSys', 'none' );
% show camera toolbar by default
cameratoolbar( 'Show' );

grid off;
xlabel('X');
ylabel('Y');
zlabel('Z');
camproj( 'orthographic' );

% apply preprocessing step of data
[ divisionProperties, cellDatas, dimData, maxT, numCellsPerTimeStep,...
  centerPosPerTimeStep, totalMinAxes, totalMaxAxes, cellFileMap ] =...
  prepareData( dataStr, startData, endData, numData, 'Ellipses',...
  renderMasterFile, cView, excludeOutliers, 0 );

% drawing handles
CONTOUR = [];
L = [];
SD = [];
ELLIP = [];
ELLIPPATCH = [];

% path to image output
imageDir = strcat( 'images/Deformation/' );
mkdir( char(imageDir) );

% output format of values
format longG

% set colormap and colorbar depending on the number of cell files
cm = jet( 256 );
colormap( cm );

% gca is the current axes handle
set( gca,'nextplot','replacechildren' );
% gcf is the current figure handle
%lighting phong
axis on
set( gcf, 'Renderer', 'zbuffer' );
lighting gouraud
shading interp
set( gcf, 'Renderer', 'OpenGL' );
set( gcf,'nextplot','replacechildren' );
if switchColor == 1
  set( gcf, 'color', [ 0 0 0 ] );
  set( gca, 'color', [ 0 0 0 ] );
  set( gca, 'XColor', [ 0.01 0.01 0.01 ] );
  set( gca, 'YColor', [ 0.01 0.01 0.01 ] );
else
  set( gcf, 'color', [ 1 1 1 ] );
  set( gca, 'color', [ 1 1 1 ] );
end
set( gca,'xtick',[],'ytick',[] )

% initialize 2D grid
[ rows, columns ] =...
  generate2DGrid(...
  [totalMinAxes(1) totalMinAxes(2)],...
  [totalMaxAxes(1) totalMaxAxes(2)], resGrid );

% start measuring of elapsed time
tic

% the variables of two successive time steps are added a 'C' (Current)
% or 'N' (Next) at the end
% number of cells for the current and next time step
numCellsC = zeros( numData, 1 );
numCellsN = zeros( numData, 1 );

% matrix of positions for current and next time step
matPosC = cell( numData );
matPosN = cell( numData );

% vector of object ids in order to access the cell
% file information in the cell file map
cellIdsC = cell( numData );
cellIdsN = cell( numData );

% vector of strings storing the precursors
cellPrecursorsN = cell( numData );

% delaunay triangulation
triC = cell( numData );
triN = cell( numData );

% number of current time step for each data
curTC = zeros( numData, 1 );
curTN = zeros( numData, 1 );

% unique edges in the triangulation
uniqueEdgesC = cell( numData );
uniqueEdgesN = cell( numData );

% number of total links in current and next time step
numTotalLinksC = zeros( numData, 1 );
numTotalLinksN = zeros( numData, 1 );
numTotalAveragedLinks = zeros( numData, 1 );

% number of total cells per step
contourCounter = 1;
lineDeformationCounter = 1;
nucleiCounter = 1;
averageDeltaT = 0;

if cView == 1
  globalMin = 0.2696;
  globalMax = 4.4663;
elseif cView == 2
  globalMin = 0.15823;
  globalMax = 4.042;
elseif cView == 3
  globalMin = 0.27102;
  globalMax = 2.7117;
end

% convert values to real units in log( µm / min )
globalMin = convertToReal( globalMin );
globalMax = convertToReal( globalMax );

% loop over all registered time steps
for curI=startI:deltaI:endI-1
  % check handles and if they exist hide them for redrawing
  hideHandle( CONTOUR );
  hideHandle( L );
  hideHandle( SD );
  hideHandle( ELLIP );
  hideHandle( ELLIPPATCH );
  
  % initialize tile grid for storing start and end points of longest
  % deformation of a cell between two subsequent time steps which includes
  % the orientation and the magnitude of deformation
  tileGrid = cell( rows*columns, 1 );
  tileGridM = cell( rows*columns, 1 );
  
  curCellsC = zeros( numData, 1 );
  curCellsN = zeros( numData, 1 );
  allCellsC = zeros( numData, 1 );
  allCellsN = zeros( numData, 1 );
  
  % store all cell positions for all data sets such that
  % the average contour is generated for all these cells
  allCells = [];
  allMasterCells = [];
  timeStepDiff = 0;
  
  % loop over all data sets
  for dataIndex=startData:endData
    % init transformation of data sets to ensure a specific view and
    % registration of data sets
    [ u, v, dir, planePos, TF ] =...
      initTransformations( dataStr( 1, dataIndex ), cView, registerBase );
    
    % get the alpha shape radii for all time steps
    if triangulationType == 2
      alphaRadiiVector = getAlphaRadius( dataStr( 1, dataIndex ) );
    end
    
    % get the corresponding time steps for the registered steps
    [ curTC(dataIndex), numNormCellsC ] = getCorrespondingTimeStep( curI, minI, maxI,...
      maxT(dataIndex), numCellsPerTimeStep{dataIndex},...
      numCellsPerTimeStep{dataIndex}(maxT(dataIndex),1), considerAllCells );
    
    [ curTN(dataIndex), numNormCellsN ] = getCorrespondingTimeStep( curI+deltaI, minI, maxI,...
      maxT(dataIndex), numCellsPerTimeStep{dataIndex},...
      numCellsPerTimeStep{dataIndex}(maxT(dataIndex),1), considerAllCells );
    
    timeStepDiff = timeStepDiff + curTN(dataIndex) - curTC(dataIndex);
    
    % update cell information
    % number of cells for current and next time step and data set
    numCellsC(dataIndex) = numCellsPerTimeStep{dataIndex}(curTC(dataIndex), 1);
    numCellsN(dataIndex) = numCellsPerTimeStep{dataIndex}(curTN(dataIndex), 1);
    
    % matrix of positions for current time step
    matPos1 = zeros( numCellsC(dataIndex), 3 );
    matPos2 = zeros( numCellsN(dataIndex), 3 );
    cellIds1 = zeros( numCellsC(dataIndex), 1 );
    cellIds2 = zeros( numCellsN(dataIndex), 1 );
    cellPrecursors2 = cell( numCellsN(dataIndex),1 );
    
    nc1 = 1;
    nc2 = 1;
    for j=1:dimData( dataIndex )
      % this is a special case for the data set 130508 for which we
      % ignore the two cells that arise in the master cell file with
      % lineage ID 5
      if excludeOutliers == 1
        if strcmp( dataStr( 1, dataIndex ), '130508_raw' ) &&...
            cellDatas{dataIndex}{j, 6} == 5
          continue;
        end
        %       if strcmp( dataStr( 1, dataIndex ), '130607_raw' ) &&...
        %           cellDatas{ dataIndex }{ j, 6 } == 6 &&...
        %           ( cellDatas{ dataIndex }{ j, 1 } == 179 ||...
        %           cellDatas{ dataIndex }{ j, 1 } == 118 )
        %         continue;
        %       end
        %
        %       if strcmp( dataStr( 1, dataIndex ), '121211_raw' ) &&...
        %           cellDatas{ dataIndex }{ j, 6 } == 13 &&...
        %           cellDatas{ dataIndex }{ j, 1 } == 413
        %         continue;
        %       end
      end
      if cellDatas{ dataIndex }{ j, 5 } == curTC(dataIndex)
        pos = [ cellDatas{ dataIndex }{ j, 2 }...
          cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
        matPos1(nc1, :) = pos;
        cellIds1(nc1, :) = cellDatas{ dataIndex }{ j, 1 };
        nc1 = nc1 + 1;
      end
      if cellDatas{ dataIndex }{ j, 5 } == curTN(dataIndex)
        pos = [ cellDatas{ dataIndex }{ j, 2 }...
          cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
        matPos2(nc2, :) = pos;
        cellIds2(nc2, :) = cellDatas{ dataIndex }{ j, 1 };
        cellPrecursors2{nc2} = cellDatas{ dataIndex }{ j, 8 };
        nc2 = nc2 + 1;
      end
    end
    
    matPosC{dataIndex} = matPos1;
    matPosN{dataIndex} = matPos2;
    cellIdsC{dataIndex} = cellIds1;
    cellIdsN{dataIndex} = cellIds2;
    cellPrecursorsN{dataIndex} = cellPrecursors2;
    
    % if at least three cells exists
    if numCellsC(dataIndex) > 3 && numCellsN(dataIndex) > 3
      if triangulationType == 1
        % delaunay triangulation
        triC{dataIndex} = delaunayTriangulation( matPosC{dataIndex}(:,1), matPosC{dataIndex}(:,2), matPosC{dataIndex}(:,3) );
        uniqueEdgesC{dataIndex} = edges( triC{dataIndex} );
        triN{dataIndex} = delaunayTriangulation( matPosN{dataIndex}(:,1), matPosN{dataIndex}(:,2), matPosN{dataIndex}(:,3) );
        uniqueEdgesN{dataIndex} = edges( triN{dataIndex} );
      else
        % alpha shape triangulation
        [VolC,ShapeC] = alphavol( [ matPosC{dataIndex}(:,1) matPosC{dataIndex}(:,2) matPosC{dataIndex}(:,3) ], sqrt( alphaRadiiVector( curTC(dataIndex), 1 )) );
        triC{dataIndex} = ShapeC.tri;
        uniqueEdgesC{dataIndex} = getUniqueEdges( triC{dataIndex} );
        [VolN,ShapeN] = alphavol( [ matPosN{dataIndex}(:,1) matPosN{dataIndex}(:,2) matPosN{dataIndex}(:,3) ], sqrt( alphaRadiiVector( curTN(dataIndex), 1 )) );
        triN{dataIndex} = ShapeN.tri;
        uniqueEdgesN{dataIndex} = getUniqueEdges( triN{dataIndex} );
      end
      numTotalLinksC(dataIndex) = size( uniqueEdgesC, 1 );
      numTotalLinksN(dataIndex) = size( uniqueEdgesN, 1 );
      numTotalAveragedLinks(dataIndex) = ( numTotalLinksC(dataIndex) + numTotalLinksN(dataIndex) )/2.;
      
      % render contour lines of nuclei located in master cell file only
      cellsInMaster = [];
      for m=1:size(matPosC{dataIndex}, 1)
        if cellFileMap{dataIndex}( cellIdsC{dataIndex}(m) ) == 0
          p = matPosC{dataIndex}(m,:) - centerPosPerTimeStep{dataIndex}(curTC(dataIndex),:);
          p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
          cellsInMaster = [ cellsInMaster ; p ];
        end
      end
      dimC = size( cellsInMaster, 1 );
      curCellsC(dataIndex) = dimC;
      cellsInMaster = [];
      for m=1:size(matPosN{dataIndex}, 1)
        if cellFileMap{dataIndex}( cellIdsN{dataIndex}(m) ) == 0
          p = matPosN{dataIndex}(m,:) - centerPosPerTimeStep{dataIndex}(curTN(dataIndex),:);
          p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
          cellsInMaster = [ cellsInMaster ; p ];
          allMasterCells = [ allMasterCells ; p ];
        end
        cp = matPosN{dataIndex}(m,:) - centerPosPerTimeStep{dataIndex}(curTN(dataIndex),:);
        cp = applyTransformations( cp, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
        allCells = [ allCells ; cp ];
      end
      dimN = size( cellsInMaster, 1 );
      curCellsN(dataIndex) = dimN;
      
      if renderSingleContours == 1
        %if dimN > 2
          if triangulationType == 1
            K = convhull( allCells(:,1), allCells(:,2) );
            
            n = size( K, 1 );
            points = cell( n, 1 );
            for l=1:n
              points{ l } = [ allCells( K(l), 1 ), allCells( K(l), 2 ) ];
            end
            hobbysplines( points, 'linestyle', {'linewidth', lineWidth},...
              'color', colors( dataIndex, : ) );
            
%             CONTOUR(contourCounter) = line(...
%               cellsInMaster( K, 1 ), cellsInMaster( K, 2 ),...
%               'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
          else
            [VolC,ShapeC] = alphavol( [ cellsInMaster(:, 1), cellsInMaster(:, 2) ],...
              sqrt( alphaRadiiVector( curTC(dataIndex), 1 )) );
            K = ShapeC.bnd(:,1);
            dimK = size( K, 1 );
            if dimK > 1
              K(dimK+1,:) = K(1,:);
              CONTOUR(contourCounter) = line( cellsInMaster( K, 1 ), cellsInMaster( K, 2 ),...
                cellsInMaster( K, 3 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
            end
          end
          contourCounter = contourCounter + 1;
        %end
      end
      
      if renderNuclei == 1
        for m=1:size(matPosN{dataIndex}, 1)
          % only render nuclei located in the master cell file
          %if cellFileMap{dataIndex}( cellIdsN{dataIndex}(m) ) == 0
            p = matPosN{dataIndex}(m,:) - centerPosPerTimeStep{dataIndex}(curTN(dataIndex),:);
            p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
            [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
              drawEllipse3d( p(1), p(2), p(3)+dataIndex*zOffset+0.2, radEllip, radEllip, 0, 0 );
            set( ELLIP(nucleiCounter), 'color', colors( dataIndex, : ), 'LineWidth', lineWidth );
            set( ELLIPPATCH(nucleiCounter), 'FaceColor', colors( dataIndex, : ), 'FaceLighting', 'none' );
            nucleiCounter = nucleiCounter+1;
          %end
        end
      end
      
      % if the normalized steps are too accurate resulting in the same time
      % step for the current and next step then just continue with the next
      % data set and do not perform any computations for this steps and data
      if curTC(dataIndex) == curTN(dataIndex)
        continue;
      end
      
      allCellsC( dataIndex, 1 ) = numCellsPerTimeStep{dataIndex}(curTC(dataIndex),1);
      allCellsN( dataIndex, 1 ) = numCellsPerTimeStep{dataIndex}(curTN(dataIndex),1);
      %allCellsC( dataIndex, 1 ) = numNormCellsC;
      %allCellsN( dataIndex, 1 ) = numNormCellsN;
      
      % check if numNormCells is greater or equal to the maximal number of
      % total cells
      if numNormCellsC > numCellsPerTimeStep{dataIndex}(maxT(dataIndex),1)
        continue;
      end
      
      % determine deltaT
      deltaT = curTN(dataIndex) - curTC(dataIndex);
      
      % compute time evolution for current deltaT and time step
      [ lineColorIndex, linePos, minMaxEigenValueIndex,...
        positiveEigenvalueVector, minMaxSemiAxisVector, centerEllipse,...
        timePositions, indexColorSet, contributions, magnitudes,...
        projectedStartEndCellPositions ]...
        = computeTimeEvolution(...
        uniqueEdgesC{dataIndex},...
        uniqueEdgesN{dataIndex},...
        cellIdsC{dataIndex},...
        cellIdsN{dataIndex},...
        numCellsN(dataIndex),...
        triC{dataIndex},...
        triN{dataIndex},...
        matPosC{dataIndex},...
        matPosN{dataIndex},...
        cellPrecursorsN{dataIndex},...
        triangulationType,...
        'All',...
        dataStr( 1, dataIndex ),...
        planePos, u, v, TF,...
        deltaT,...
        centerPosPerTimeStep{dataIndex},...
        curTC(dataIndex),...
        curTN(dataIndex),...
        renderMasterFile,...
        cellFileMap{dataIndex},...
        cView );
      
      positiveEVVector{dataIndex} = positiveEigenvalueVector;
      contribut{dataIndex} = contributions;
      
      if overlapping == 1
        zOffset = 0.5;
      else
        zOffset = 0.;
      end
      
      dimL = size( centerEllipse, 1 );
      for l=1:dimL
        c = centerEllipse( l, : );
        if onlyLongestDeformation == 1
          % get position of POINT on projected ellipse which has the
          % longest distance to the center of the ellipse
          maxSemiPoint = [minMaxSemiAxisVector( l, 4 )...
            minMaxSemiAxisVector( l, 5 )...
            minMaxSemiAxisVector( l, 6 )];
          
          % line coords of major semi axis
          lineX = [ maxSemiPoint(1), c(1) + c(1)-maxSemiPoint(1) ];
          lineY = [ maxSemiPoint(2), c(2) + c(2)-maxSemiPoint(2) ];
          
          % determine direction of deformation depending on the current and
          % the next position of the current cell; cellStartIndex is either 1
          % or 2 which denotes the start position of the deformation vector
          cellStartIndex = determineDeformationDirection(...
            projectedStartEndCellPositions(l,:), lineX, lineY );
          
          if cellStartIndex == 1
            lineDirection = [ lineX(1) lineY(1) lineX(2) lineY(2) ];
          else
            lineDirection = [ lineX(2) lineY(2) lineX(1) lineY(1) ];
          end
        else
          pcD = [];
          for pc=1:3
            % current principal component
            index = (pc-1)*6;
            PCVectorX = linePos(l, [ 1+index 4+index ] );
            PCVectorY = linePos(l, [ 2+index 5+index ] );
            % determine direction of each PC deformation depending on the current and
            % the next position of the current cell; cellStartIndex is either 1
            % or 2 which denotes the start position of the deformation vector
            cellStartIndex = determineDeformationDirection(...
              projectedStartEndCellPositions(l,:), PCVectorX, PCVectorY );
            
            if cellStartIndex == 1
              lDirT = [ linePos(l,1+index) linePos(l,2+index)...
                linePos(l,4+index) linePos(l,5+index) ];
            else
              lDirT = [ linePos(l,4+index) linePos(l,5+index)...
                linePos(l,1+index) linePos(l,2+index) ];
            end
            pcD = [ pcD ; lDirT ];
          end
          
          % average start and end points of the three principal components
          % of deformation in 3D projected onto the viewing plane
          %averagePCDir = [ (linePos(l,1:3)+linePos(l,7:9)+linePos(l,13:15))/3.,...
          %  (linePos(l,4:6)+linePos(l,10:12)+linePos(l,16:18))/3. ];
          averagePCDir = [ (pcD(1,1:2)+pcD(2,1:2)+pcD(3,1:2))/3.,...
            (pcD(1,3:4)+pcD(2,3:4)+pcD(3,3:4))/3. ];
          
          lineDirection = [ averagePCDir(1) averagePCDir(2)...
            averagePCDir(3) averagePCDir(4) ];
          
          % line coords of major semi axis
          %lineX = [ averagePCDir(1), c(1) + c(1)-averagePCDir(1) ];
          %lineY = [ averagePCDir(2), c(2) + c(2)-averagePCDir(2) ];
        end
        
        % determine tileIndex of current ellipse position and add it to the
        % tile grid in order to average all deformations occurring in each
        % tile
        tileIndex = getTileIndex( c, [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns );
        tileGrid{ tileIndex } = [ tileGrid{ tileIndex }; lineDirection ];
        lineLength = norm( [ lineDirection(3)-lineDirection(1) lineDirection(4)-lineDirection(2) ] );
        tileGridM{ tileIndex } = [ tileGridM{ tileIndex }; lineLength ];
      end
      
      disp( strcat( {'Between current time '}, num2str(curTC(dataIndex, 1)),...
        {' and next time '}, num2str(curTN(dataIndex, 1)),...
        {' for data '} , dataStr( 1, dataIndex ) ) );
      
    else
      % if too few cells are available for triangulation then just continue
      continue;
    end
  end
  
  % consider the average
  timeStepDiff = timeStepDiff/(endData-startData+1);
  averageDeltaT = averageDeltaT + timeStepDiff;
  
  % render average contour of all data sets
  if renderContour == 1 || renderMasterContour == 1
    if renderContour == 1
      cpos = allCells;
    else
      cpos = allMasterCells;
    end
    dimA = size( cpos, 1 );
    % increase the number of cells by adding four points to each nuclei
    % postion
    for c=1:dimA
      p1 = [ cpos(c, 1)+contourEps cpos(c, 2) 0. ];
      p2 = [ cpos(c, 1)-contourEps cpos(c, 2) 0. ];
      p3 = [ cpos(c, 1) cpos(c, 2)+contourEps 0. ];
      p4 = [ cpos(c, 1) cpos(c, 2)-contourEps 0. ];
      cpos = [ cpos ; p1 ; p2 ; p3 ; p4 ];
    end
    
    if dimA > 2
      if triangulationType == 1
        K = convhull( cpos(:,1), cpos(:,2) );
        n = size( K, 1 );
        points = cell( n, 1 );
        for l=1:n
          points{ l } = [ cpos( K(l), 1 ), cpos( K(l), 2 ) ];
        end
        if switchColor == 1
          co = [ 1 1 1 ];
        else
          co = [ 0 0 0 ];
        end
        hobbysplines( points, 'linestyle', {'linewidth', lineWidth},...
          'color', co );
%         CONTOUR(contourCounter) = line( cpos( K, 1 ), cpos( K, 2 ),...
%           'Color', [ 0 0 0 ], 'LineWidth', lineWidth );
      else
        [VolC,ShapeC] = alphavol( [ cpos(:, 1), cpos(:, 2) ],...
          sqrt( alphaRadiiVector( curTC(dataIndex), 1 )) );
        K = ShapeC.bnd(:,1);
        dimK = size( K, 1 );
        if dimK > 1
          K(dimK+1,:) = K(1,:);
          CONTOUR(contourCounter) = line( cpos( K, 1 ), cpos( K, 2 ),...
            cpos( K, 3 ), 'Color', [ 0 0 0 ], 'LineWidth', lineWidth );
        end
      end
      contourCounter = contourCounter + 1;
    end
  end
  
  % compute the average and min and max values of the magnitudes
%   minM = 10000;
%   maxM = 0;
%   for gt=1:rows*columns
%     numLines = size( tileGridM{ gt }, 1 );
%     if numLines ~= 0
%       average = 0;
%       for l=1:numLines
%         average = average + tileGridM{ gt }(l, 1);
%       end
%       average = average / numLines;
%       
%       if log(average) < minM
%         minM = log(average);
%       end
%       
%       if log(average) >= maxM
%         maxM = log(average);
%       end
%     end
%   end
  
%   % update global min value
%   if minM < globalMin
%     globalMin = minM;
%   end
%   
%   % update global max value
%   if maxM >= globalMax
%     globalMax = maxM;
%   end
  
  % after all data is processed determine the average lines in the grid
  if strcmp( lineStr( 1, renderLineType ), 'renderOnlyAveragedDeformationsInGrid' )
    for gt=1:rows*columns
      numLines = size( tileGrid{ gt }, 1 );
      if numLines ~= 0
        averageSlope = determineAverageSlope( tileGrid{ gt } );
        cValue = determineColorValue( tileGridM{ gt }, globalMin, globalMax );
        cValue = cValue * size( cm, 1 );
        cValue = round(cValue);
        if cValue == 0
          cValue = 1;
        end
        color = cm( cValue, : );
        L(lineDeformationCounter) = drawAverageLines( averageSlope, gt, [totalMinAxes(1) totalMinAxes(2)],...
          [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, color );
        sd = determineSDDirections( tileGrid{ gt } );
        if sd ~= 0
          SD(lineDeformationCounter) = drawStandardDeviationArea( sd, averageSlope, gt, [totalMinAxes(1) totalMinAxes(2)],...
            [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, color );
          lineDeformationCounter = lineDeformationCounter + 1;
        end
        lineDeformationCounter = lineDeformationCounter + 1;
      end
    end
  % render all deformations as lines in the corresponding tiles
  elseif strcmp( lineStr( 1, renderLineType ), 'renderAllDeformationsInGrid' )
    for gt=1:rows*columns
      numLines = size( tileGrid{ gt }, 1 );
      for l=1:numLines
        startPos = tileGrid{ gt }(l, 1:2);
        endPos = tileGrid{ gt }(l, 3:4);
        L(lineDeformationCounter) = drawGridArrow( startPos, endPos, gt, [totalMinAxes(1) totalMinAxes(2)],...
          [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 0 0 ] );
        lineDeformationCounter = lineDeformationCounter + 1;
      end
    end
  end
  
  % render only first and last frame
  if curI == startI || curI == endI
    render = 1;
  else
    render = 0;
  end
  
  hold off;
  set( f,'nextplot','replacechildren' );
  viewOffset = 100;
  axis( [ totalMinAxes(1) totalMaxAxes(1)...
    totalMinAxes(2) totalMaxAxes(2)...
    totalMinAxes(3)-viewOffset totalMaxAxes(3)+viewOffset ...
    totalMinAxes(3) totalMaxAxes(3) ] );
  %axis on
  daspect( [ 1 1 1 ] );
  cb = colorbar;
  set( gca, 'CLim', [ globalMin globalMax ] );
  set( cb, 'YTick', [ globalMin globalMax ] );
  set( cb, 'YTickLabel', { num2str(globalMin, '%.1f'), num2str(globalMax, '%.1f') } );
  xlabel(cb, 'log(µm/min)' );
  if switchColor == 1
    set( cb, 'color', [ 1 1 1 ] );
  else
    set( cb, 'color', [ 0 0 0 ] );
  end
  
%   % legend
%   hold on;
%   %leg = legend( '120830', '121204', '121211', '130508', '130607', '131203' );
%   leg = legend( '120830', '121204', '121211', '130508', '130607' );
%   set(leg, 'Location', 'NorthWestOutside');
%   linecolors = { [ 1 0 1 ], [ 0 0 0 ], [ 1 0 0 ], [ 0 1 0 ], [ 0 0 1 ] };
%   legendlinestyles( leg, {}, {}, linecolors );
  
  grid off;
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  title( strcat( 'Between Normalized Step ', num2str(curI), 'and',...
    num2str(curI+deltaI), '-Cells', num2str(numNormCellsN) ) );
  
  % image output options
  if curI < 10
    digit = strcat( viewStr( 1, cView ), '_00' );
  elseif curI < 100
    digit = strcat( viewStr( 1, cView ), '_0' );
  else
    digit = strcat( viewStr( 1, cView ), '_' );
  end
    
  % output with number of cells
  filePath = strcat( imageDir, digit, num2str(curI), '-Cells', num2str(numNormCellsN), '.svg' );
  %saveas( gcf, char(filePath) );
  %export_fig( gcf, char(filePath), '-m2', '-png' );
  axis off
  plot2svg( char(filePath), gcf, 'png' );
  
  % clear current figure window
  if curI ~= endI-1
    clf
    hold on;
    generate2DGrid( [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid );
    % hide axis ticks and labels
    if switchColor == 1
      set( gcf, 'color', [ 0 0 0 ] );
      set( gca, 'color', [ 0 0 0 ] );
      set( gca, 'XColor', [ 0.01 0.01 0.01 ] );
      set( gca, 'YColor', [ 0.01 0.01 0.01 ] );
    else
      set( gcf, 'color', [ 1 1 1 ] );
      set( gca, 'color', [ 1 1 1 ] );
    end
    set( gca,'xtick',[],'ytick',[] )
  end  
end

% average deltaT value for all data sets and reg time steps
%averageDeltaT = averageDeltaT/(endI-startI);
%disp( strcat( {'Average DeltaT '}, num2str(averageDeltaT) ) );
%disp( strcat( {'Global Min and Max '}, num2str(globalMin), ' and ', num2str(globalMax) ) );

% print elapsed time
toc