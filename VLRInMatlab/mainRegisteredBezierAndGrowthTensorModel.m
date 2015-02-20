geomPath = strcat( pwd, '/geom3d' );
bezierPath = strcat( pwd, '/BezierPatchSurface' );
addpath( geomPath );
addpath( bezierPath );

setenv('LC_ALL','C')

%%%%% setting of properties %%%%%%
% color range for different data sets
% { 'm' 'y' 'r' 'g' 'b' 'c' 'dr' };
colors = [ [ 1 0 1 ]; [ 1 1 0 ]; [ 1 0 0 ]; [ 0 1 0 ]; [ 0 0 1 ]; [ 0 1 1 ]; [ 0.5 0 0 ] ];
% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
cView = 2;
% startIndex
startI = 10;
% endIndex
endI = 10;
% step size
deltaI = 1;
% min and max index
minI = 1;
maxI = 10;
% draw point set as ellipses at curI+deltaI
renderNuclei = 1;
% draw contour of cell nuclei at curI+deltaI
renderContour = 0;
% export as image file
exportImages = 1;
% exclude nuclei outliers
excludeOutliers = 1;
% render only master file?
renderMasterFile = 1;
% render bezier surface and control points
renderBezier = 1;
renderAllLines = 0;
% line width of ellipses and semi axes
lineWidth = 1.2;
% enable z overlapping
overlapping = 1;
% use triangulation based on delaunay or alpha shape
% 1 -> delaunay
% 2 -> alpha shape
triangulationType = 1;
% term that should be included in the time evolution
renderTermType = 3;
termTypeStr = { 'B' 'T' 'All' };
% vector of data strings
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' '131203_raw' };
pureDataStr = { '120830' '121204' '121211' '130508' '130607' };
% vector of view strings
viewStr = { 'Top' 'Side' 'Radial' '3D' };
% number of subdivisions for ellipsoids
nEllip = 6;
% radius of ellipses
radEllip = 3;
% data Index:
% 1 -> 120830
% 2 -> 121204
% 3 -> 121211
% 4 -> 130508
% 5 -> 130607
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

% properties of bezier patches
numPatches = 4;
bezierOffset = 0;

% figure properties
f = figure( 'Name', 'Mesh Deformation', 'Position', [100 100 1600 1200] );
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
ELLIP = [];
ELLIPPATCH = [];
BEZIER = [];
BEZIERPOINTS = [];
BEZIERPOINTSPATCH = [];

% path to image output
imageDir = strcat( 'images/AllData/' );
mkdir( char(imageDir) );

% output format of values
format longG

% set colormap and colorbar depending on the number of cell files
cm = hsv( 6 );
colormap( cm );

% gca is the current axes handle
set( gca,'nextplot','replacechildren' );
% gcf is the current figure handle
%lighting phong
set( gcf, 'Renderer', 'zbuffer' );
lighting gouraud
set( gcf, 'Renderer', 'OpenGL' );
set( gcf,'nextplot','replacechildren' );
set( gcf, 'color', [ 1 1 1 ] );

% initialize 2D grid
[ rows, columns ] =...
  generate2DGrid(...
  [totalMinAxes(1) totalMinAxes(2)],...
  [totalMaxAxes(1) totalMaxAxes(2)], resGrid );

% line render index for each normalized step
lineRenderIndex = 1;

% start measuring of elapsed time
tic

% initialization of bezier surface
initialMinPos = [ -170 -30 0 ];
initialMaxPos = [ 175 0 0];
lastMinPos = [ -280 -50 0 ];
lastMaxPos = [ 300 65 0];
[ curS, finalS, Q ] = initializeBezierSurface( numPatches, initialMinPos,...
  initialMaxPos, bezierOffset );
% at start initialS = curS
initialS = curS;

% path to image output
surfaceDir = strcat( 'bezierGrowthSurfaces/' );
mkdir( char(surfaceDir) );

% export initial bezier surface with header file
fileName = strcat( surfaceDir, 'header.txt' );
fileId = fopen( char(fileName), 'w' );
% number of surfaces +1 because the intial surface is also included
fprintf( fileId, '%1d\n', endI-startI+2 );
% number of patches for each surface
fprintf( fileId, '%1d\n', numPatches );
exportBezierSurface( 0, curS, surfaceDir );

if considerAllCells == 0
  regLastTimeStep = [ 269 277 230 344 213 ];
else
  % TODO
  regLastTimeStep = [];
end

% this variable is set to zero after the first traversal of the loop
% since it indicated the begin of the loop; this is required in order to
% set the information of the next time step to the infos of the current
% one in each time step which prevent redundant computations, but this is
% valid for all loop traversals except the first one
begin = 1;

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

% loop over all registered time steps
for curI=startI:deltaI:endI
  % check handles and if they exist hide them for redrawing
  hideHandle( CONTOUR );
  hideHandle( ELLIP );
  hideHandle( ELLIPPATCH );
  hideHandle( BEZIER );
  hideHandle( BEZIERPOINTS );
  hideHandle( BEZIERPOINTSPATCH );
  
  % TODO
  % initialize tile grid for storing start and end points of longest
  % deformation of a cell between two subsequent time steps which includes
  % the orientation and the magnitude of deformation
  
  % tile grid for T term when there is no positive/negative issue
  tileGrid = cell( rows*columns, 1 );
  % tile grid for only positive axes (for B term)
  tileGridP = cell( rows*columns, 1 );
  % tile grid for only negative axes (for B term)
  tileGridN = cell( rows*columns, 1 );
  % tile grid for storing the average values of B and T terms related to
  % their sum
  tileGridC = cell( rows*columns, 1 );
  % tile grid for storing the average values of magnitudes
  tileGridM = cell( rows*columns, 1 );
  
  % counter for changing control points
  numCP = 1;
  % number of total cells per step
  numTotalCells = 1;
  
  allCells = zeros( numData, 1 );
  allCurCells = zeros( numData, 1 );
  
  % min and max values for positions to initialize the bezier surface
  %minPos = [ 5000 5000 5000 ];
  %maxPos = [ -5000 -5000 -5000 ];
  
  % loop over all data sets
  for dataIndex=startData:endData
    % init transformation of data sets to ensure a specific view and
    % registration of data sets
    [ u, v, dir, planePos, TF ] =...
      initTransformations( dataStr( 1, dataIndex ), cView, registerBase );
    
    % get the corresponding time steps for the registered step
    if begin ~= 1
      curTC(dataIndex) = curTN(dataIndex);
      numNormCellsC = numNormCellsN;
      [ curTN(dataIndex), numNormCellsN ] = getCorrespondingTimeStep( curI+deltaI, minI, maxI,...
        maxT(dataIndex), numCellsPerTimeStep{dataIndex},...
        numCellsPerTimeStep{dataIndex}(maxT(dataIndex),1), considerAllCells );
      % initial curI at start
    else
      [ curTC(dataIndex), numNormCellsC ] = getCorrespondingTimeStep( curI, minI, maxI,...
        maxT(dataIndex), numCellsPerTimeStep{dataIndex},...
        numCellsPerTimeStep{dataIndex}(maxT(dataIndex),1), considerAllCells );
      
      [ curTN(dataIndex), numNormCellsN ] = getCorrespondingTimeStep( curI+deltaI, minI, maxI,...
        maxT(dataIndex), numCellsPerTimeStep{dataIndex},...
        numCellsPerTimeStep{dataIndex}(maxT(dataIndex),1), considerAllCells );
    end
    
    % update cell information
    if begin ~= 1
      % number of cells for current and next time step and data set
      numCellsC(dataIndex) = numCellsN(dataIndex);
      numCellsN(dataIndex) = numCellsPerTimeStep{dataIndex}(curTN(dataIndex), 1);
      matPosC{dataIndex} = matPosN{dataIndex};
      matPosN{dataIndex} = [];
      cellIdsC{dataIndex} = cellIdsN{dataIndex};
      cellIdsN{dataIndex} = [];
      
      % matrix of positions for current time step
      matPos2 = zeros( numCellsN(dataIndex), 3 );
      cellIds2 = zeros( numCellsN(dataIndex), 3 );
      cellPrecursors2 = cell( numCellsN(dataIndex),1 );
      
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
        if cellDatas{ dataIndex }{ j, 5 } == curTN(dataIndex)
          pos = [ cellDatas{ dataIndex }{ j, 2 }...
            cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
          matPos2(nc2, :) = pos;
          cellIds2(nc2, :) = cellDatas{ dataIndex }{ j, 1 };
          cellPrecursors2{nc2} = cellDatas{ dataIndex }{ j, 8 };
          nc2 = nc2 + 1;
        end
      end
    else
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
    end
    
    % if at least three cells exists
    if numCellsC(dataIndex) > 3 && numCellsN(dataIndex) > 3
      if begin ~= 1
        uniqueEdgesC{dataIndex} = uniqueEdgesN{dataIndex};
        triC{dataIndex} = triN{dataIndex};
        if triangulationType == 1
          % delaunay triangulation
          triN{dataIndex} = delaunayTriangulation( matPosN{dataIndex}(:,1), matPosN{dataIndex}(:,2), matPosN{dataIndex}(:,3) );
          uniqueEdgesN{dataIndex} = edges(triN{dataIndex});
        else
          % alpha shape triangulation
          [VolN,ShapeN] = alphavol( [ matPosN{dataIndex}(:,1) matPosN{dataIndex}(:,2) matPosN{dataIndex}(:,3) ], sqrt( alphaRadiiVector( curTN(dataIndex), 1 )) );
          triN{dataIndex} = ShapeN.tri;
          uniqueEdgesN{dataIndex} = getUniqueEdges( triN{dataIndex} );
        end
        numTotalLinksC(dataIndex) = numTotalLinksN(dataIndex);
        numTotalLinksN(dataIndex) = size( uniqueEdgesN{dataIndex}, 1 );
        numTotalAveragedLinks(dataIndex) = ( numTotalLinksC(dataIndex) + numTotalLinksN(dataIndex) )/2.;
      else
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
      end
      
      if renderNuclei == 1
        hold on
        if overlapping == 1
          zOffset = 0.5;
        else
          zOffset = 0.;
        end
        for m=1:size(matPosN{dataIndex}, 1)
          % only render nuclei located in the master cell file
          if cellFileMap{dataIndex}( cellIdsN{dataIndex}(m) ) == 0
            p = matPosN{dataIndex}(m,:) - centerPosPerTimeStep{dataIndex}(curTN(dataIndex),:);
            p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
            [ ELLIP(numTotalCells), ELLIPPATCH(numTotalCells) ] =...
              drawEllipse3d( p(1), p(2), p(3)+dataIndex*zOffset+0.2, radEllip, radEllip, 0, 0 );
            set( ELLIP(numTotalCells), 'color', colors( dataIndex, : ), 'LineWidth', lineWidth );
            set( ELLIPPATCH(numTotalCells), 'FaceColor', colors( dataIndex, : ), 'FaceLighting', 'none' );
            numTotalCells = numTotalCells+1;
          end
        end
      end
      
      % render contour lines of nuclei located in master cell file only
      if renderContour == 1
        hold on
        cellsInMaster = [];
        for m=1:size(matPosN{dataIndex}, 1)
          if cellFileMap{dataIndex}( cellIdsN{dataIndex}(m) ) == 0
            p = matPosN{dataIndex}(m,:) - centerPosPerTimeStep{dataIndex}(curTN(dataIndex),:);
            p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
            cellsInMaster = [ cellsInMaster ; p ];
          end
        end
        dimC = size( cellsInMaster, 1 );
        if dimC > 2
          K = convhull( cellsInMaster(:,1), cellsInMaster(:,2) );
          CONTOUR(curI + maxI * dataIndex) = line(...
            cellsInMaster( K, 1 ), cellsInMaster( K, 2 ),...
            'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
        end
      end
      
      % if the normalized steps are too accurate resulting in the same time
      % step for the current and next step then just continue with the next
      % data set and do not perform any computations for this steps and data
      if curTC(dataIndex) == curTN(dataIndex)
        continue;
      end
      
      allCells( dataIndex, 1 ) = numCellsPerTimeStep{dataIndex}(curTC(dataIndex),1);
      
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
        timePositions, indexColorSet, contributions, magnitudes ]...
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
        termTypeStr( 1, renderTermType ),...
        dataStr( 1, dataIndex ),...
        planePos, u, v, TF,...
        deltaT,...
        centerPosPerTimeStep{dataIndex},...
        curTC(dataIndex),...
        curTN(dataIndex),...
        renderMasterFile,...
        cellFileMap{dataIndex} );
      
      if overlapping == 1
        zOffset = 0.5;
      else
        zOffset = 0.;
      end
      
      dimL = size( centerEllipse, 1 );
      for l=1:dimL
        maxSemiPoint = [minMaxSemiAxisVector( l, 4 )...
          minMaxSemiAxisVector( l, 5 )...
          minMaxSemiAxisVector( l, 6 )];
        c = centerEllipse( l, : );
        
        % line coords of major semi axis
        lineX = [ maxSemiPoint(1), c(1) + c(1)-maxSemiPoint(1) ];
        lineY = [ maxSemiPoint(2), c(2) + c(2)-maxSemiPoint(2) ];
        
        % store the line positions in such a way that the position with
        % the smaller x value is always the first entry
        if lineX(1) <= lineX(2)
          lineDirection = [ lineX(1) lineY(1) lineX(2) lineY(2) ];
        else
          lineDirection = [ lineX(2) lineY(2) lineX(1) lineY(1) ];
        end
        
        % determine tileIndex of current ellipse position and add it to the
        % tile grid in order to average all deformations occurring in each
        % tile
        tileIndex = getTileIndex( c, [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns );
        if strcmp( termTypeStr( 1, renderTermType ), 'B' )
          if indexColorSet( l, 2 ) == 1
            tileGridP{ tileIndex } = [ tileGridP{ tileIndex }; lineDirection ];
          else
            tileGridN{ tileIndex } = [ tileGridN{ tileIndex }; lineDirection ];
          end
        else
          tileGrid{ tileIndex } = [ tileGrid{ tileIndex }; lineDirection ];
        end
      end
    else
      % if too few cells are available for triangulation then just continue
      continue;
    end
  end
  
  % after all data are processed determine the average visualization
  if renderAllLines == 0
    for gt=1:rows*columns
      if strcmp( termTypeStr( 1, renderTermType ), 'B' )
        % ignore empty tiles
        if size( tileGridP{ gt }, 1 ) ~= 0
          averageSlopeP = determineAverageSlope( tileGridP{ gt } );
          LP(lineRenderIndex) = drawAverageLines( averageSlopeP, gt, [totalMinAxes(1) totalMinAxes(2)],...
            [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 1 0 0 ], 0 );
          sd = determineSDDirections( tileGridP{ gt } );
          if sd ~= 0
            SDP(lineRenderIndex) = drawStandardDeviationArea( sd, averageSlopeP, gt, [totalMinAxes(1) totalMinAxes(2)],...
              [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 1 0 0 ] );
          end
        end
        if size( tileGridN{ gt }, 1 ) ~= 0
          averageSlopeN = determineAverageSlope( tileGridN{ gt } );
          LN(lineRenderIndex) = drawAverageLines( averageSlopeN, gt, [totalMinAxes(1) totalMinAxes(2)],...
            [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 0 1 ], 0 );
          sd = determineSDDirections( tileGridN{ gt } );
          if sd ~= 0
            SDN(lineRenderIndex) = drawStandardDeviationArea( sd, averageSlopeN, gt, [totalMinAxes(1) totalMinAxes(2)],...
              [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 0 1 ] );
          end
        end
      else
        % ignore empty tiles
        if size( tileGrid{ gt }, 1 ) ~= 0
          averageSlope = determineAverageSlope( tileGrid{ gt } );
          L(lineRenderIndex) = drawAverageLines( averageSlope, gt, [totalMinAxes(1) totalMinAxes(2)],...
            [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 0 0 ], 0 );
          sd = determineSDDirections( tileGrid{ gt } );
          if sd ~= 0
            SD(lineRenderIndex) = drawStandardDeviationArea( sd, averageSlope, gt, [totalMinAxes(1) totalMinAxes(2)],...
              [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 0 0 ] );
          end
        end
      end
      lineRenderIndex = lineRenderIndex + 1;
    end
  else
    for gt=1:rows*columns
      if strcmp( termTypeStr( 1, renderTermType ), 'B' )
        numLines = size( tileGridP{ gt }, 1 );
        for l=1:numLines
          % first compute the slope
          startPos = tileGridP{ gt }(l, 1:2);
          endPos = tileGridP{ gt }(l, 3:4);
          slope = (endPos(2)-startPos(2))/(endPos(1)-startPos(1));
          LP(lineRenderIndex) = drawAverageLines( slope, gt, [totalMinAxes(1) totalMinAxes(2)],...
            [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 1 0 0 ], 0 );
          lineRenderIndex = lineRenderIndex + 1;
        end
        numLines = size( tileGridN{ gt }, 1 );
        for l=1:numLines
          % first compute the slope
          startPos = tileGridN{ gt }(l, 1:2);
          endPos = tileGridN{ gt }(l, 3:4);
          slope = (endPos(2)-startPos(2))/(endPos(1)-startPos(1));
          LN(lineRenderIndex) = drawAverageLines( slope, gt, [totalMinAxes(1) totalMinAxes(2)],...
            [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 0 1 ], 0 );
          lineRenderIndex = lineRenderIndex + 1;
        end
      else
        numLines = size( tileGrid{ gt }, 1 );
        for l=1:numLines
          % first compute the slope
          startPos = tileGrid{ gt }(l, 1:2);
          endPos = tileGrid{ gt }(l, 3:4);
          slope = (endPos(2)-startPos(2))/(endPos(1)-startPos(1));
          L(lineRenderIndex) = drawAverageLines( slope, gt, [totalMinAxes(1) totalMinAxes(2)],...
            [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 0 0 ], 0 );
          lineRenderIndex = lineRenderIndex + 1;
        end
      end
    end
  end
  
  % update boundary control points of bezier surface
  [Q, curS] = updateBezierSurfaceBoundary( curI/maxI, initialS, finalS, curS, numPatches );
  
  % export bezier surface
  exportBezierSurface( curI, curS, surfaceDir );
  
  if renderBezier == 1
    hold on;
    for p=1:numPatches
      BEZIER(p + numPatches*curI) = mesh( Q(:,:,1,p), Q(:,:,2,p), Q(:,:,3,p) );
    end
    for p=1:numPatches
      cp = 1;
      for ci=1:4
        for cj=1:4
          [ BEZIERPOINTS(numCP), BEZIERPOINTSPATCH(numCP) ] =...
            drawEllipse3d( curS( ci, cj, 1, p ), curS( ci, cj, 2, p ), 0.2, radEllip, radEllip, 0, 0 );
          set( BEZIERPOINTS(numCP), 'color', [ 0 0 0 ], 'LineWidth', lineWidth );
          set( BEZIERPOINTSPATCH(numCP), 'FaceColor', [ 0 0 0 ], 'FaceLighting', 'none' );
          numCP = numCP + 1;
        end
      end
      cp = cp + 1;
    end
  end
  
  % render only first and last frame
  if curI == startI || curI == endI
    render = 1;
  else
    render = 0;
  end
  
  if exportImages == 1 || render == 1
    exportResult( totalMinAxes, totalMaxAxes, bezierOffset, numData, pureDataStr,...
      curTC, allCurCells, allCells, curI, numNormCells,...
      viewStr( 1, cView ), imageDir, f );
  else
    % else just print the current registered step
    if mod(curI,10) == 0
      disp( strcat( 'Normalized Step ', num2str(curI), '-Cells', num2str(numNormCells) ) );
    end
  end
  
end

% print elapsed time
toc

