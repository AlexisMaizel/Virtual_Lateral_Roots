geomPath = strcat( pwd, '/geom3d' );
bezierPath = strcat( pwd, '/BezierPatchSurface' );
exportLibPath = strcat( pwd, '/export_fig' );
addpath( geomPath );
addpath( bezierPath );
addpath( exportLibPath );

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
startI = 1;
% endIndex
endI = 20;
% step size
deltaI = 1;
% min and max index
minI = 1;
maxI = 20;
% rendering stuff at curI and curI+deltaI
renderNuclei = 0;
renderPastNuclei = 0;
renderContour = 0;
% render bezier surface and control points
renderBezier = 1;
% export as image file
exportImages = 1;
% exclude nuclei outliers
excludeOutliers = 1;
% render only master file?
renderMasterFile = 1;
renderLineType = 2;
lineStr = { 'renderAllDeformationsInGrid'...
  'renderOnlyAveragedDeformationsInGrid' };
% line width of ellipses and semi axes
lineWidth = 1.2;
% enable z overlapping
overlapping = 1;
% use triangulation based on delaunay or alpha shape
% 1 -> delaunay
% 2 -> alpha shape
triangulationType = 2;
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
resGrid = 50;
% register data sets based on base instead of dome tip
registerBase = 0;
% apply the registration to all existing cell ranges and not only the one
% that all data sets share (cells in [18,143])
considerAllCells = 0;

% properties of bezier patches or surface
% use of exactly one bezier surface or a set of merged bezier patches
useBezierSurface = 1;
% dim of control points for the bezier surface -> dim x dim control points
dimCP = 7;
bezierSteps = 10;
% number of bezier patches (if selected) each with 15 control points
numPatches = 1;
bezierOffset = 0;
% add next to the height deformation based on the growth tensor information
% also an interpolated value of the bezier control points between the
% boundary control points (but only for the height)
interpolatedHeightGrowth = 1;
% properties of growth tensor
magnitudeScaling = 1.;
% manually increase the height of the upper control points
% such that the dome tip development is better captured in the model
emphasizeDomeTip = 1;
% represent the deformation only based on the longest deformation or all
% principal components of the deformation
onlyLongestDeformation = 0;

% figure properties
f = figure( 'Name', 'Mesh Deformation', 'Position', [0 0 1200 1600] );
% activate orbit rotation by default
cameratoolbar( 'SetMode', 'orbit' );
% activate none coord system by default for not resetting the camera up
% vector when rotating
cameratoolbar( 'SetCoordSys', 'none' );
% show camera toolbar by default
cameratoolbar( 'Show' );

% four plots for the current and next time step with bezier plots, the
% third one for the averaged deformations in the grid and the forth one for
% the visualization of the real 3D deformations
numPlots = 4;

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

if emphasizeDomeTip == 1
  totalMaxAxes(2) = totalMaxAxes(2)*2.;
end

% drawing handles
CONTOUR = [];
ELLIP = [];
ELLIPPATCH = [];
BEZIER = [];
BEZIERPOINTS = [];
BEZIERPOINTSPATCH = [];
L = [];
LP = [];
LN = [];
SD = [];
SDP = [];
SDN = [];
PC1 = [];
PC2 = [];
PC3 = [];

% path to image output
imageDir = strcat( 'images/' );
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

for pl=1:numPlots
  hold on;
  subplot( numPlots, 1, pl );
  % initialize 2D grid
  [ rows, columns ] =...
    generate2DGrid(...
    [totalMinAxes(1) totalMinAxes(2)],...
    [totalMaxAxes(1) totalMaxAxes(2)], resGrid );
end

% start measuring of elapsed time
tic

% bezier surface properties
initialMinPos = [ -170 -30 0 ];
initialMaxPos = [ 175 0 0];
lastMinPos = [ -280 -50 0 ];
lastMaxPos = [ 300 65 0];
% initialization of bezier surface
if useBezierSurface == 0
  [ curS, finalS, Q ] = initializeBezierPatches( numPatches, initialMinPos,...
    initialMaxPos, bezierOffset, emphasizeDomeTip );
else
  [ curS, finalS, Q ] = initializeBezierSurface( dimCP, initialMinPos,...
  initialMaxPos, bezierOffset, bezierSteps, emphasizeDomeTip );
end
% at start initialS = curS
initialS = curS;
surfaceDir = strcat( 'bezierGrowthSurfaces/' );
mkdir( char(surfaceDir) );

% export initial bezier surface with header file
fileName = strcat( surfaceDir, 'header.txt' );
fileId = fopen( char(fileName), 'w' );
% number of surfaces +1 because the intial surface is also included
fprintf( fileId, '%1d\n', endI-startI+2 );
% number of patches for each surface
fprintf( fileId, '%1d\n', numPatches );
% number of control points for each patch dimension
fprintf( fileId, '%1d\n', dimCP );
if useBezierSurface == 0
  exportBezierPatches( 0, initialS, surfaceDir );
else
  exportBezierSurface( 0, initialS, surfaceDir, dimCP );
end

if considerAllCells == 0
  regLastTimeStep = [ 269 277 230 344 213 ];
else
  % TODO
  regLastTimeStep = [];
end

% this variable is set to zero after the first traversal of the loop
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

% number of total cells per step
nucleiCounter = 1;
contourCounter = 1;
bezierCounter = 1;
bezierPatchCounter = 1;
lineDeformationCounter = 1;
pcCounter = 1;

% loop over all registered time steps
for curI=startI:deltaI:endI-1
  % check handles and if they exist hide them for redrawing
  hideHandle( CONTOUR );
  hideHandle( ELLIP );
  hideHandle( ELLIPPATCH );
  hideHandle( BEZIER );
  hideHandle( BEZIERPOINTS );
  hideHandle( BEZIERPOINTSPATCH );
  hideHandle( L );
  hideHandle( LP );
  hideHandle( LN );
  hideHandle( SD );
  hideHandle( SDP );
  hideHandle( SDN );
  hideHandle( PC1 );
  hideHandle( PC2 );
  hideHandle( PC3 );
  
  % initialize tile grid for storing start and end points of longest
  % deformation of a cell between two subsequent time steps which includes
  % the orientation and the magnitude of deformation
  
  % tile grid for direction vectors of deformations
  tileGrid = cell( rows*columns, 1 );
  
  curCellsC = zeros( numData, 1 );
  curCellsN = zeros( numData, 1 );
  allCellsC = zeros( numData, 1 );
  allCellsN = zeros( numData, 1 );
  
  %if begin == 1 && curI > 1
    % if curI is beginning at a later time step then update the bezier
    % surface depending on the current curI
    if useBezierSurface == 0
      [Q, curS] = updateBezierPatchesBoundary( curI/maxI, initialS, finalS, curS, numPatches );
      % export initial bezier surface
      exportBezierPatches( curI, curS, surfaceDir );
    else
      [Q, curS] = updateBezierSurfaceBoundary( curI/maxI, initialS, finalS, curS, dimCP, bezierSteps );
      % export initial bezier surface
      exportBezierSurface( curI, curS, surfaceDir, dimCP );
    end
  %end
  
  % min and max values for positions to initialize the bezier surface
  %minPos = [ 5000 5000 5000 ];
  %maxPos = [ -5000 -5000 -5000 ];
  
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
      
      if renderNuclei == 1
        if overlapping == 1
          zOffset = 0.5;
        else
          zOffset = 0.;
        end
        % render current nuclei in first plot
        for m=1:size(matPosC{dataIndex}, 1)
          subplot( numPlots, 1, 1 );
          hold on;
          % only render nuclei located in the master cell file
          if cellFileMap{dataIndex}( cellIdsC{dataIndex}(m) ) == 0
            p = matPosC{dataIndex}(m,:) - centerPosPerTimeStep{dataIndex}(curTC(dataIndex),:);
            p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
            [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
              drawEllipse3d( p(1), p(2), p(3)+dataIndex*zOffset+0.2, radEllip, radEllip, 0, 0 );
            set( ELLIP(nucleiCounter), 'color', colors( dataIndex, : ), 'LineWidth', lineWidth );
            set( ELLIPPATCH(nucleiCounter), 'FaceColor', colors( dataIndex, : ), 'FaceLighting', 'none' );
            nucleiCounter = nucleiCounter+1;
            
            % also draw as the past nuclei position if desired
            if renderPastNuclei == 1
              subplot( numPlots, 1, 2 );
              hold on;
              [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
                drawEllipse3d( p(1), p(2), p(3)+dataIndex*zOffset+0.2, radEllip, radEllip, 0, 0 );
              set( ELLIP(nucleiCounter), 'color', [ 1 1 1 ], 'LineWidth', lineWidth );
              set( ELLIPPATCH(nucleiCounter), 'FaceColor', [ 0 0 0 ], 'FaceAlpha', 0.3, 'FaceLighting', 'none' );
              nucleiCounter = nucleiCounter+1;
            end
          end
        end
        % and render nuclei at next time step in second plot
        subplot( numPlots, 1, 2 );
        hold on;
        for m=1:size(matPosN{dataIndex}, 1)
          % only render nuclei located in the master cell file
          if cellFileMap{dataIndex}( cellIdsN{dataIndex}(m) ) == 0
            p = matPosN{dataIndex}(m,:) - centerPosPerTimeStep{dataIndex}(curTN(dataIndex),:);
            p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
            [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
              drawEllipse3d( p(1), p(2), p(3)+dataIndex*zOffset+0.2, radEllip, radEllip, 0, 0 );
            set( ELLIP(nucleiCounter), 'color', colors( dataIndex, : ), 'LineWidth', lineWidth );
            set( ELLIPPATCH(nucleiCounter), 'FaceColor', colors( dataIndex, : ), 'FaceLighting', 'none' );
            nucleiCounter = nucleiCounter+1;
          end
        end
      end
      
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
      if dimC > 2 && renderContour == 1
        % render contour of current cells in first plot
        subplot( numPlots, 1, 1 );
        hold on;
        if triangulationType == 1
          K = convhull( cellsInMaster(:,1), cellsInMaster(:,2) );
          CONTOUR(contourCounter) = line(...
            cellsInMaster( K, 1 ), cellsInMaster( K, 2 ),...
            'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
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
      end
      cellsInMaster = [];
      for m=1:size(matPosN{dataIndex}, 1)
        if cellFileMap{dataIndex}( cellIdsN{dataIndex}(m) ) == 0
          p = matPosN{dataIndex}(m,:) - centerPosPerTimeStep{dataIndex}(curTN(dataIndex),:);
          p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
          cellsInMaster = [ cellsInMaster ; p ];
        end
      end
      dimN = size( cellsInMaster, 1 );
      curCellsN(dataIndex) = dimN;
      if dimN > 2 && renderContour == 1
        % and render contour of cells at next time step in second plot
        subplot( numPlots, 1, 2 );
        hold on;
        if triangulationType == 1
          K = convhull( cellsInMaster(:,1), cellsInMaster(:,2) );
          CONTOUR(contourCounter) = line(...
            cellsInMaster( K, 1 ), cellsInMaster( K, 2 ),...
            'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
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
        termTypeStr( 1, renderTermType ),...
        dataStr( 1, dataIndex ),...
        planePos, u, v, TF,...
        deltaT,...
        centerPosPerTimeStep{dataIndex},...
        curTC(dataIndex),...
        curTN(dataIndex),...
        renderMasterFile,...
        cellFileMap{dataIndex} );
      
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
      end
    else
      % if too few cells are available for triangulation then just continue
      continue;
    end
    
    % render the real 3D deformations
    subplot( numPlots, 1, 4 );
    hold on;
    % get index of all elongation types
    for l=1:dimL
      c = centerEllipse( l, : );
      for el=1:3
        index = minMaxEigenValueIndex( l, el );
        lineX = [ linePos( l, (index-1)*6 +1 ), linePos( l, (index-1)*6 +4 ) ];
        lineY = [ linePos( l, (index-1)*6 +2 ), linePos( l, (index-1)*6 +5 ) ];
        lineZ = [ 0.2, 0.2 ];
        color = [ lineColorIndex( l, (index-1)*3 +1 )...
          lineColorIndex( l, (index-1)*3 +2 )...
          lineColorIndex( l, (index-1)*3 +3 ) ];
        if el == 1
          PC1(pcCounter) = line( lineX, lineY, lineZ,...
            'Color', color, 'LineWidth', 1 );
        elseif el == 2
          PC2(pcCounter) = line( lineX, lineY, lineZ,...
            'Color', color, 'LineWidth', 1 );
        else
          PC3(pcCounter) = line( lineX, lineY, lineZ,...
            'Color', color, 'LineWidth', 1 );
        end
        pcCounter = pcCounter + 1;
      end
    end
    
  end
  
  % after all data is processed determine the average lines in the grid
  if strcmp( lineStr( 1, renderLineType ), 'renderOnlyAveragedDeformationsInGrid' )
    subplot( numPlots, 1, 3 );
    hold on;
    for gt=1:rows*columns
      numLines = size( tileGrid{ gt }, 1 );
      startPos = [ 0 0 ];
      endPos = [ 0 0 ];
      if numLines ~= 0
        for l=1:numLines
          startPos = startPos + tileGrid{ gt }(l, 1:2);
          endPos = endPos + tileGrid{ gt }(l, 3:4);
        end
        startPos = startPos./numLines;
        endPos = endPos./numLines;
        L(lineDeformationCounter) = drawGridArrow( startPos, endPos, gt, [totalMinAxes(1) totalMinAxes(2)],...
          [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 0 0 ] );
        %sd = determineSDDirections( tileGrid{ gt } );
        %         if sd ~= 0
        %           SD(lineDeformationCounter) = drawStandardDeviationArea( sd, averageSlope, gt, [totalMinAxes(1) totalMinAxes(2)],...
        %             [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 0 0 ] );
        %         end
      end
      lineDeformationCounter = lineDeformationCounter + 1;
    end
  % render all deformations as lines in the corresponding tiles
  elseif strcmp( lineStr( 1, renderLineType ), 'renderAllDeformationsInGrid' )
    subplot( numPlots, 1, 3 );
    hold on;
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
  
  % render the current bezier surface at curI
  if renderBezier == 1
    subplot( numPlots, 1, 1 );
    hold on;
    if useBezierSurface == 0
      for p=1:numPatches
        BEZIER(bezierPatchCounter) = mesh( Q(:,:,1,p), Q(:,:,2,p), Q(:,:,3,p), Q(:,:,3,p) );
        bezierPatchCounter = bezierPatchCounter + 1;
      end
      for p=1:numPatches
        for ci=1:4
          for cj=1:4
            [ BEZIERPOINTS(bezierCounter), BEZIERPOINTSPATCH(bezierCounter) ] =...
              drawEllipse3d( curS( ci, cj, 1, p ), curS( ci, cj, 2, p ), 0.2, radEllip, radEllip, 0, 0 );
            set( BEZIERPOINTS(bezierCounter), 'color', [ 0 0 0 ], 'LineWidth', lineWidth );
            set( BEZIERPOINTSPATCH(bezierCounter), 'FaceColor', [ 0 0 0 ], 'FaceLighting', 'none' );
            bezierCounter = bezierCounter + 1;
          end
        end
      end
    else
      BEZIER(bezierPatchCounter) = mesh( Q(:,:,1), Q(:,:,2), Q(:,:,3), Q(:,:,3) );
      bezierPatchCounter = bezierPatchCounter + 1;
      for ci=1:dimCP
        for cj=1:dimCP
          [ BEZIERPOINTS(bezierCounter), BEZIERPOINTSPATCH(bezierCounter) ] =...
            drawEllipse3d( curS( ci, cj, 1 ), curS( ci, cj, 2 ), 0.2, radEllip, radEllip, 0, 0 );
          set( BEZIERPOINTS(bezierCounter), 'color', [ 0 0 0 ], 'LineWidth', lineWidth );
          set( BEZIERPOINTSPATCH(bezierCounter), 'FaceColor', [ 0 0 0 ], 'FaceLighting', 'none' );
          bezierCounter = bezierCounter + 1;
        end
      end
    end
  end

  % update boundary control points of bezier surface for the next reg. time step
  if useBezierSurface == 0
    [Q, curS] = updateBezierPatchesBoundary( (curI+deltaI)/maxI, initialS, finalS, curS, numPatches );
    % update the inner control points of the bezier surface
    if strcmp( termTypeStr( 1, renderTermType ), 'B' )
      [Q, curS] = updateBezierPatches( tileGrid, curS, totalMinAxes, totalMaxAxes,...
        resGrid, rows, columns, magnitudeScaling, numPatches, interpolatedHeightGrowth );
    else
      [Q, curS] = updateBezierPatches( tileGrid, curS, totalMinAxes, totalMaxAxes,...
        resGrid, rows, columns, magnitudeScaling, numPatches, interpolatedHeightGrowth );
    end
    exportBezierPatches( curI+deltaI, curS, surfaceDir );
  else
    [Q, curS] = updateBezierSurfaceBoundary( (curI+deltaI)/maxI, initialS, finalS, curS, dimCP, bezierSteps );
    % update the inner control points of the bezier surface
    if strcmp( termTypeStr( 1, renderTermType ), 'B' )
      [Q, curS] = updateBezierSurface( tileGrid, curS, totalMinAxes, totalMaxAxes,...
        resGrid, rows, columns, magnitudeScaling, dimCP, bezierSteps, interpolatedHeightGrowth );
    else
      [Q, curS] = updateBezierSurface( tileGrid, curS, totalMinAxes, totalMaxAxes,...
        resGrid, rows, columns, magnitudeScaling, dimCP, bezierSteps, interpolatedHeightGrowth );
    end
    exportBezierSurface( curI+deltaI, curS, surfaceDir, dimCP );
  end
  
  % render the changing bezier surface at curI+deltaI
  if renderBezier == 1
    subplot( numPlots, 1, 2 );
    hold on;
    if useBezierSurface == 0
      for p=1:numPatches
        BEZIER(bezierPatchCounter) = mesh( Q(:,:,1,p), Q(:,:,2,p), Q(:,:,3,p), Q(:,:,3,p) );
        bezierPatchCounter = bezierPatchCounter + 1;
      end
      for p=1:numPatches
        for ci=1:4
          for cj=1:4
            [ BEZIERPOINTS(bezierCounter), BEZIERPOINTSPATCH(bezierCounter) ] =...
              drawEllipse3d( curS( ci, cj, 1, p ), curS( ci, cj, 2, p ), 0.2, radEllip, radEllip, 0, 0 );
            set( BEZIERPOINTS(bezierCounter), 'color', [ 0 0 0 ], 'LineWidth', lineWidth );
            set( BEZIERPOINTSPATCH(bezierCounter), 'FaceColor', [ 0 0 0 ], 'FaceLighting', 'none' );
            bezierCounter = bezierCounter + 1;
          end
        end
      end
    else
      BEZIER(bezierPatchCounter) = mesh( Q(:,:,1), Q(:,:,2), Q(:,:,3), Q(:,:,3) );
      bezierPatchCounter = bezierPatchCounter + 1;
      for ci=1:dimCP
        for cj=1:dimCP
          [ BEZIERPOINTS(bezierCounter), BEZIERPOINTSPATCH(bezierCounter) ] =...
            drawEllipse3d( curS( ci, cj, 1 ), curS( ci, cj, 2 ), 0.2, radEllip, radEllip, 0, 0 );
          set( BEZIERPOINTS(bezierCounter), 'color', [ 0 0 0 ], 'LineWidth', lineWidth );
          set( BEZIERPOINTSPATCH(bezierCounter), 'FaceColor', [ 0 0 0 ], 'FaceLighting', 'none' );
          bezierCounter = bezierCounter + 1;
        end
      end
    end
  end

  % render only first and last frame
  if curI == startI || curI == endI
    render = 1;
  else
    render = 0;
  end
  
  if exportImages == 1 || render == 1
    exportDeformationResult( totalMinAxes, totalMaxAxes, bezierOffset, numData, pureDataStr,...
      curTC, curTN, curCellsC, curCellsN, allCellsC, allCellsN, curI, deltaI,...
      viewStr( 1, cView ), imageDir, f, numPlots, renderTermType,...
      contribut, positiveEVVector, startData, endData );
  else
    % else just print the current registered step
    disp( strcat( {'Deformation between normalized step '}, num2str(curI),...
  {' and '} , num2str(curI+deltaI) ) );
  end
  % the first traversal is done
  begin = 0;
end

% print elapsed time
toc

