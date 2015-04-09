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
% line width of ellipses and semi axes
lineWidth = 1.2;
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

% dim of control points for the bezier surface -> dim x dim control points
dimCP = 7;
bezierSteps = 10;
% number of bezier patches (if selected) each with 15 control points
numPatches = 1;
bezierOffset = 0;
% add next to the height deformation based on the growth tensor information
% also an interpolated value of the bezier control points between the
% boundary control points (but only for the height)
interpolatedHeightGrowth = 0;
% properties of growth tensor
magnitudeScaling = 5.;
% manually increase the height of the upper control points
% such that the dome tip development is better captured in the model
emphasizeDomeTip = 1;
% affect also the neighborhood tiles depending on the displacement
% information such that the surounding tiles may affect the control points
affectNeighborHood = 1;

% figure properties
f = figure( 'Name', 'Mesh Deformation', 'Position', [0 0 600 800] );
% activate orbit rotation by default
cameratoolbar( 'SetMode', 'orbit' );
% activate none coord system by default for not resetting the camera up
% vector when rotating
cameratoolbar( 'SetCoordSys', 'none' );
% show camera toolbar by default
cameratoolbar( 'Show' );

% four plots for the current and next time step with bezier plots, the
% third one for all displacement directions and the forth one for
% the averaged direction of displacements
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
  1, cView, excludeOutliers, 0 );

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

% start measuring elapsed time
tic

% bezier surface properties
initialMinPos = [ -170 -30 0 ];
initialMaxPos = [ 175 0 0];
lastMinPos = [ -280 -50 0 ];
lastMaxPos = [ 300 65 0];
% initialization of bezier surface
[ curS, finalS, Q ] = initializeBezierSurface( dimCP, initialMinPos,...
  initialMaxPos, bezierOffset, bezierSteps, emphasizeDomeTip );
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
exportBezierSurface( 0, initialS, surfaceDir, dimCP );

if considerAllCells == 0
  regLastTimeStep = [ 269 277 230 344 213 ];
else
  % TODO
  regLastTimeStep = [];
end

% the variables of two successive time steps are added a 'C' (Current)
% or 'N' (Next) at the end
% number of cells for the current and next time step
numCellsC = zeros( numData, 1 );
numCellsN = zeros( numData, 1 );

% matrix of positions for current and next time step
curPos = cell( numData );
nextPos = cell( numData );

% number of current time step for each data
curTC = zeros( numData, 1 );
curTN = zeros( numData, 1 );

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
  
  % initialize tile grid for storing start and end points of nuclei
  % displacement of a cell between two subsequent time steps which includes
  % the orientation and the magnitude
  
  % tile grid for direction vectors of displacements
  tileGrid = cell( rows*columns, 1 );
  
  curCellsC = zeros( numData, 1 );
  curCellsN = zeros( numData, 1 );
  allCellsC = zeros( numData, 1 );
  allCellsN = zeros( numData, 1 );
  
  % if curI is beginning at a later time step then update the bezier
  % surface depending on the current curI
  [Q, curS] = updateBezierSurfaceBoundary( curI/maxI, initialS, finalS, curS, dimCP, bezierSteps );
  % export initial bezier surface
  exportBezierSurface( curI, curS, surfaceDir, dimCP );
  
  % loop over all data sets
  for dataIndex=startData:endData
    % init transformation of data sets to ensure a specific view and
    % registration of data sets
    [ u, v, dir, planePos, TF ] =...
      initTransformations( dataStr( 1, dataIndex ), cView, registerBase );
    
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
    
    % determine the subsequent position that linked with each other
    [ curPos{dataIndex}, nextPos{dataIndex} ] =...
      determineSubsequentPositionsAndTracking( dimData(dataIndex),...
      curTC(dataIndex), curTN(dataIndex), cellDatas{dataIndex},...
      dataStr(1, dataIndex), centerPosPerTimeStep{dataIndex},...
      excludeOutliers, planePos, u, v, TF );
    
    if renderNuclei == 1
      zOffset = 0.5;
      % render current nuclei in first plot
      for m=1:size(curPos{dataIndex}, 1)
        subplot( numPlots, 1, 1 );
        hold on;
        p = curPos{dataIndex}(m,:);
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
      % and render nuclei at next time step in second plot
      subplot( numPlots, 1, 2 );
      hold on;
      for m=1:size(nextPos{dataIndex}, 1)
        p = nextPos{dataIndex}(m,:);
        [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
          drawEllipse3d( p(1), p(2), p(3)+dataIndex*zOffset+0.2, radEllip, radEllip, 0, 0 );
        set( ELLIP(nucleiCounter), 'color', colors( dataIndex, : ), 'LineWidth', lineWidth );
        set( ELLIPPATCH(nucleiCounter), 'FaceColor', colors( dataIndex, : ), 'FaceLighting', 'none' );
        nucleiCounter = nucleiCounter+1;
      end
    end
    
    % render contour lines of nuclei
    dimC = size(curPos{dataIndex}, 1);
    % store the number of unique positions as the number of cells in the
    % current time step because there may be at least one division in the
    % next time steps; however curPos and nextPos share the same dimension
    % to easily "read" the tracking information
    curCellsC(dataIndex) = size( unique(curPos{dataIndex}, 'rows'), 1 );
    curCellsN(dataIndex) = dimC;
    if dimC > 2 && renderContour == 1
      % render contour of current cells in first plot
      subplot( numPlots, 1, 1 );
      hold on;
      K = convhull( curPos{dataIndex}(:,1), curPos{dataIndex}(:,2) );
      CONTOUR(contourCounter) = line( curPos{dataIndex}(K, 1), curPos{dataIndex}(K, 2),...
        'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
      contourCounter = contourCounter + 1;
      % and render contour of cells at next time step in second plot
      subplot( numPlots, 1, 2 );
      hold on;
      K = convhull( nextPos{dataIndex}(:,1), nextPos{dataIndex}(:,2) );
      CONTOUR(contourCounter) = line( nextPos{dataIndex}(K, 1), nextPos{dataIndex}(K, 2),...
        'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
      contourCounter = contourCounter + 1;
    end
    
    % if the normalized steps are too accurate resulting in the same time
    % step for the current and next step then just continue with the next
    % data set and do not perform any computations for this steps and data
    if curTC(dataIndex) == curTN(dataIndex)
      continue;
    end
    
    % check if numNormCells is greater or equal to the maximal number of
    % total cells
    if numNormCellsC > numCellsPerTimeStep{dataIndex}(maxT(dataIndex),1)
      continue;
    end

    % fill the tileGrid with the displacement information
    for l=1:dimC
      c = curPos{dataIndex}(l,:);
      n = nextPos{dataIndex}(l,:);
      displacement = [ c(1) c(2) n(1) n(2) ];
      % determine tileIndex of current position and add it to the
      % tile grid for averaging it afterwards for all data sets
      tileIndex = getTileIndex( c, [totalMinAxes(1) totalMinAxes(2)],...
        [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns );
      tileGrid{tileIndex} = [ tileGrid{tileIndex}; displacement ];
      
      % also affect the neighbor tiles for each displacement in the center
      % of tileIndex multiplied by an attenuation factor
      if affectNeighborHood == 1
        aFactor = 0.5;
        for i=1:8
          nIndex = getNeighborTileIndex( i, tileIndex, columns );
          tileGrid{nIndex} = [ tileGrid{nIndex}; displacement*aFactor ];
        end
      end
      
      % render the individual displacements for each data set
      subplot( numPlots, 1, 3 );
      hold on;
      L(lineDeformationCounter) = drawGridArrow( displacement(1:2),...
        displacement(3:4), tileIndex, [totalMinAxes(1) totalMinAxes(2)],...
        [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, colors( dataIndex, : ) );
      lineDeformationCounter = lineDeformationCounter + 1;
    end
  end
  
  % after all data is processed determine the average displacement direction
  % in each tile of the grid
  subplot( numPlots, 1, 4 );
  hold on;
  for gt=1:rows*columns
    numLines = size( tileGrid{gt}, 1 );
    startPos = [ 0 0 ];
    endPos = [ 0 0 ];
    if numLines ~= 0
      for l=1:numLines
        startPos = startPos + tileGrid{gt}(l, 1:2);
        endPos = endPos + tileGrid{gt}(l, 3:4);
      end
      startPos = startPos./numLines;
      endPos = endPos./numLines;
      L(lineDeformationCounter) = drawGridArrow( startPos, endPos, gt,...
        [totalMinAxes(1) totalMinAxes(2)],...
        [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 0 0 ] );
      lineDeformationCounter = lineDeformationCounter + 1;
    end
  end
  
  % render the current bezier surface at curI
  if renderBezier == 1
    subplot( numPlots, 1, 1 );
    hold on;
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
  
  % update boundary control points of bezier surface for the next reg. time step
  [Q, curS] = updateBezierSurfaceBoundary( (curI+deltaI)/maxI, initialS, finalS, curS, dimCP, bezierSteps );
  % update the inner control points of the bezier surface
  [Q, curS] = updateBezierSurface( tileGrid, curS, totalMinAxes, totalMaxAxes,...
    resGrid, rows, columns, magnitudeScaling, dimCP, bezierSteps, interpolatedHeightGrowth );
  exportBezierSurface( curI+deltaI, curS, surfaceDir, dimCP );
  
  % render the changing bezier surface at curI+deltaI
  if renderBezier == 1
    subplot( numPlots, 1, 2 );
    hold on;
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
  
  % render only first and last frame
  if curI == startI || curI == endI
    render = 1;
  else
    render = 0;
  end
  
  if exportImages == 1 || render == 1
    exportDisplacementResult( totalMinAxes, totalMaxAxes, bezierOffset, numData, pureDataStr,...
      curTC, curTN, curCellsC, curCellsN, numCellsC, numCellsN, curI, deltaI,...
      viewStr( 1, cView ), imageDir, f, numPlots );
  else
    % else just print the current registered step
    disp( strcat( {'Deformation between normalized step '}, num2str(curI),...
      {' and '} , num2str(curI+deltaI) ) );
  end
end

% print elapsed time
toc

