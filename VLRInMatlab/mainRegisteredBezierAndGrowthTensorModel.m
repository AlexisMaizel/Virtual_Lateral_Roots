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
startI = 1;
% endIndex
endI = 100;
% step size
deltaI = 1;
% min and max index
minI = 1;
maxI = 100;
% draw point set as ellipses at curI
drawNuclei = 0;
% draw contour of cell nuclei at curI
drawContour = 1;
% export as image file
exportImages = 0;
% exclude nuclei outliers
excludeOutliers = 1;
% render only master file?
renderMasterFile = 1;
% render bezier surface and control points
renderBezier = 1;
% line width of ellipses and semi axes
lineWidth = 1.2;
% enable z overlapping
overlapping = 1;
% use triangulation based on delaunay or alpha shape
% 1 -> delaunay
% 2 -> alpha shape
triangulationType = 1;
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
resGrid = 20;
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
[ divisionProperties, cellDatas, dimData, maxT, numCellsPerTimeStep, centerPosPerTimeStep, totalMinAxes, totalMaxAxes, cellFileMap ] =...
  prepareData( dataStr, startData, endData, numData, 'Ellipses', renderMasterFile, cView, 0 );

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

% number of total cells per step
numTotalCells = 1;
links = 1;

% start measuring of elapsed time
tic

% number of cells in the last registered step
maxNumCellsAtLastStep = 0;

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

% loop over all registered time steps
for curI=startI:deltaI:endI
  
  % check handles and if they exist hide them for redrawing
  hideHandle( CONTOUR );
  hideHandle( ELLIP );
  hideHandle( ELLIPPATCH );
  hideHandle( BEZIER );
  hideHandle( BEZIERPOINTS );
  hideHandle( BEZIERPOINTSPATCH );
  
  % store the number of total cells at reg step
  numTotalOfCurCells = 0;
  
  % counter for changing control points
  numCP = 1;
  
  allCurT = zeros( numData, 1 );
  allCells = zeros( numData, 1 );
  allCurCells = zeros( numData, 1 );
  
  % min and max values for positions to initialize the bezier surface
  minPos = [ 5000 5000 5000 ];
  maxPos = [ -5000 -5000 -5000 ];
  
  % loop over all data sets
  for dataIndex=startData:endData
    % init transformation of data sets to ensure a specific view and
    % registration of data sets
    [ u, v, dir, planePos, TF ] =...
      initTransformations( dataStr( 1, dataIndex ), cView, registerBase );
    
    % get the corresponding time steps for the registered step
    [ curT, numNormCells ] = getCorrespondingTimeStep( curI, minI, maxI,...
      maxT(dataIndex), numCellsPerTimeStep{dataIndex},...
      numCellsPerTimeStep{dataIndex}(maxT(dataIndex),1), considerAllCells );
    
    % store the chosen time step for each data set
    allCurT( dataIndex, 1 ) = curT;
    allCells( dataIndex, 1 ) = numCellsPerTimeStep{dataIndex}(curT,1);
    
    % check if numNormCells is greater or equal to the maximal number of
    % total cells
    if numNormCells > numCellsPerTimeStep{dataIndex}(maxT(dataIndex),1)
      continue;
    end
    
    [ nextT, numNextNormCells ]  = getCorrespondingTimeStep( curI+1, minI, maxI,...
      maxT(dataIndex), numCellsPerTimeStep{dataIndex},...
      numCellsPerTimeStep{dataIndex}(maxT(dataIndex),1), considerAllCells );
    
    % required vectors to store current positions and cell IDs
    curPos = [];
    cellIDs = [];
    
    for j=1:dimData( dataIndex )
      if cellDatas{ dataIndex }{ j, 5 } == curT
        if renderMasterFile == 1 &&...
            cellDatas{ dataIndex }{ j, 7 } ~= 0
          continue;
        end
        % this is a special case for the data set 130508 for which we
        % ignore the two cells that arise in the master cell file with
        % lineage ID 5
        if excludeOutliers == 1
          if strcmp( dataStr( 1, dataIndex ), '130508_raw' ) &&...
              cellDatas{ dataIndex }{ j, 6 } == 5
            continue;
          end
          
          %         if strcmp( dataStr( 1, dataIndex ), '130607_raw' ) &&...
          %             cellDatas{ dataIndex }{ j, 6 } == 6 &&...
          %             ( cellDatas{ dataIndex }{ j, 1 } == 179 ||...
          %             cellDatas{ dataIndex }{ j, 1 } == 118 )
          %           continue;
          %         end
          %
          %         if strcmp( dataStr( 1, dataIndex ), '121211_raw' ) &&...
          %             cellDatas{ dataIndex }{ j, 6 } == 13 &&...
          %             cellDatas{ dataIndex }{ j, 1 } == 413
          %           continue;
          %         end
        end
        
        % get position of current cell
        p = [ cellDatas{ dataIndex }{ j, 2 }...
          cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
        p = p - centerPosPerTimeStep{dataIndex}(regLastTimeStep(dataIndex),:);
        p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
        curPos = [ curPos ; p ];
        numTotalOfCurCells = numTotalOfCurCells +1;
        
        % update min and max values of positions
        [minPos, maxPos] = updateMinMax( p, minPos, maxPos );
        
        if drawNuclei == 1
          hold on
          if overlapping == 1
            zOffset = 0.5;
          else
            zOffset = 0.;
          end
          [ ELLIP(numTotalCells), ELLIPPATCH(numTotalCells) ] =...
            drawEllipse3d( p(1), p(2), p(3)+dataIndex*zOffset+0.2, radEllip, radEllip, 0, 0 );
          set( ELLIP(numTotalCells), 'color', colors( dataIndex, : ), 'LineWidth', lineWidth );
          set( ELLIPPATCH(numTotalCells), 'FaceColor', colors( dataIndex, : ), 'FaceLighting', 'none' );
          numTotalCells = numTotalCells+1;
        end
        
        cellIDs = [ cellIDs ; cellDatas{ dataIndex }{ j, 1 } ];
      end
    end
    % not that elegenat but necessary; go over the whole loop again to
    % match successor cells based on lineage information
    nextPos = [];
    daughterPos = [];
    if curT ~= nextT
      for j=1:dimData( dataIndex )
        % also check next time step
        if cellDatas{ dataIndex }{ j, 5 } == nextT
          if renderMasterFile == 1 &&...
              cellDatas{ dataIndex }{ j, 7 } ~= 0
            continue;
          end
          % this is a special case for the data set 130508 for which we
          % ignore the two cells that arise in the master cell file with
          % lineage ID 5
          if excludeOutliers == 1
            if strcmp( dataStr( 1, dataIndex ), '130508_raw' ) &&...
                cellDatas{ dataIndex }{ j, 6 } == 5
              continue;
            end
            
            %         if strcmp( dataStr( 1, dataIndex ), '130607_raw' ) &&...
            %             cellDatas{ dataIndex }{ j, 6 } == 6 &&...
            %             ( cellDatas{ dataIndex }{ j, 1 } == 179 ||...
            %             cellDatas{ dataIndex }{ j, 1 } == 118 )
            %           continue;
            %         end
            %
            %         if strcmp( dataStr( 1, dataIndex ), '121211_raw' ) &&...
            %             cellDatas{ dataIndex }{ j, 6 } == 13 &&...
            %             cellDatas{ dataIndex }{ j, 1 } == 413
            %           continue;
            %         end
          end
          
          % get position of next cell
          p = [ cellDatas{ dataIndex }{ j, 2 }...
            cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
          p = p - centerPosPerTimeStep{dataIndex}(regLastTimeStep(dataIndex),:);
          p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
          
          % ID of precursor
          cellID = cellDatas{ dataIndex }{ j, 1 };
          %precurID = cellDatas{ dataIndex }{ j, 8 };
          %T = cellDatas{ dataIndex }{ j, 5 }
          [ precurIDList, numEntries ] = getPrecursorIDList( cellDatas{ dataIndex }{ j, 8 } );
          
          % if the cell is a daughter cell of a non-division cell
          % then the pos of the next time step is located at the same index
          % position as in the previous time step
          
          [row,col] = find( cellIDs == cellID );
          if size( row, 1 ) > 0
            nextPos( row, :) = p;
            % else there is a division and we store both positions of
            % the two daughter cells to later generate a cell that is associated with
            % the cell at the previous time step
          else
            pre = numEntries;
            while pre > 0
              precurID = precurIDList( pre, 1 );
              [row,col] = find( cellIDs == precurID );
              if size( row, 1 ) > 0
                daughterPos = [ daughterPos; row p ];
                break;
              else
                pre = pre - 1;
              end
            end
          end
        end
      end
    else
      nextPos = curPos;
    end
    
    % store the total number of cells
    numCells = size( curPos, 1 );
    allCurCells( dataIndex, 1 ) = numCells;
    
    if drawContour == 1
      if numCells > 2
        if triangulationType == 1
          K = convhull( curPos( :, 1 ), curPos( :, 2 ) );
          CONTOUR(curI + maxI * dataIndex) =...
            line( curPos( K, 1 ), curPos( K, 2 ), curPos( K, 3 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
        else
          [VolC,ShapeC] = alphavol( [ curPos( :, 1 ), curPos( :, 2 ) ], sqrt( alphaRadiiVector( curT, 1 )) );
          K = ShapeC.bnd(:,1);
          dimK = size( K, 1 );
          if dimK > 1
            K(dimK+1,:) = K(1,:);
            CONTOUR(curI + maxI * dataIndex) = line( curPos( K, 1 ), curPos( K, 2 ), curPos( K, 3 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
          end
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
  allCurT, allCurCells, allCells, curI, numNormCells,...
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

