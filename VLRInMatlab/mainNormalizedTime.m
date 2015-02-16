%%%% The mesh deformation is based on the paper of Graner et al,
%%%% Discrete rearranging disordered patterns, part I: Robust
%%%% statistical tools in two or three dimensions, Eur. Phys. J.,
%%%% pp 349 -- 369, 2008

geomPath = strcat( pwd, '/geom3d' );
addpath( geomPath );

setenv('LC_ALL','C')

%%%%% setting of properties %%%%%%
% color range for different data sets
% { 'm' 'k' 'r' 'g' 'b' 'c' };
colors = [ [ 1 0 1 ]; [ 0 0 0 ]; [ 1 0 0 ]; [ 0 1 0 ]; [ 0 0 1 ] ];
% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
cView = 2;
% decide which term should be included in the time evolution
renderTermType = 3;
termTypeStr = { 'B' 'T' 'All' };
% render average lines or not
renderAverage = 1;
renderAveragePerTimeStep = 1;
renderAllLines = 0;
% take the average over the data set or over the time
averageOverData = 1;
% render division orientations
renderDivisions = 0;
% startIndex
startI = 19;
% endIndex
endI = 20;
% min and max index
minI = 1;
maxI = 20;
% deltaT value based on the paper mentioned above
% have to be a divider of the max index
deltaI = 1;
% draw delaunay tri?
drawDelaunay = 0;
% if set to one then only a single time step
% is rendered given by startT
exportType = 2;
% vector of data strings
exportTypeStr = { 'SingleFigure' 'AsImages' };
% render only master file?
renderMasterFile = 1;
% render principal components
renderPrincipalComponents = 0;
% render contributions of B and T terms
renderContributions = 0;
% line width of ellipses and semi axes
lineWidth = 1.2;
% enable z overlapping
overlapping = 1;
% only render specific ranges of cell numbers
% render contour (=convexHull) instead of ellipses
visualizationType = { 'Ellipsoids' 'Ellipses' };
visType = 2;
% offset for cell ranges
epsilon = 3;
lineRenderType = 4;
% vector of line render types
% 1. draw only those lines with the largest elongation of the ellipoids in 3D
% 2: render all three major lines of elongation of the ellipsoids in 3D
% 3. draw only those lines with the largest elongation of the ellipoids
%    projected onto the generated 2D plane
% 4: render all two major lines of elongation of the ellipses in 2D
lineStr = { 'renderLargest3DElongation'...
            'renderAll3DElongation'...
            'renderLargest2DElongation'...
            'renderAll2DElongation' };
% use triangulation based on delaunay or alpha shape
% 1 -> delaunay
% 2 -> alpha shape
triangulationType = 2;
% vector of data strings
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' '131203_raw' };
pureDataStr = { '120830', '121204', '121211', '130508', '130607', '131203' };
% vector of view strings
viewStr = { 'Top' 'Side' 'Radial' '3D' };
% number of subdivisions for ellipsoids
nEllip = 10;
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
  prepareData( dataStr, startData, endData, numData, visualizationType( 1, visType ), renderMasterFile, cView, 0 );

if strcmp( visualizationType( 1, visType ), 'Ellipsoids' )
  % surface instance
  S = [];
  % projected surface instance
  PS = [];
elseif strcmp( visualizationType( 1, visType ), 'Ellipses' )
  % semi axes instances
  MIN = [];
  MID = [];
  MAX = [];
  % ellipse instance
  ELLIP = [];
  ELLIPPATCH = [];
end

if renderAverage == 1
  % contour instance
  CONTOUR = [];
  % average line instance
  L = [];
  LP = [];
  LN = [];
  SD = [];
  DIV = [];
  SDD = [];
  SDP = [];
  SDN = [];
end

% render contributions for B and T term related to the sum of both
if renderContributions == 1
  BC = [];
  TC = [];
  BM = [];
  TM = [];
  COLORBAR1 = [];
  COLORBAR2 = [];
  COLORBAR3 = [];
end

if renderPrincipalComponents == 1
  % PC instance
  P = [];
end

% path to image output
imageDir = strcat( 'images/AllData/' );

% create directory if required
if strcmp( exportTypeStr( 1, exportType ), 'AsImages' )
  mkdir( char(imageDir) );
end

% output format of values
format longG

% gca is the current axes handle
set( gca,'nextplot','replacechildren' );
% gcf is the current figure handle
%lighting phong
set( gcf, 'Renderer', 'zbuffer' );
lighting gouraud
set( gcf, 'Renderer', 'OpenGL' );
set( gcf,'nextplot','replacechildren' );
set( gcf, 'color', [ 1 1 1 ] );

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

% storage for dividing cells plus division direction
cellDivisions = cell( numData, 1 );

cpuT = cputime;

hold on;
[ rows, columns ] = generate2DGrid( [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid );

% initialize grid size consisting of tiles
if averageOverData == 0
  tileGrid = cell( numData, rows*columns );
end

% line render index for each normalized step
lineRenderIndex = 1;

globalBTerm = zeros( maxI-minI, 1 );
globalTTerm = zeros( maxI-minI, 1 );

% loop over all normalized steps
curI=startI;
while curI < endI-deltaI+1
  % if averaging over all data sets then we initialize the tile grid 
  % for each normalized step
  if averageOverData == 1
    % tile grid for T term when there is no positive/negative issue
    tileGrid = cell( rows*columns, 1 );
    % tile grid for only positive axes (for B term)
    tileGridP = cell( rows*columns, 1 );
    % tile grid for only negative axes (for B term)
    tileGridN = cell( rows*columns, 1 );
    % tile grid for division directions
    tileGridD = cell( rows*columns, 1 );
    % tile grid for storing the average values of B and T terms related to
    % their sum
    tileGridC = cell( rows*columns, 1 );
    % tile grid for storing the average values of magnitudes
    tileGridM = cell( rows*columns, 1 );
  end
  
  if strcmp( visualizationType( 1, visType ), 'Ellipsoids' )
    if ishandle( S )
      set( S, 'Visible', 'off' );
    end
    if ishandle( PS )
      set( PS, 'Visible', 'off' );
    end
  elseif strcmp( visualizationType( 1, visType ), 'Ellipses' )
    if ishandle( MAX )
      set( MAX, 'Visible', 'off' );
    end
    if ishandle( MID )
      set( MID, 'Visible', 'off' );
    end
    if ishandle( MIN )
      set( MIN, 'Visible', 'off' );
    end
    if ishandle( ELLIP )
      set( ELLIP, 'Visible', 'off' );
    end
    if ishandle( ELLIPPATCH )
      set( ELLIPPATCH, 'Visible', 'off' );
    end
  end
  
  if renderAverage == 1
    if ishandle( L )
      set( L, 'Visible', 'off' );
    end
    if ishandle( LP )
      set( LP, 'Visible', 'off' );
    end
    if ishandle( LN )
      set( LN, 'Visible', 'off' );
    end
    if ishandle( DIV )
      set( DIV, 'Visible', 'off' );
    end
    if ishandle( SD )
      set( SD, 'Visible', 'off' );
    end
    if ishandle( SDD )
      set( SDD, 'Visible', 'off' );
    end
    if ishandle( SDP )
      set( SDP, 'Visible', 'off' );
    end
    if ishandle( SDN )
      set( SDN, 'Visible', 'off' );
    end
    if ishandle( CONTOUR )
      set( CONTOUR, 'Visible', 'off' );
    end
  end
  
  if renderContributions == 1
    if ishandle( BC )
      set( BC, 'Visible', 'off' );
    end
    if ishandle( TC )
      set( TC, 'Visible', 'off' );
    end
    if ishandle( BM )
      set( BM, 'Visible', 'off' );
    end
    if ishandle( TM )
      set( TM, 'Visible', 'off' );
    end
    if ishandle( COLORBAR1 )
      set( COLORBAR1, 'Visible', 'off' );
    end
    if ishandle( COLORBAR2 )
      set( COLORBAR2, 'Visible', 'off' );
    end
    if ishandle( COLORBAR3 )
      set( COLORBAR3, 'Visible', 'off' );
    end
  end
  
  if renderPrincipalComponents == 1
    if ishandle( P )
      set( P, 'Visible', 'off' );
    end
  end
  
  % min and max values of contributions for the color of the rectangles
  minContr = 10000.;
  maxContr = 0.;
  minMag = 10000.;
  maxMag = -10000.;
  
  % loop over all data sets
  for dataIndex=startData:endData
    % get stored eigenvectors for the last time step to set the same
    % direction view for each time step
    coeff = getNormalizedPrincipalComponents( dataStr( 1, dataIndex ), 1 );
    
    % set PC depending on the viewing direction
    if cView == 1
      dir = coeff(:,2);
      u = coeff(:,1);
      v = coeff(:,3);
    elseif cView == 2
      dir = coeff(:,3);
      u = coeff(:,1);
      v = coeff(:,2);
      if strcmp( dataStr( 1, dataIndex ), '121211_raw' )
        v = -v;
      end
    elseif cView == 3
      dir = -coeff(:,1);
      u = -coeff(:,3);
      v = coeff(:,2);
      if strcmp( dataStr( 1, dataIndex ), '121211_raw' )
        v = -v;
      end
    end
    
    % set plane position
    planePos = dir * 1;
    
    plane = [ planePos(1) planePos(2) planePos(3)...
      u(1) u(2) u(3)...
      v(1) v(2) v(3) ];
    TF = createBasisTransform3d( 'g', plane );
    
    % get the alpha shape radii for all time steps
    if triangulationType == 2
      alphaRadiiVector = getAlphaRadius( dataStr( 1, dataIndex ) );
    end
    
    % NEW: only continue with this time step that is synchronized in
    % number of cells for each data set
    if begin ~= 1
      numNormCellsN = getNormalizedCellNumber( curI+deltaI, 18, 143, minI, maxI );
      curTC(dataIndex) = curTN(dataIndex);
      curTN(dataIndex) = 1;
      for j=1:maxT(dataIndex)
        if numNormCellsN - epsilon < numCellsPerTimeStep{dataIndex}(j,1)...
            && numNormCellsN + epsilon > numCellsPerTimeStep{dataIndex}(j,1)
          curTN(dataIndex) = j;
          break;
        end
      end
    else
      numNormCellsC = getNormalizedCellNumber( curI, 18, 143, minI, maxI );
      numNormCellsN = getNormalizedCellNumber( curI+deltaI, 18, 143, minI, maxI );
      curTC(dataIndex) = 1;
      for j=1:maxT(dataIndex)
        if numNormCellsC - epsilon < numCellsPerTimeStep{dataIndex}(j,1)...
            && numNormCellsC + epsilon > numCellsPerTimeStep{dataIndex}(j,1)
          curTC(dataIndex) = j;
          break;
        end
      end
      curTN(dataIndex) = 1;
      for j=1:maxT(dataIndex)
        if numNormCellsN - epsilon < numCellsPerTimeStep{dataIndex}(j,1)...
            && numNormCellsN + epsilon > numCellsPerTimeStep{dataIndex}(j,1)
          curTN(dataIndex) = j;
          break;
        end
      end
    end
    
    % if the normalized steps are too accurate resulting in the same time
    % step for the current and next step then just continue with the next
    % data set and do not perform any computations for this steps and data
    if curTC(dataIndex) == curTN(dataIndex)
      continue;
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
        if cellDatas{ dataIndex }{ j, 5 } == curTN(dataIndex)
          pos = [ cellDatas{ dataIndex }{ j, 2 } cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
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
      cellIds1 = zeros( numCellsC(dataIndex), 3 );
      cellIds2 = zeros( numCellsN(dataIndex), 3 );
      cellPrecursors2 = cell( numCellsN(dataIndex),1 );
      
      nc1 = 1;
      nc2 = 1;
      for j=1:dimData( dataIndex )
        if cellDatas{ dataIndex }{ j, 5 } == curTC(dataIndex)
          pos = [ cellDatas{ dataIndex }{ j, 2 } cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
          matPos1(nc1, :) = pos;
          cellIds1(nc1, :) = cellDatas{ dataIndex }{ j, 1 };
          nc1 = nc1 + 1;
        end
        if cellDatas{ dataIndex }{ j, 5 } == curTN(dataIndex)
          pos = [ cellDatas{ dataIndex }{ j, 2 } cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
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
    
    % determine all cell divisions in between curTC and curTN
    divs = getDivisionsInTimeStepRange( curTC(dataIndex), curTN(dataIndex),...
      renderMasterFile, cellFileMap{dataIndex}, divisionProperties{dataIndex} );
    
    %cellDivisions{dataIndex} = divs;
    
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
      
      % draw principal components
      if renderPrincipalComponents == 1
        start = mean( cellFileMat );
        %coeff = pca( cellFileMat )
        if strcmp( visualizationType( 1, visType ), 'Ellipses' )
          start = applyTransformations( start, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
        end
        arrowLength = 150;
        % set colormap
        cm = hsv( 3 );
        %colormap( cm );
        for a=1:3
          if strcmp( visualizationType( 1, visType ), 'Ellipses' )
            endP = applyTransformations( coeff(:, a), planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
          elseif strcmp( visualizationType( 1, visType ), 'Ellipsoids' )
            endP = coeff(:, a);
          end
          hold on;
          P(a) = quiver3( start(1), start(2), start(3),...
            endP(1), endP(2), endP(3),...
            arrowLength, 'LineWidth', 3,...
            'Color', cm( a, : ) );
        end
      end
      
      if overlapping == 1
        zOffset = 0.5;
      else
        zOffset = 0.;
      end
      
      dimL = size( centerEllipse, 1 );
      if renderAverage == 1
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
          if averageOverData == 0
            tileGrid{ dataIndex, tileIndex } = [ tileGrid{ dataIndex, tileIndex }; lineDirection ];
          else
            if strcmp( termTypeStr( 1, renderTermType ), 'All' ) &&...
                renderContributions == 1
              tileGridC{ tileIndex } = [ tileGridC{ tileIndex }; [ contributions( l, 1 ) contributions( l, 2 ) ] ];
              tileGridM{ tileIndex } = [ tileGridM{ tileIndex }; [ magnitudes( l, 1 ) magnitudes( l, 2 ) ] ];
            end
            
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
        end
        
        % determine tile positions of division centers and assign division
        % orientation to the corresponding tile
        if renderDivisions == 1
          numDivs = size( divs, 1 );
          for d=1:numDivs
            time = divs( d, 2 ) + 1;
            divCenter = divs( d, 3:5 ) - centerPosPerTimeStep{dataIndex}(time, :);
            divCenter = applyTransformations( divCenter, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
            div1 = divs( d, 6:8 ) - centerPosPerTimeStep{dataIndex}(time, :);
            div2 = divs( d, 9:11 ) - centerPosPerTimeStep{dataIndex}(time, :);
            div1 = applyTransformations( div1, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
            div2 = applyTransformations( div2, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile );
            tileIndex = getTileIndex( divCenter, [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns );
            if div1(1) <= div2(1)
              divDir = [ div1(1) div1(2) div2(1) div2(2) ];
            else
              divDir = [ div2(1) div2(2) div1(1) div1(2) ];
            end
            tileGridD{ tileIndex } = [ tileGridD{ tileIndex }; divDir ];
          end
        end
        % render contour
%         if dimL > 2
%           if triangulationType == 1
%             K = convhull( timePositions( :, 4 ), timePositions( :, 5 ) );
%             CONTOUR(dataIndex) = line( timePositions( K, 4 ), timePositions( K, 5 ), timePositions( K, 6 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
%           else
%             [VolC,ShapeC] = alphavol( [ timePositions( :, 4 ), timePositions( :, 5 ) ], sqrt( alphaRadiiVector( curTC(dataIndex), 1 )) );
%             K = ShapeC.bnd(:,1);
%             dimK = size( K, 1 );
%             if dimK > 1
%               K(dimK+1,:) = K(1,:);
%               CONTOUR(dataIndex) = line( timePositions( K, 4 ), timePositions( K, 5 ), timePositions( K, 6 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
%             end
%           end
%         end
      else
        for l=1:dimL
          minSemiPoint = [minMaxSemiAxisVector( l, 1 )...
            minMaxSemiAxisVector( l, 2 )...
            minMaxSemiAxisVector( l, 3 )];
          maxSemiPoint = [minMaxSemiAxisVector( l, 4 )...
            minMaxSemiAxisVector( l, 5 )...
            minMaxSemiAxisVector( l, 6 )];
          c = centerEllipse( l, : );
          minLength = norm( minSemiPoint - c );
          maxLength = norm( maxSemiPoint - c );
          
          % line of major semi axis
          lineMaxX = [ maxSemiPoint(1), c(1) + c(1)-maxSemiPoint(1) ];
          lineMaxY = [ maxSemiPoint(2), c(2) + c(2)-maxSemiPoint(2) ];
          lineMaxZ = [ maxSemiPoint(3), c(3) + c(3)-maxSemiPoint(3) ];
          
          % line of minor semi axis
          % direction have to be here the normal of the x-y plane
          minorSemiAxes = cross( [ 0 0 1 ], maxSemiPoint-c );
          minorSemiAxes = normalizeVector3d( minorSemiAxes );
          lineMinX = [ c(1) - minLength*minorSemiAxes(1), c(1) + minLength*minorSemiAxes(1) ];
          lineMinY = [ c(2) - minLength*minorSemiAxes(2), c(2) + minLength*minorSemiAxes(2) ];
          %lineMinZ = [ c(3) - minLength*minorSemiAxes(3), c(3) + minLength*minorSemiAxes(3) ];
          lineMinZ = [ 0., 0. ];
          
          theta = acos( dot(maxSemiPoint - c, [ 1 0 0 ])/norm(maxSemiPoint - c) );
          theta = theta*180/pi;
          
          % check y values because this is required for
          % a correct rotation of the ellipses around the z axis
          if lineMaxY(1) <= lineMaxY(2)
            theta = -theta;
          end
          
          % draw ellipse
          hold on;
          ind = l+(dataIndex-1)*dimL;
          [ ELLIP(ind) ELLIPPATCH(ind) ] = drawEllipse3d( c(1), c(2), c(3)+ind*zOffset+0.2, maxLength, minLength, 0, theta );
          set( ELLIP(ind), 'color', colors( dataIndex, : ), 'LineWidth', lineWidth );
          set( ELLIPPATCH(ind), 'FaceColor', [ 1 1 1 ], 'FaceLighting', 'none' );
          
          if strcmp( lineStr( 1, lineRenderType ), 'renderLargest2DElongation' )
            colorIndex = indexColorSet( l, 2 );
            if colorIndex == 0
              color = [ 0 0 1 ];
            else
              color = [ 1 0 0 ];
            end
            
            if strcmp( termTypeStr( 1, renderTermType ), 'T' ) ||...
                strcmp( termTypeStr( 1, renderTermType ), 'All' )
              color = colors( dataIndex, : );
            end
            MAX(ind) = line( lineMaxX, lineMaxY, lineMaxZ+ind*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
            hold on;
          elseif strcmp( lineStr( 1, lineRenderType ), 'renderAll2DElongation' )
            colorIndex = indexColorSet( l, 1 );
            draw = 1;
            % negative eigenvalue
            if colorIndex == 0
              if strcmp( termTypeStr( 1, renderTermType ), 'T' )
                % do not draw a line for negative eigenvalue
                draw = 0;
              else
                color = [ 0 0 1 ];
              end
              % positive eigenvalue
            else
              if strcmp( termTypeStr( 1, renderTermType ), 'T' )
                color = [ 0 0 0 ];
              else
                color = [ 1 0 0 ];
              end
            end
            if draw == 1
              MAX(ind) = line( lineMaxX, lineMaxY, lineMaxZ+ind*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
            end
            hold on;
            colorIndex = indexColorSet( l, 2 );
            draw = 1;
            % negative eigenvalue
            if colorIndex == 0
              if strcmp( termTypeStr( 1, renderTermType ), 'T' )
                % do not draw a line for negative eigenvalue
                draw = 0;
              else
                color = [ 0 0 1 ];
              end
              % positive eigenvalue
            else
              if strcmp( termTypeStr( 1, renderTermType ), 'T' )
                color = [ 0 0 0 ];
              else
                color = [ 1 0 0 ];
              end
            end
            if draw == 1
              MIN(ind) = line( lineMinX, lineMinY, lineMinZ+ind*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
            end
            hold on;
          elseif strcmp( lineStr( 1, lineRenderType ), 'renderLargest3DElongation' )
            % get index of longest elongation
            index = minMaxEigenValueIndex( l, 3 );
            lineMaxX = [ linePos( l, (index-1)*6 +1 ), linePos( l, (index-1)*6 +4 ) ];
            lineMaxY = [ linePos( l, (index-1)*6 +2 ), linePos( l, (index-1)*6 +5 ) ];
            lineMaxZ = [ linePos( l, (index-1)*6 +3 ), linePos( l, (index-1)*6 +6 ) ];
            color = [ lineColorIndex( l, (index-1)*3 +1 ) lineColorIndex( l, (index-1)*3 +2 ) lineColorIndex( l, (index-1)*3 +3 ) ];
            if strcmp( termTypeStr( 1, renderTermType ), 'T' ) ||...
                strcmp( termTypeStr( 1, renderTermType ), 'All' )
              color = [ 0 0 0 ];
            end
            MAX(ind) = line( lineMaxX, lineMaxY, lineMaxZ+ind*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
            hold on;
          elseif strcmp( lineStr( 1, lineRenderType ), 'renderAll3DElongation' )
            % get index of all elongation types
            for el=1:3
              index = minMaxEigenValueIndex( l, el );
              lineX = [ linePos( l, (index-1)*6 +1 ), linePos( l, (index-1)*6 +4 ) ];
              lineY = [ linePos( l, (index-1)*6 +2 ), linePos( l, (index-1)*6 +5 ) ];
              lineZ = [ linePos( l, (index-1)*6 +3 ), linePos( l, (index-1)*6 +6 ) ];
              color = [ lineColorIndex( l, (index-1)*3 +1 ) lineColorIndex( l, (index-1)*3 +2 ) lineColorIndex( l, (index-1)*3 +3 ) ];
              if el == 1
                MIN(ind) = line( lineX, lineY, lineZ+ind*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
              elseif el == 2
                MID(ind) = line( lineX, lineY, lineZ+ind*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
              else
                MAX(ind) = line( lineX, lineY, lineZ+ind*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
              end
              hold on;
            end
          end
        end
      end
    else
      % if too few cells are available for PCA then just continue
      curI = curI + deltaI;
      continue;
    end
  end
  
  % determining the min and max values of contributions and magnitudes
  contributionsVector = zeros( rows*columns, 2 );
  magnitudesVector = zeros( rows*columns, 2 );
  for gt=1:rows*columns
    if size( tileGrid{ gt }, 1 ) ~= 0
      [ averageBTerm, averageTTerm ] = determineAverageTerm( tileGridC{ gt } );
      contributionsVector( gt, 1 ) = averageBTerm;
      contributionsVector( gt, 2 ) = averageTTerm;
      [ averageBMagnitude, averageTMagnitude ] = determineAverageTerm( tileGridM{ gt } );
      magnitudesVector( gt, 1 ) = averageBMagnitude;
      magnitudesVector( gt, 2 ) = averageTMagnitude;
      
      numLines = size( tileGridC{ gt }, 1 );
      for i=1:numLines
        globalBTerm( curI, 1 ) = globalBTerm( curI, 1 ) + tileGridC{ gt }( i, 1 );
        globalTTerm( curI, 1 ) = globalTTerm( curI, 1 ) + tileGridC{ gt }( i, 2 );
      end
    end
  end
  
  [ minContr, maxContr ] = determineMinMax( contributionsVector );
  [ minMag, maxMag ] = determineMinMax( magnitudesVector );

  % setting for colorbar
  %COLORBAR1 = createColorbar( 0 );
  %COLORBAR2 = createColorbar( 1 );
  %COLORBAR3 = createColorbar( 2 );
  
  % after all data are processed determine the average visualization
  if renderAverage == 1 && renderAveragePerTimeStep == 1 && averageOverData == 0
    for dataIndex=startData:endData
      for gt=1:rows*columns
        % ignore empty tiles
        if size( tileGrid{ gt }, 1 ) == 0
          continue;
        end
        averageSlope = determineAverageSlope( tileGrid{ gt } );
        L(gt) = drawAverageLines( averageSlope, gt, [totalMinAxes(1) totalMinAxes(2)],...
          [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, colors( dataIndex, : ), 0 );
      end
    end
  elseif renderAverage == 1 && renderAveragePerTimeStep == 1 && averageOverData == 1
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
            % compute the average of all contributions assigned to a
            % specific tile
            if strcmp( termTypeStr( 1, renderTermType ), 'All' ) &&...
                renderContributions == 1
              [ averageBTerm, averageTTerm ] = determineAverageTerm( tileGridC{ gt } );
              [ averageBMagnitude, averageTMagnitude ] = determineAverageTerm( tileGridM{ gt } );
              colorB = determineInterpolatedColor( averageBTerm, minContr, maxContr, 0 );
              colorT = determineInterpolatedColor( averageTTerm, minContr, maxContr, 0 );
              if averageBMagnitude > 0
                colorBM = determineInterpolatedColor( averageBMagnitude, 0, maxMag, 1 );
              elseif averageBMagnitude == 0.
                colorBM = [ 1 1 1 ];
              else
                colorBM = determineInterpolatedColor( averageBMagnitude, minMag, 0, 2 );
              end
              if averageTMagnitude > 0
                colorTM = determineInterpolatedColor( averageTMagnitude, 0, maxMag, 1 );
              elseif averageTMagnitude == 0.
                colorTM = [ 1 1 1 ];
              else
                colorTM = determineInterpolatedColor( averageTMagnitude, minMag, 0, 2 );
              end
              
              BC(lineRenderIndex) = drawContributionRectangle( gt, [totalMinAxes(1) totalMinAxes(2)],...
                [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, colorB, 0 );
              TC(lineRenderIndex) = drawContributionRectangle( gt, [totalMinAxes(1) totalMinAxes(2)],...
                [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, colorT, 1 );
              BM(lineRenderIndex) = drawContributionRectangle( gt, [totalMinAxes(1) totalMinAxes(2)],...
                [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, colorBM, 2 );
              TM(lineRenderIndex) = drawContributionRectangle( gt, [totalMinAxes(1) totalMinAxes(2)],...
                [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, colorTM, 3 );
            end
            L(lineRenderIndex) = drawAverageLines( averageSlope, gt, [totalMinAxes(1) totalMinAxes(2)],...
              [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 0 0 ], 0 );
            sd = determineSDDirections( tileGrid{ gt } );
            if sd ~= 0
              SD(lineRenderIndex) = drawStandardDeviationArea( sd, averageSlope, gt, [totalMinAxes(1) totalMinAxes(2)],...
                [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 0 0 ] );
            end
          end
        end
        
        % rendering of division orientations
        if renderDivisions == 1
          if size( tileGridD{ gt }, 1 ) ~= 0
            averageSlopeD = determineAverageSlope( tileGridD{ gt } );
            DIV(lineRenderIndex) = drawAverageLines( averageSlopeD, gt, [totalMinAxes(1) totalMinAxes(2)],...
              [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 1 0 ], 0 );
            sd = determineSDDirections( tileGridD{ gt } );
            if sd ~= 0
              SDD(lineRenderIndex) = drawStandardDeviationArea( sd, averageSlopeD, gt, [totalMinAxes(1) totalMinAxes(2)],...
                [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 1 0 ] );
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
        
        % rendering of division orientations
        if renderDivisions == 1
          numLines = size( tileGridD{ gt }, 1 );
          for l=1:numLines
            % first compute the slope
            startPos = tileGridD{ gt }(l, 1:2);
            endPos = tileGridD{ gt }(l, 3:4);
            slope = (endPos(2)-startPos(2))/(endPos(1)-startPos(1));
            DIV(lineRenderIndex) = drawAverageLines( slope, gt, [totalMinAxes(1) totalMinAxes(2)],...
              [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, [ 0 1 0 ], 0 );
            lineRenderIndex = lineRenderIndex + 1;
          end
        end
        
      end
    end
  end
  
  hold off;
  set( f,'nextplot','replacechildren' );
  viewOffset = 100;
  axis( [ totalMinAxes(1) totalMaxAxes(1)...
    totalMinAxes(2) totalMaxAxes(2)...
    totalMinAxes(3)-viewOffset totalMaxAxes(3)+viewOffset ...
    totalMinAxes(3) totalMaxAxes(3) ] );
  axis on
  daspect( [ 1 1 1 ] );
  
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
  
  if strcmp( visualizationType( 1, visType ), 'Ellipsoids' )
    % first delete last lightsource
    delete(findall(gcf, 'Type', 'light'))
    camlight headlight;
  end
  
  % if images instead of a video should be exported
  if strcmp( exportTypeStr( 1, exportType ), 'AsImages' ) &&...
      renderAveragePerTimeStep == 1
    if curI < 10
      digit = strcat( viewStr( 1, cView ), '_00' );
    elseif curI < 100
      digit = strcat( viewStr( 1, cView ), '_0' );
    else
      digit = strcat( viewStr( 1, cView ), '_' );
    end
    
    % output with number of cells
    filePath = strcat( imageDir, digit, num2str(curI), '-Cells', num2str(numNormCellsN), '.png' );
    
    saveas( gcf, char(filePath) );
  end
  % at last increment the current index loop parameter
  curI = curI + deltaI;
end

% after all data are processed determine the average visualization
if renderAverage == 1 && renderAveragePerTimeStep == 0 && averageOverData == 0
  hold off;
  set( f,'nextplot','replacechildren' );
  viewOffset = 100;
  axis( [ totalMinAxes(1) totalMaxAxes(1)...
    totalMinAxes(2) totalMaxAxes(2)...
    totalMinAxes(3)-viewOffset totalMaxAxes(3)+viewOffset ...
    totalMinAxes(3) totalMaxAxes(3) ] );
  axis on
  daspect( [ 1 1 1 ] );
  
  % legend
  hold on;
  %leg = legend( '120830', '121204', '121211', '130508', '130607', '131203' );
  leg = legend( '120830', '121204', '121211', '130508', '130607' );
  set(leg, 'Location', 'NorthWestOutside');
  linecolors = { [ 1 0 1 ], [ 0 0 0 ], [ 1 0 0 ], [ 0 1 0 ], [ 0 0 1 ] };
  legendlinestyles( leg, {}, {}, linecolors );
  
  for dataIndex=startData:endData
    for gt=1:rows*columns
      % ignore empty tiles
      if size( tileGrid{ dataIndex, gt }, 1 ) == 0
        continue;
      end
      averageSlope = determineAverageSlope( tileGrid{ dataIndex, gt } );
        L(gt) = drawAverageLines( averageSlope, gt, [totalMinAxes(1) totalMinAxes(2)],...
          [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, colors( dataIndex, : ), 0 );
    end
  end
  
  if strcmp( exportTypeStr( 1, exportType ), 'AsImages' )
    filePath = strcat( imageDir, pureDataStr{dataIndex}, '-AverageLines.png' );
    saveas( gcf, char(filePath) );
  end
end

ElapsedTimeIndexLoop = cputime - cpuT
