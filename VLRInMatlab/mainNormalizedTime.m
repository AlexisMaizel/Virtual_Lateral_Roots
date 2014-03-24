%%%% The mesh deformation is based on the paper of Graner et al,
%%%% Discrete rearranging disordered patterns, part I: Robust
%%%% statistical tools in two or three dimensions, Eur. Phys. J.,
%%%% pp 349 -- 369, 2008

addpath( '/home/necrolyte/Data/VLR/Virtual_Lateral_Roots/VLRInMatlab/geom3d' );

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
renderTermType = 1;
termTypeStr = { 'B' 'T' 'All' };
% startIndex
startI = 1;
% endIndex
endI = 20;
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
renderSingleCellFile = 1;
% render principal components
renderPrincipalComponents = 0;
% vector of view offsets for x/Y/Z coordinates
viewOffsets = [ 100 100 100 250 250 250 ];
% line width of ellipses and semi axes
lineWidth = 1.2;
% enable z overlapping
overlapping = 1;
% only render specific ranges of cell numbers
% render contour (=convexHull) instead of ellipses
visualizationType = { 'Ellipsoids' 'Ellipses' 'Contour' };
visType = 2;
% offset for cell ranges
epsilon = 3;
lineRenderType = 3;
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
% master cell file information taken by the picture made of Daniel located in dropbox
% and the trackGroup information of the raw data sets
%masterCellFile = [ 4 3 4 2 3 0 ];
masterCellFile = [ 4 5 4 3 3 0 ];
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
% number of cell files
numCellFiles = 6;%double( maxCF-minCF+1 );
% resolution of grid
resGrid = 30;
% render average lines or not
renderAverage = 1;
renderAveragePerTimeStep = 0;

% figure properties
f = figure( 'Name', 'Mesh Deformation', 'Position', [100 100 800 800] );
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
[ cellDatas, dimData, maxT, numCellsPerTimeStep, centerPosPerTimeStep, totalMinAxes, totalMaxAxes, cellFileMap ] =...
  prepareData( dataStr, startData, endData, numData, visualizationType( 1, visType ), renderSingleCellFile, cView );

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

% set colormap and colorbar depending on the number of cell files
cm = hsv( numCellFiles );
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

cpuT = cputime;

hold on;
[ rows, columns ] = generate2DGrid( [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid );
% initialize grid size consisting of tiles
tileGrid = cell( numData, rows*columns );

% loop over all normalized steps
curI=startI;
while curI < endI-deltaI+1
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
%     if ishandle( L )
%       set( L, 'Visible', 'off' );
%     end
    if ishandle( CONTOUR )
      set( CONTOUR, 'Visible', 'off' );
    end
  end
  
  if renderPrincipalComponents == 1
    if ishandle( P )
      set( P, 'Visible', 'off' );
    end
  end
  
  % loop over all data sets
  for dataIndex=startData:endData
    % get stored eigenvectors for the last time step to set the same
    % direction view for each time step
    coeff = getNormalizedPrincipalComponents( dataStr( 1, dataIndex ), renderSingleCellFile );
    
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
      if strcmp( dataStr( 1, dataIndex ), '120830_raw' ) &&...
          renderSingleCellFile == 0
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
    
    % render the corresponding master cell file
    singleCellFile = masterCellFile( 1, dataIndex );
    
    % get the alpha shape radii for all time steps
    if triangulationType == 2
      alphaRadiiVector = getAlphaRadius( dataStr( 1, dataIndex ) );
    end
    
    % NEW: only continue with this time step that is synchronized in
    % number of cells for each data set
    if begin ~= 1
      numNormCellsN = getNormalizedCellNumber( curI+deltaI, 18, 143, 1, 20 );
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
      numNormCellsC = getNormalizedCellNumber( curI, 18, 143, 1, 20 );
      numNormCellsN = getNormalizedCellNumber( curI+deltaI, 18, 143, 1, 20 );
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
        elseif cellDatas{ dataIndex }{ j, 5 } == curTN(dataIndex)
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
        timePositions, indexColorSet ]...
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
        renderSingleCellFile,...
        singleCellFile,...
        cellFileMap{dataIndex} );
      
%      
%       TODO in computeTimeEvolution for speedup
%       cellFileMat = zeros( numConsideredCells, 3 );
%       linePos = zeros( numConsideredCells, 18 );
%       minMaxSemiAxisVector = zeros( numConsideredCells, 6 );
%       centerEllipse = zeros( numConsideredCells, 3 );
%       
      % draw principal components
      if renderPrincipalComponents == 1
        start = mean( cellFileMat );
        %coeff = pca( cellFileMat )
        if strcmp( visualizationType( 1, visType ), 'Ellipses' ) ||...
            strcmp( visualizationType( 1, visType ), 'Contour' )
          start = applyTransformations( start, planePos, u, v, TF, dataStr( 1, dataIndex ) );
        end
        arrowLength = 150;
        for a=1:3
          if strcmp( visualizationType( 1, visType ), 'Ellipses' ) ||...
              strcmp( visualizationType( 1, visType ), 'Contour' )
            endP = applyTransformations( coeff(:, a), planePos, u, v, TF, dataStr( 1, dataIndex ) );
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
      
      if renderAverage == 1
        dimL = size( centerEllipse, 1 );
        for l=1:dimL
          maxSemiPoint = [minMaxSemiAxisVector( l, 4 )...
            minMaxSemiAxisVector( l, 5 )...
            minMaxSemiAxisVector( l, 6 )];
          c = centerEllipse( l, : );
          
          % line coords of major semi axis
          lineX = [ maxSemiPoint(1), c(1) + c(1)-maxSemiPoint(1) ];
          lineY = [ maxSemiPoint(2), c(2) + c(2)-maxSemiPoint(2) ];
          
          % check which start/end point of the line is nearest to pStart or
          % pEnd of the two cell positions that are compared between the
          % two time steps/indexes deltaT/deltaI
          dist1 = distancePoints3d( [ timePositions(l,1) timePositions(l,2) 0.0 ],...
            [ lineX(1) lineY(1) 0.0 ] );
          dist2 = distancePoints3d( [ timePositions(l,1) timePositions(l,2) 0.0 ],...
            [ lineX(2) lineY(2) 0.0 ] );
          
          if dist1 < dist2
            lineDirection = [ lineX(2)-lineX(1) lineY(2)-lineY(1) 0.0 ];
          else
            lineDirection = [ lineX(1)-lineX(2) lineY(1)-lineY(2) 0.0 ];
          end
          
          % determine tileIndex of current ellipse position and add it to the
          % tile grid in order to average all deformations occurring in each
          % tile
          tileIndex = getTileIndex( c, [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns );
          tileGrid{ dataIndex, tileIndex } = [ tileGrid{ dataIndex, tileIndex }; lineDirection ];
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
        dimL = size( centerEllipse, 1 );
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
  
  % after all data are processed determine the average visualization
  if renderAverage == 1 && renderAveragePerTimeStep == 1
    for dataIndex=startData:endData
      for gt=1:rows*columns
        averageDirection = determineAverageDirection( tileGrid{ dataIndex, gt } );
      if all(averageDirection == 0)
        continue;
      else
        L(gt) = drawAverageLines( averageDirection, gt, [totalMinAxes(1) totalMinAxes(2)],...
          [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, colors( dataIndex, : ) );
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
  
  grid off;
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  title( strcat( 'Normalized Step ', num2str(curI+deltaI), '-Cells', num2str(numNormCellsN) ) );
  
  if strcmp( visualizationType( 1, visType ), 'Ellipsoids' )
    % first delete last lightsource
    delete(findall(gcf, 'Type', 'light'))
    camlight headlight;
  end
  
  % if images instead of a video should be exported
  if strcmp( exportTypeStr( 1, exportType ), 'AsImages' )
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
if renderAverage == 1 && renderAveragePerTimeStep == 0
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
  
  for dataIndex=startData:endData
    for gt=1:rows*columns
      averageDirection = determineAverageDirection( tileGrid{ dataIndex, gt } );
      if all(averageDirection == 0)
        continue;
      else
        L(gt) = drawAverageLines( averageDirection, gt, [totalMinAxes(1) totalMinAxes(2)],...
          [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, colors( dataIndex, : ) );
      end
    end
  end
  
  if strcmp( exportTypeStr( 1, exportType ), 'AsImages' )
    filePath = strcat( imageDir, pureDataStr{dataIndex}, '-AverageLines.png' );
    saveas( gcf, char(filePath) );
  end
end

ElapsedTimeIndexLoop = cputime - cpuT
