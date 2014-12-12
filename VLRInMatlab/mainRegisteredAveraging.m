geomPath = strcat( pwd, '/geom3d' );
addpath( geomPath );

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
startI = 149;
% endIndex
endI = 150;
% min and max index
minI = 1;
maxI = 150;
% draw delaunay tri?
drawDelaunay = 0;
% draw average delaunay tri?
drawAverageDelaunay = 0;
% draw point set as ellipses
drawNuclei = 0;
% draw contour of cell nuclei
drawContour = 0;
% draw average contour of cell nuclei
drawAverageContour = 1;
% include contour points
includeContourPoints = 1;
% index for color of average stuff
averageColorIndex = 6;
% render average nuclei positions
drawAverageNuclei = 1;
% render average nuclei positions at nextT
drawNextAverageNuclei = 0;
% draw linking between cells at curT and nextT
drawLinking = 0;
% if set to one then only a single time step
% is rendered given by startT
exportType = 2;
% vector of data strings
exportTypeStr = { 'SingleFigure' 'AsImages' };
% exclude nuclei outliers
excludeOutliers = 1;
% render only master file?
renderMasterFile = 1;
% render principal components
renderPrincipalComponents = 0;
% line width of ellipses and semi axes
lineWidth = 1.2;
% enable z overlapping
overlapping = 1;
renderCellRanges = 0;
% only generate output images if the number of cells are within a desired
% range
cellRange = [ 20 40 60 80 100 120 140 ];
% use triangulation based on delaunay or alpha shape
% 1 -> delaunay
% 2 -> alpha shape
triangulationType = 1;
% vector of data strings
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' '131203_raw' };
pureDataStr = { '120830' '121204' '121211' '130508' '130607' 'Average' };
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
% number of cell files
numCellFiles = 6;%double( maxCF-minCF+1 );
% resolution of grid
resGrid = 50;

% define offset for boundary points which is the distance between the
% acutal location of the contour points and their initial position within
% the model
conOffset = 20;
% distance from contour point to nearest nuclei position
cellDist = 25;
% generate contour points automatically or use the points set manually
autoContour = 0;

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
  prepareData( dataStr, startData, endData, numData, 'Ellipses', renderMasterFile, cView, 1 );

if drawContour == 1 || drawAverageContour == 1
  CONTOUR = [];
end
if drawNuclei == 1
  % ellipse instance
  ELLIP = [];
  ELLIPPATCH = [];
end

if drawLinking == 1
  LINKING = [];
end

if drawAverageNuclei == 1
  AV = [];
  AVPATCH = [];
end

if drawNextAverageNuclei == 1
  AVN = [];
  AVNPATCH = [];
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

[ rows, columns ] = generate2DGrid( [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid );

if drawAverageNuclei == 1 || drawAverageContour == 1 || drawAverageDelaunay == 1
  if autoContour == 1
    fileName = strcat( '/tmp/triangulation-', 'Average', '_auto.txt' );
  else
    fileName = strcat( '/tmp/triangulation-', 'Average', '.txt' );
    % contour points for the first and last time step
    cPointsFirst = generateContourPoints( 'Average', true, conOffset );
    cPointsLast = generateContourPoints( 'Average', false, conOffset );
  end
  fileId = fopen( char(fileName), 'w' );
  % first write the maximum number of time steps
  fprintf( fileId, '%1d\n', startI );
  fprintf( fileId, '%1d\n', maxI );
end

% number of total cells per step
numTotalCells = 1;
links = 1;

cpuT = cputime;
% loop over all normalized steps
for curI=startI:endI
  
  % initialize tiles in the grid
  tileGrid = cell( rows*columns, 1 );
  nextTileGrid = cell( rows*columns, 1 );
  
  if drawContour == 1 || drawAverageContour == 1
    if ishandle( CONTOUR )
      set( CONTOUR, 'Visible', 'off' );
    end
  end
  
  if drawNuclei == 1
    if ishandle( ELLIP )
      set( ELLIP, 'Visible', 'off' );
    end
    if ishandle( ELLIPPATCH )
      set( ELLIPPATCH, 'Visible', 'off' );
    end
  end
  
  if drawLinking == 1
    if ishandle( LINKING )
      set( LINKING, 'Visible', 'off' );
    end
  end
  
  if drawAverageNuclei == 1
    if ishandle( AV )
      set( AV, 'Visible', 'off' );
    end
    if ishandle( AVPATCH )
      set( AVPATCH, 'Visible', 'off' );
    end
  end
  
  if drawNextAverageNuclei == 1
    if ishandle( AVN )
      set( AVN, 'Visible', 'off' );
    end
    if ishandle( AVNPATCH )
      set( AVNPATCH, 'Visible', 'off' );
    end
  end
  
  if renderPrincipalComponents == 1
    if ishandle( P )
      set( P, 'Visible', 'off' );
    end
  end
  
  % loop over all data sets
  allCurT = zeros( numData, 1 );
  allCells = zeros( numData, 1 );
  allCurCells = zeros( numData, 1 );
  %averagePos = [];
  %nextAveragePos = [];
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
    
    if strcmp( dataStr( 1, dataIndex ), '130607_raw' )
      TF = createRotationOz( -pi/36. ) * TF;
    end
    
    % get the alpha shape radii for all time steps
    if triangulationType == 2
      alphaRadiiVector = getAlphaRadius( dataStr( 1, dataIndex ) );
    end
    
    % get the corresponding time steps for the reg. step
    [ curT, numNormCells ] = getCorrespondingTimeStep( curI, minI, maxI,...
      maxT(dataIndex), numCellsPerTimeStep{dataIndex} );
    
    % store the chosen time step for each data set
    allCurT( dataIndex, 1 ) = curT;
    allCells( dataIndex, 1 ) = numCellsPerTimeStep{dataIndex}(curT,1);
    
    [ nextT, numNextNormCells ]  = getCorrespondingTimeStep( curI+1, minI, maxI,...
      maxT(dataIndex), numCellsPerTimeStep{dataIndex} );
    
    % matrix of positions for current and next time step
    centerEllipse = [];
    averagePos = [];
    nextAveragePos = [];
    
    % determine cell position and cell file information
    for j=1:dimData( dataIndex )
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
      
      if cellDatas{ dataIndex }{ j, 5 } == curT
        % get position of current cell
        p = [ cellDatas{ dataIndex }{ j, 2 }...
          cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
        p = p - centerPosPerTimeStep{dataIndex}(curT,:);
        p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile, curI );
        centerEllipse = [ centerEllipse ; p ];
        % determine tileIndex of current nuclei position and add it to the
        % tile grid in order to average all positions later
        tileIndex = getTileIndex( p, [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns );
        tileGrid{ tileIndex } = [ tileGrid{ tileIndex }; [ p(1) p(2) p(3) dataIndex ] ];
        averagePos = [ averagePos ; p(1), p(2) 0 ];
        
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
      end
      if cellDatas{ dataIndex }{ j, 5 } == nextT
        % get position of next cell
        p = [ cellDatas{ dataIndex }{ j, 2 }...
          cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
        p = p - centerPosPerTimeStep{dataIndex}(nextT,:);
        p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile, curI+1 );
        % determine tileIndex of next nuclei position and add it to the
        % tile grid in order to average all positions later
        nextTileIndex = getTileIndex( p, [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns );
        nextTileGrid{ nextTileIndex } = [ nextTileGrid{ nextTileIndex }; [ p(1) p(2) p(3) dataIndex ] ];
        nextAveragePos = [ nextAveragePos ; p(1), p(2) 0 ];
      end
    end
    
    numCells = size( centerEllipse, 1 );
    % store the total number of cells
    allCurCells( dataIndex, 1 ) = numCells;
    
    % if at least three cells exists
    if numCells > 3 && drawDelaunay == 1
      if triangulationType == 1
        % delaunay triangulation
        tri = delaunayTriangulation( centerEllipse(:,1), centerEllipse(:,2), centerEllipse(:,3) );
        tetramesh( tri, 'FaceColor', colors( dataIndex, : ), 'FaceAlpha', 0.9 );
      else
        % alpha shape triangulation
        [Vol,Shape] = alphavol( [ centerEllipse(:,1) centerEllipse(:,2) centerEllipse(:,3) ], sqrt( alphaRadiiVector( curT, 1 )) );
        tri = Shape.tri;
        trisurf( tri, centerEllipse(:,1), centerEllipse(:,2), centerEllipse(:,3),...
          'FaceColor', colors( dataIndex, : ), 'FaceAlpha', 0. );
      end
    end
    
    % draw principal components
    if renderPrincipalComponents == 1
      start = mean( centerEllipse );
      %coeffMat = pca( centerEllipse )
      start = applyTransformations( start, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile, curI );
      arrowLength = 150;
      for a=1:3
        %endP = applyTransformations( coeff(:, a), planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile, curI );
        endP = coeff( :, a );
        P(a) = quiver3( start(1), start(2), start(3),...
          endP(1), endP(2), endP(3),...
          arrowLength, 'LineWidth', 3,...
          'Color', cm( a, : ) );
      end
    end
    
    if drawContour == 1
      dimL = size( centerEllipse, 1 );
      if dimL > 2
        if triangulationType == 1
          K = convhull( centerEllipse( :, 1 ), centerEllipse( :, 2 ) );
          CONTOUR(curI) = line( centerEllipse( K, 1 ), centerEllipse( K, 2 ), centerEllipse( K, 3 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
        else
          [VolC,ShapeC] = alphavol( [ centerEllipse( :, 1 ), centerEllipse( :, 2 ) ], sqrt( alphaRadiiVector( curT, 1 )) );
          K = ShapeC.bnd(:,1);
          dimK = size( K, 1 );
          if dimK > 1
            K(dimK+1,:) = K(1,:);
            CONTOUR(curI) = line( centerEllipse( K, 1 ), centerEllipse( K, 2 ), centerEllipse( K, 3 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
          end
        end
      end
    end
  end
  
  % compute the average of nuclei positions in each tile
  if overlapping == 1
    zOffset = 0.5;
  else
    zOffset = 0.;
  end
  averagePos = [];
  nextAveragePos = [];
  for t=1:rows*columns
    hold on
    % next step
    if curI ~= maxI
      if size( nextTileGrid{ t }, 1 ) ~= 0
        % compute the average of all positions in the tile
        sum = zeros( numData, 3 );
        sizes = zeros( numData, 1 );
        for i=1:size( nextTileGrid{ t }, 1 )
          dIndex = nextTileGrid{ t }(i, 4);
          sum( dIndex, : ) = sum(dIndex, : ) + [ nextTileGrid{ t }(i, 1) nextTileGrid{ t }(i, 2) nextTileGrid{ t }(i, 3)];
          sizes( dIndex, 1 ) = sizes( dIndex, 1 ) + 1;
        end
        for d=1:numData
          if sizes( d, 1 ) ~= 0
            pos = sum( d, : )/sizes( d, 1 );
            nextAveragePos = [ nextAveragePos ; pos(1), pos(2) 0 ];
            if drawNextAverageNuclei == 1
              [ AVN(numTotalCells), AVNPATCH(numTotalCells) ] = drawEllipse3d( pos(1), pos(2), pos(3)+zOffset+0.2, radEllip, radEllip, 0, 0 );
              set( AVN(numTotalCells), 'color', colors( averageColorIndex+1, : ), 'LineWidth', lineWidth );
              set( AVNPATCH(numTotalCells), 'FaceColor', colors( averageColorIndex+1, : ), 'FaceLighting', 'none' );
            end
          end
        end        
      end
    end
  end
  
  numNextAverageCells = size( nextAveragePos, 1 );
  for t=1:rows*columns
    hold on
    if size( tileGrid{ t }, 1 ) ~= 0
      % compute the average of all positions in the tile
      sum = zeros( numData, 3 );
      sizes = zeros( numData, 1 );
      for i=1:size( tileGrid{ t }, 1 )
        dIndex = tileGrid{ t }(i, 4);
        sum( dIndex, : ) = sum(dIndex, : ) + [ tileGrid{ t }(i, 1) tileGrid{ t }(i, 2) tileGrid{ t }(i, 3)];
        sizes( dIndex, 1 ) = sizes( dIndex, 1 ) + 1;
      end
      for d=1:numData
        if sizes( d, 1 ) ~= 0
          pos = sum( d, : )/sizes( d, 1 );
          averagePos = [ averagePos ; pos(1), pos(2) 0 ];
          if size( averagePos, 1 ) == numNextAverageCells
            break;
          end
          if drawAverageNuclei == 1
            [ AV(numTotalCells), AVPATCH(numTotalCells) ] = drawEllipse3d( pos(1), pos(2), pos(3)+zOffset+0.2, radEllip, radEllip, 0, 0 );
            set( AV(numTotalCells), 'color', colors( averageColorIndex, : ), 'LineWidth', lineWidth );
            set( AVPATCH(numTotalCells), 'FaceColor', colors( averageColorIndex, : ), 'FaceLighting', 'none' );
            numTotalCells = numTotalCells+1;
          end
        end
      end
      if size( averagePos, 1 ) == numNextAverageCells
        break;
      end
    end
  end
  
  numAverageCells = size( averagePos, 1 );
  
  % else apply a new tracking approach to map cells from time step curT
  % to cells at time step nextT such that the number of cells at both
  % time steps are identical; in fact there are three cases to handle:
  
  % first case: the number of cells are identical -> we apply a mapping
  % from one cell at curT to its nearest neighbor in the time step nextT
  if numAverageCells == numNextAverageCells
    nextLinkedPos = zeros( numNextAverageCells, 3 );
    linkedPos = averagePos;
    IDX = knnsearch( nextAveragePos, averagePos );
    for l=1:size( IDX, 1 )
      nextLinkedPos( l, : ) = nextAveragePos( IDX( l, 1 ), : );
      if drawLinking == 1
        linePos = [ averagePos( l, 1 ) averagePos( l, 2 ) ; nextAveragePos( IDX( l, 1 ), 1 ) nextAveragePos( IDX( l, 1 ), 2 ) ];
        LINKING(links) = line( linePos( :, 1 ), linePos( :, 2 ), 'Color', colors( averageColorIndex+1, : ), 'LineWidth', lineWidth );
        links = links+1;
      end
    end
  % second case: we merge and average nearest neighbor at time step nextT until
  % the number of cells are identical
  elseif numAverageCells < numNextAverageCells
    nextLinkedPos = zeros( numAverageCells, 3 );
    linkedPos = averagePos;
    newNumNext = numNextAverageCells;
    while numAverageCells ~= newNumNext
      IDX = knnsearch( nextAveragePos, averagePos, 'K', 2 );
      allDist = zeros( size( IDX, 1 ), 1 );
      for l=1:size( IDX, 1 )
        allDist( l, : ) = distancePoints3d( nextAveragePos( IDX( l, 1 ), : ), nextAveragePos( IDX( l, 2 ), : ) );
      end
      [ dist, I ] = sort( allDist );
      indexOfSmallestEntry = I(1,1);
      k1 = IDX( indexOfSmallestEntry, 1 );
      k2 = IDX( indexOfSmallestEntry, 2 );
      % compute the average of the points with the smallest distance
      newMeanPos = (nextAveragePos( k1, : ) + nextAveragePos( k2, : ))/2.;
      % set the new averaged position and delete the second entry
      nextAveragePos( k1, : ) = newMeanPos;
      nextAveragePos( k2, : ) = [];
      newNumNext = size( nextAveragePos, 1 );
    end
    % finally find the nearest neighbor to apply the linking
    IDX = knnsearch( nextAveragePos, averagePos );
    for l=1:size( IDX, 1 )
      nextLinkedPos( l, : ) = nextAveragePos( IDX( l, 1 ), : );
      if drawLinking == 1
        linePos = [ averagePos( l, 1 ) averagePos( l, 2 ) ; nextAveragePos( IDX( l, 1 ), 1 ) nextAveragePos( IDX( l, 1 ), 2 ) ];
        LINKING(links) = line( linePos( :, 1 ), linePos( :, 2 ), 'Color', colors( averageColorIndex+1, : ), 'LineWidth', lineWidth );
        links = links+1;
      end
    end
  % third case: (numAverageCells > numNextAverageCells) we merge and
  % average nearest neighbor at time step curT until the number of cells 
  % are identical
  else
    if numNextAverageCells ~= 0
      disp( strcat( 'There are fewer cells in the next time step from ', num2str(numAverageCells), ' > ', num2str(numNextAverageCells) ) );
      %linkedPos = averagePos;
      %nextLinkedPos = nextAveragePos;
      linkedPos = zeros( numNextAverageCells, 3 );
      nextLinkedPos = nextAveragePos;
      newNum = numAverageCells;
      while numNextAverageCells ~= newNum
        IDX = knnsearch( averagePos, nextAveragePos, 'K', 2 );
        allDist = zeros( size( IDX, 1 ), 1 );
        for l=1:size( IDX, 1 )
          allDist( l, : ) = distancePoints3d( averagePos( IDX( l, 1 ), : ), averagePos( IDX( l, 2 ), : ) );
        end
        [ dist, I ] = sort( allDist );
        indexOfSmallestEntry = I(1,1);
        k1 = IDX( indexOfSmallestEntry, 1 );
        k2 = IDX( indexOfSmallestEntry, 2 );
        % compute the average of the points with the smallest distance
        newMeanPos = (averagePos( k1, : ) + averagePos( k2, : ))/2.;
        % set the new averaged position and delete the second entry
        averagePos( k1, : ) = newMeanPos;
        averagePos( k2, : ) = [];
        newNum = size( averagePos, 1 );
      end
      % finally find the nearest neighbor to apply the linking
      IDX = knnsearch( averagePos, nextAveragePos );
      for l=1:size( IDX, 1 )
        linkedPos( l, : ) = averagePos( IDX( l, 1 ), : );
        if drawLinking == 1
          linePos = [ nextAveragePos( l, 1 ) nextAveragePos( l, 2 ) ; averagePos( IDX( l, 1 ), 1 ) averagePos( IDX( l, 1 ), 2 ) ];
          LINKING(links) = line( linePos( :, 1 ), linePos( :, 2 ), 'Color', colors( averageColorIndex+1, : ), 'LineWidth', lineWidth );
          links = links+1;
        end
      end
    end
  end
  
  % draw contour of data and the single marks
  if includeContourPoints == 1 &&...
      (drawAverageNuclei == 1 || drawAverageContour == 1 || drawAverageDelaunay == 1)
    if autoContour == 0
      factor = curI/maxI;
      cPoints = (1-factor) * cPointsFirst + factor * cPointsLast;
    else
      cPoints = generateAutomaticContourPoints( linkedPos, cellDist, conOffset,...
        curI == startI, 'Average' );
    end
  end
  
  % include and render contour points
  if includeContourPoints == 1
    for cc=1:size( cPoints, 1 )
      hold on
      [ AV(rows*columns+cc), AVPATCH(rows*columns+cc) ] = drawEllipse3d( cPoints(cc,1), cPoints(cc,2), cPoints(cc,3)+zOffset+0.2,...
        radEllip, radEllip, 0, 0 );
      set( AV(rows*columns+cc), 'color', 'k', 'LineWidth', lineWidth );
      set( AVPATCH(rows*columns+cc), 'FaceColor', 'k', 'FaceLighting', 'none' );
      
      linkedPos = [ linkedPos; cPoints( cc, : ) ];
    end
  end
  
  % generate triangulation of average point set
  if triangulationType == 1
    % delaunay triangulation
    averageTri = delaunayTriangulation( linkedPos(:,1), linkedPos(:,2) );
  else
    % alpha shape triangulation
    [Vol,Shape] = alphavol( [ linkedPos(:,1) linkedPos(:,2) ], 85 );
    averageTri = Shape.tri;
  end
  
  % if at least three cells exists
  if drawAverageDelaunay == 1 && size( linkedPos, 1 ) > 2
    if triangulationType == 1
      triplot( averageTri, 'b' );
      % draw triangle labels in the center of each triangle
      %TEXT = drawTriangleLabels( averageTri );
    else
      trisurf( averageTri, linkedPos(:,1), linkedPos(:,2), linkedPos(:,3),...
        'FaceColor', 'blue', 'FaceAlpha', 0. );
    end
  end
  
  if drawAverageContour == 1
    dimL = size( linkedPos, 1 );
    if dimL > 2
      if triangulationType == 1
        K = convhull( linkedPos( :, 1 ), linkedPos( :, 2 ) );
        CONTOUR(curI) = line( linkedPos( K, 1 ), linkedPos( K, 2 ), 'Color', colors( averageColorIndex, : ), 'LineWidth', lineWidth );
      else
        [VolC,ShapeC] = alphavol( [ linkedPos( :, 1 ), linkedPos( :, 2 ) ], 85 );
        K = ShapeC.bnd(:,1);
        dimK = size( K, 1 );
        if dimK > 1
          K(dimK+1,:) = K(1,:);
          CONTOUR(curI) = line( linkedPos( K, 1 ), linkedPos( K, 2 ), 'Color', colors( averageColorIndex, : ), 'LineWidth', lineWidth );
        end
      end
    end
  end
  
  if includeContourPoints == 1
    % export triangulation properties
    exportTriangulation( averageTri, linkedPos, curI, 'Average', triangulationType, autoContour );
    
    if curI ~= maxI
      if autoContour == 0
        fac = (curI+1)/maxI;
        nextInterPoints = (1-fac) * cPointsFirst + fac * cPointsLast;
        exportNewPosOfTriangulation( nextInterPoints, nextLinkedPos, 'Average', autoContour );
      else
        newCPoints = generateAutomaticContourPoints( nextLinkedPos, cellDist, conOffset,...
          false, 'Average' );
        exportNewPosOfTriangulation( newCPoints, nextLinkedPos, 'Average', autoContour );
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
  
  % legend
  stringData = [];
  for d=1:numData
    stringData = [ stringData ; strcat( pureDataStr( 1, d ),...
      ' - T', num2str(allCurT(d, 1)),...
      ' - C', num2str(allCurCells(d, 1)), '/', num2str(allCells(d, 1)) ) ];
  end
  stringData = [ stringData ; strcat( 'Average - C', num2str( numAverageCells ) ) ];
  str = cellstr( stringData );
  %leg = legend( '120830', '121204', '121211', '130508', '130607', 'Average' );
  leg = legend( str );
  set(leg, 'Location', 'NorthWestOutside');
  linecolors = { [ 1 0 1 ], [ 1 1 0 ], [ 1 0 0 ], [ 0 1 0 ], [ 0 0 1 ], [ 0 1 1 ] };
  legendlinestyles( leg, {}, {}, linecolors );
  
  grid off;
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  title( strcat( 'Normalized Step ', num2str(curI), '-Cells', num2str(numNormCells) ) );
  
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
    filePath = strcat( imageDir, digit, num2str(curI), '-Cells', num2str(numNormCells), '.png' );
    
    saveas( gcf, char(filePath) );
  end
end
ElapsedTimeIndexLoop = cputime - cpuT;
disp( strcat( 'Elapsed time: ', num2str(ElapsedTimeIndexLoop), ' sec.' ) );

