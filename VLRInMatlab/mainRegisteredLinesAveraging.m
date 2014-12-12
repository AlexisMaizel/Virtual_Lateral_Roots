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
startI = 90;
% endIndex
endI = 90;
% min and max index
minI = 1;
maxI = 100;
% draw delaunay tri?
drawDelaunay = 0;
% draw average delaunay tri?
drawAverageDelaunay = 1;
% draw point set as ellipses
drawNuclei = 0;
% draw contour of cell nuclei
drawContour = 0;
% draw average contour of cell nuclei
drawAverageContour = 0;
% include contour points
includeContourPoints = 1;
% index for color of average stuff
averageColorIndex = 6;
% render average lineage lines
drawAverageLines = 1;
% render average nuclei positions
drawAverageNuclei = 1;
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
conOffset = 10;
% distance from contour point to nearest nuclei position
cellDist = 25;
% generate contour points automatically or use the points set manually
autoContour = 0;

yScale = 2.;

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

if drawContour == 1 || drawAverageContour == 1
  CONTOUR = [];
end
if drawNuclei == 1
  % ellipse instance
  ELLIP = [];
  ELLIPPATCH = [];
end

if includeContourPoints == 1
  AV = [];
  AVPATCH = [];
end

if drawAverageLines == 1
  LINKING = [];
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

[ rows, columns ] = generate2DGrid( [totalMinAxes(1) totalMinAxes(2)*yScale], [totalMaxAxes(1) totalMaxAxes(2)*yScale], resGrid );

if drawAverageLines == 1 || drawAverageContour == 1 || drawAverageDelaunay == 1
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

% number of cells in the last registered step
maxNumCellsAtLastStep = 0;

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
  
  if includeContourPoints == 1
    if ishandle( AV )
      set( AV, 'Visible', 'off' );
    end
    if ishandle( AVPATCH )
      set( AVPATCH, 'Visible', 'off' );
    end
  end
  
  if drawAverageLines == 1
    if ishandle( LINKING )
      set( LINKING, 'Visible', 'off' );
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
  lineagePos = cell( numData );
  
  % store the number of total cells at reg step
  numTotalOfCurCells = 0;
  
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
    
    % required vectors
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
        
        if cellDatas{ dataIndex }{ j, 5 } == curT
          % get position of current cell
          p = [ cellDatas{ dataIndex }{ j, 2 }...
            cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
          p = p - centerPosPerTimeStep{dataIndex}(curT,:);
          p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile, curI );
          curPos = [ curPos ; p ];
          numTotalOfCurCells = numTotalOfCurCells +1;
          
          if drawNuclei == 1
            hold on
            if overlapping == 1
              zOffset = 0.5;
            else
              zOffset = 0.;
            end
            [ ELLIP(numTotalCells), ELLIPPATCH(numTotalCells) ] =...
              drawEllipse3d( p(1), p(2)*yScale, p(3)+dataIndex*zOffset+0.2, radEllip, radEllip, 0, 0 );
            set( ELLIP(numTotalCells), 'color', colors( dataIndex, : ), 'LineWidth', lineWidth );
            set( ELLIPPATCH(numTotalCells), 'FaceColor', colors( dataIndex, : ), 'FaceLighting', 'none' );
            numTotalCells = numTotalCells+1;
          end
        end
        
        cellIDs = [cellIDs; cellDatas{ dataIndex }{ j, 1 }];
      end
    end
    % not that elegenat but necessary; go over the whole loop again to
    % match successor cells based on lineage information
    nextPos = [];
    daughterPos = [];
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
        p = p - centerPosPerTimeStep{dataIndex}(nextT,:);
        p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile, curI+1 );
        
        % ID of precursor
        cellID = cellDatas{ dataIndex }{ j, 1 };
        %precurID = cellDatas{ dataIndex }{ j, 8 };
        %T = cellDatas{ dataIndex }{ j, 5 }
        [ precurIDList, numEntries ] = getPrecursorIDList( cellDatas{ dataIndex }{ j, 8 } );
        
        % if the cell is a daughter cell of a non-division cell
        % then the pos of the next time step is located at the same index
        % position as in the previous time step
        if curT ~= nextT
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
          % if the time steps are identical then just copy the same
          % position
        else
          nextPos = [ nextPos ; p ];
        end
      end
    end
    
    % compute the mid for all daughter cells that followed quite after a
    % division such that nextPos and curPos have the same dimension for
    % applying the triangulation in t to the point set in t+1
    if size( daughterPos, 1 ) > 0
      sortrows( daughterPos, 1 );
    end
    d = 1;
    while d<size( daughterPos, 1 )
      index = daughterPos( d, 1 );
      pos1 = daughterPos( d, 2:4 );
      pos2 = daughterPos( d+1, 2:4 );
      newPos = (pos1 + pos2)/2;
      nextPos( index, :) = newPos;
      d = d + 2;
    end
    
    numCells = size( curPos, 1 );
    
    for c=1:numCells
      % determine tileIndex of mid position for a cell between two
      % subsequent time steps and add this position to the
      % tile grid in order to average all lineage lines later
      midPos = (curPos( c, : ) + nextPos( c, : ))/2.;
      tileIndex = getTileIndex( midPos, [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns );
      tileGrid{ tileIndex } = [ tileGrid{ tileIndex }; [ curPos( c, 1 ) curPos( c, 2 ) curPos( c, 3 )...
        nextPos( c, 1 ) nextPos( c, 2 ) nextPos( c, 3 ) ] ];
    end
    
    % store the total number of cells
    allCurCells( dataIndex, 1 ) = numCells;
    
    % if at least three cells exists
    if numCells > 3 && drawDelaunay == 1
      if triangulationType == 1
        % delaunay triangulation
        tri = delaunayTriangulation( curPos(:,1), curPos(:,2), curPos(:,3) );
        tetramesh( tri, 'FaceColor', colors( dataIndex, : ), 'FaceAlpha', 0.9 );
      else
        % alpha shape triangulation
        [Vol,Shape] = alphavol( [ curPos(:,1) curPos(:,2) curPos(:,3) ], sqrt( alphaRadiiVector( curT, 1 )) );
        tri = Shape.tri;
        trisurf( tri, curPos(:,1), curPos(:,2), curPos(:,3),...
          'FaceColor', colors( dataIndex, : ), 'FaceAlpha', 0. );
      end
    end
    
    % draw principal components
    if renderPrincipalComponents == 1
      start = mean( curPos );
      %coeffMat = pca( centerEllipse )
      %start = applyTransformations( start, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile, curI );
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
      dimL = size( curPos, 1 );
      if dimL > 2
        if triangulationType == 1
          K = convhull( curPos( :, 1 ), curPos( :, 2 ) );
          CONTOUR(curI) = line( curPos( K, 1 ), curPos( K, 2 ), curPos( K, 3 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
        else
          [VolC,ShapeC] = alphavol( [ curPos( :, 1 ), curPos( :, 2 ) ], sqrt( alphaRadiiVector( curT, 1 )) );
          K = ShapeC.bnd(:,1);
          dimK = size( K, 1 );
          if dimK > 1
            K(dimK+1,:) = K(1,:);
            CONTOUR(curI) = line( curPos( K, 1 ), curPos( K, 2 ), curPos( K, 3 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
          end
        end
      end
    end
  end
  
  % compute the average of lineage line positions in each tile
  if overlapping == 1
    zOffset = 0.5;
  else
    zOffset = 0.;
  end
  averagePos = [];
  nextAveragePos = [];
  numCellCounter = 0;
  fillRest = 0;
  drawingCounter = 1;
  for t=1:rows*columns
    hold on
    numMidPos = size( tileGrid{ t }, 1 );
    if numMidPos ~= 0
      if fillRest == 1
        % add the remaining line positions
        for i=1:numMidPos
          posFirst = [ tileGrid{ t }(i, 1) tileGrid{ t }(i, 2) tileGrid{ t }(i, 3)];
          posLast = [ tileGrid{ t }(i, 4) tileGrid{ t }(i, 5) tileGrid{ t }(i, 6)];
          averagePos = [ averagePos ; posFirst(1), posFirst(2)*yScale 0 ];
          nextAveragePos = [ nextAveragePos ; posLast(1), posLast(2)*yScale 0 ];
          
          linePos = [  tileGrid{ t }(i, 1) tileGrid{ t }(i, 2)*yScale ; tileGrid{ t }(i, 4) tileGrid{ t }(i, 5)*yScale ];
          LINKING(drawingCounter) = line( linePos( :, 1 ), linePos( :, 2 ), 'Color', colors( averageColorIndex+1, : ), 'LineWidth', lineWidth );
          drawingCounter = drawingCounter + 1;
        end
      else
        % compute the average of all mid line positions in the tile
        sumFirst = zeros( 1, 3 );
        sumLast = zeros( 1, 3 );
        numMidPos = size( tileGrid{ t }, 1 );
        counter = 0;
        for i=1:numMidPos
          sumFirst( 1, : ) = sumFirst( 1, : ) + [ tileGrid{ t }(i, 1) tileGrid{ t }(i, 2) tileGrid{ t }(i, 3)];
          sumLast( 1, : ) = sumLast( 1, : ) + [ tileGrid{ t }(i, 4) tileGrid{ t }(i, 5) tileGrid{ t }(i, 6)];
          
          linePos = [  tileGrid{ t }(i, 1) tileGrid{ t }(i, 2)*yScale ; tileGrid{ t }(i, 4) tileGrid{ t }(i, 5)*yScale ];
          LINKING(drawingCounter) = line( linePos( :, 1 ), linePos( :, 2 ), 'Color', colors( averageColorIndex+1, : ), 'LineWidth', lineWidth );
          drawingCounter = drawingCounter + 1;
          
          counter = counter + 1;
          numTotalOfCurCells = numTotalOfCurCells - 1;
          if numCellCounter + numTotalOfCurCells == maxNumCellsAtLastStep && curI ~= startI
            posFirst = sumFirst( 1, : )/counter;
            posLast = sumLast( 1, : )/counter;
            averagePos = [ averagePos ; posFirst(1), posFirst(2)*yScale 0 ];
            nextAveragePos = [ nextAveragePos ; posLast(1), posLast(2)*yScale 0 ];
            fillRest = 1;
            break;
          end
        end
        if fillRest ~= 1
          posFirst = sumFirst( 1, : )/numMidPos;
          posLast = sumLast( 1, : )/numMidPos;
          averagePos = [ averagePos ; posFirst(1), posFirst(2)*yScale 0 ];
          nextAveragePos = [ nextAveragePos ; posLast(1), posLast(2)*yScale 0 ];
          numCellCounter = numCellCounter + 1;
        end
      end
    end
  end

  numAverageCells = size( averagePos, 1 );
  numNextAverageCells = size( nextAveragePos, 1 );
  
%   if drawAverageLines == 1
%     for l=1:numAverageCells
%       linePos = [ averagePos(l, 1) averagePos(l, 2) ; nextAveragePos(l, 1) nextAveragePos(l, 2) ];
%       LINKING(l) = line( linePos( :, 1 ), linePos( :, 2 ), 'Color', colors( averageColorIndex+1, : ), 'LineWidth', lineWidth );
%     end
%   end
  
  % update number of required averaged cells for the next time step
  maxNumCellsAtLastStep = numAverageCells;
  
  % draw contour of data and the single marks
  if includeContourPoints == 1 &&...
      (drawAverageLines == 1 || drawAverageContour == 1 || drawAverageDelaunay == 1)
    if autoContour == 0
      factor = curI/maxI;
      cPoints = (1-factor) * cPointsFirst + factor * cPointsLast;
    else
      cPoints = generateAutomaticContourPoints( averagePos, cellDist, conOffset,...
        curI == startI, 'Average' );
    end
    
    cPoints = [ cPoints(:,1) cPoints(:,2)*yScale cPoints(:,3) ];
  end
  
  % render average nuclei positions
%   if drawAverageNuclei == 1
%     for c=1:numAverageCells
%       hold on
%       [ AV(c), AVPATCH(c) ] = drawEllipse3d( averagePos(c,1), averagePos(c,2), averagePos(c,3)+zOffset+0.2,...
%         radEllip, radEllip, 0, 0 );
%       set( AV(c), 'color', colors( averageColorIndex, : ), 'LineWidth', lineWidth );
%       set( AVPATCH(c), 'FaceColor', colors( averageColorIndex, : ), 'FaceLighting', 'none' );
%     end
%   end
  
  % include and render contour points
  if includeContourPoints == 1
    for cc=1:size( cPoints, 1 )
      hold on
      [ AV(numAverageCells+cc), AVPATCH(numAverageCells+cc) ] = drawEllipse3d( cPoints(cc,1), cPoints(cc,2), cPoints(cc,3)+zOffset+0.2,...
        radEllip, radEllip, 0, 0 );
      set( AV(numAverageCells+cc), 'color', 'k', 'LineWidth', lineWidth );
      set( AVPATCH(numAverageCells+cc), 'FaceColor', 'k', 'FaceLighting', 'none' );
      
      averagePos = [ averagePos; cPoints( cc, : ) ];
    end
  end
  
  % generate triangulation of average point set
  if triangulationType == 1
    % delaunay triangulation
    averageTri = delaunayTriangulation( averagePos(:,1), averagePos(:,2) );
  else
    % alpha shape triangulation
    [Vol,Shape] = alphavol( [ averagePos(:,1) averagePos(:,2) ], 85 );
    averageTri = Shape.tri;
  end
  
  % if at least three cells exists
  if drawAverageDelaunay == 1 && size( averagePos, 1 ) > 2
    if triangulationType == 1
      triplot( averageTri, 'b' );
      % draw triangle labels in the center of each triangle
      %TEXT = drawTriangleLabels( averageTri );
    else
      trisurf( averageTri, averagePos(:,1), averagePos(:,2), averagePos(:,3),...
        'FaceColor', 'blue', 'FaceAlpha', 0. );
    end
  end
  
  if drawAverageContour == 1
    dimL = size( averagePos, 1 );
    if dimL > 2
      if triangulationType == 1
        K = convhull( averagePos( :, 1 ), averagePos( :, 2 ) );
        CONTOUR(curI) = line( averagePos( K, 1 ), averagePos( K, 2 ), 'Color', colors( averageColorIndex, : ), 'LineWidth', lineWidth );
      else
        [VolC,ShapeC] = alphavol( [ averagePos( :, 1 ), averagePos( :, 2 ) ], 85 );
        K = ShapeC.bnd(:,1);
        dimK = size( K, 1 );
        if dimK > 1
          K(dimK+1,:) = K(1,:);
          CONTOUR(curI) = line( averagePos( K, 1 ), averagePos( K, 2 ), 'Color', colors( averageColorIndex, : ), 'LineWidth', lineWidth );
        end
      end
    end
  end
  
  if includeContourPoints == 1
    % export triangulation properties
    exportTriangulation( averageTri, averagePos, curI, 'Average', triangulationType, autoContour );
    
    if curI ~= maxI
      if autoContour == 0
        fac = (curI+1)/maxI;
        nextInterPoints = (1-fac) * cPointsFirst + fac * cPointsLast;
        nextInterPoints = [ nextInterPoints(:,1) nextInterPoints(:,2)*yScale nextInterPoints(:,3) ];
        exportNewPosOfTriangulation( nextInterPoints, nextAveragePos, 'Average', autoContour );
      else
        newCPoints = generateAutomaticContourPoints( nextAveragePos, cellDist, conOffset,...
          false, 'Average' );
        exportNewPosOfTriangulation( newCPoints, nextAveragePos, 'Average', autoContour );
      end
    end
  end
  
  hold off;
  set( f,'nextplot','replacechildren' );
  viewOffset = 100;
  axis( [ totalMinAxes(1) totalMaxAxes(1)...
    totalMinAxes(2)*yScale totalMaxAxes(2)*yScale...
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

