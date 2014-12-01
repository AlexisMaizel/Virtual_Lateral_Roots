geomPath = strcat( pwd, '/geom3d' );
addpath( geomPath );

setenv('LC_ALL','C')

%%%%% setting of properties %%%%%%
% color range for different data sets
% { 'm' 'y' 'r' 'g' 'b' 'c' };
colors = [ [ 1 0 1 ]; [ 1 1 0 ]; [ 1 0 0 ]; [ 0 1 0 ]; [ 0 0 1 ]; [ 0 1 1 ] ];
% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
cView = 2;
% startIndex
startI = 1;
% endIndex
endI = 2;
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
% index for color of average stuff
averageColorIndex = 6;
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
% offset for cell ranges
epsilon = 3;
% use triangulation based on delaunay or alpha shape
% 1 -> delaunay
% 2 -> alpha shape
triangulationType = 1;
% vector of data strings
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' '131203_raw' };
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
  prepareData( dataStr, startData, endData, numData, 'Ellipses', renderMasterFile, cView );

if drawContour == 1 || drawAverageContour == 1
  CONTOUR = [];
end
if drawNuclei == 1
  % ellipse instance
  ELLIP = [];
  ELLIPPATCH = [];
end

if drawAverageNuclei == 1
  AV = [];
  AVPATCH = [];
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

if drawAverageNuclei == 1
  if autoContour == 1
    fileName = strcat( '/tmp/triangulation-', 'Average', '_auto.txt' );
  else
    fileName = strcat( '/tmp/triangulation-', 'Average', '.txt' );
  end
  fileId = fopen( char(fileName), 'w' );
  % first write the maximum number of time steps
  fprintf( fileId, '%1d\n', startI );
  fprintf( fileId, '%1d\n', maxI );
  
  % contour points for the first and last time step
  if autoContour == 0
    cPointsFirst = generateContourPoints( 'Average', true, conOffset );
    cPointsLast = generateContourPoints( 'Average', false, conOffset );
  end
end

% number of total cells per step
numTotalCells = 1;

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
  
  if drawAverageNuclei == 1
    if ishandle( AV )
      set( AV, 'Visible', 'off' );
    end
    if ishandle( AVPATCH )
      set( AVPATCH, 'Visible', 'off' );
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
      TF = TF * createRotationOz( -pi/36. );
    end
    
    % get the alpha shape radii for all time steps
    if triangulationType == 2
      alphaRadiiVector = getAlphaRadius( dataStr( 1, dataIndex ) );
    end
    
    % only continue with this time step that is synchronized in
    % number of cells for each data set
    numNormCells = getNormalizedCellNumber( curI, 18, 143, minI, maxI );
    curT = 1;
    for j=1:maxT(dataIndex)
      if numNormCells - epsilon < numCellsPerTimeStep{dataIndex}(j,1)...
          && numNormCells + epsilon > numCellsPerTimeStep{dataIndex}(j,1)
        curT = j;
        break;
      end
    end
    
    if curI ~= endI
      numNextNormCells = getNormalizedCellNumber( curI+1, 18, 143, minI, maxI );
      nextT = 1;
      for j=1:maxT(dataIndex)
        if numNextNormCells - epsilon < numCellsPerTimeStep{dataIndex}(j,1)...
            && numNextNormCells + epsilon > numCellsPerTimeStep{dataIndex}(j,1)
          nextT = j;
          break;
        end
      end
    end
    
    % matrix of positions for current time step
    matPos = [];
    
    % determine cell position and cell file information
    nc = 1;
    % cell ids
    cellIDs = [];
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
        pos = [ cellDatas{ dataIndex }{ j, 2 }...
          cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
        matPos = [ matPos; pos ];
        cellIDs = [cellIDs; cellDatas{ dataIndex }{ j, 1 }];
      end
    end
    
    numCells = size( matPos, 1 );
    matNextPos = zeros( numCells, 3 );
    daughterPos = [];
    
    if curT ~= nextT
      if curI ~= endI
        % next loop for getting the nuclei positions of the next reg step
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
            
%             if strcmp( dataStr( 1, dataIndex ), '130607_raw' ) &&...
%                 cellDatas{ dataIndex }{ j, 6 } == 6 &&...
%                 ( cellDatas{ dataIndex }{ j, 1 } == 179 ||...
%                 cellDatas{ dataIndex }{ j, 1 } == 118 )
%               continue;
%             end
%             
%             if strcmp( dataStr( 1, dataIndex ), '121211_raw' ) &&...
%                 cellDatas{ dataIndex }{ j, 6 } == 13 &&...
%                 cellDatas{ dataIndex }{ j, 1 } == 413
%               continue;
%             end
          end
          % also check next time step
          if cellDatas{ dataIndex }{ j, 5 } == nextT
            pos = [ cellDatas{ dataIndex }{ j, 2 }...
              cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
            
            % ID of precursor
            cellID = cellDatas{ dataIndex }{ j, 1 };
            precurID = cellDatas{ dataIndex }{ j, 8 };
            
            % if the cell is a daughter cell of a non-division cell
            % then the pos of the next time step is located at the same index
            % position as in the previous time step
            if precurID == -1
              [row,col] = find( cellIDs == cellID );
              matNextPos( row, :) = pos;
              % else there is a division and we store the both positions of
              % the two daughter cells to later generate a cell that is associated with
              % the cell at the previous time step
            else
              [row,col] = find( cellIDs == precurID );
              daughterPos = [ daughterPos; row pos ];
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
          matNextPos( index, :) = newPos;
          d = d + 2;
        end
        % else both reg. steps correspond to the same time step
      end
    else
      matNextPos = matPos;
    end
    
    % check if the current time step should be skipped depending
    % on the number of cells
    %       if renderCellRanges == 1
    %         skipTimeStep = 1;
    %         for nc=1:size( cellRange, 2 )
    %           if numCells - epsilon < cellRange(nc)...
    %               && numCells + epsilon > cellRange(nc)
    %             skipTimeStep = 0;
    %             break;
    %           end
    %         end
    %         if skipTimeStep == 1
    %           continue;
    %         end
    %       end
    
    centerEllipse = zeros( numCells, 3 );
    % draw an ellipsoid for each cell
    for c=1:numCells
      hold on
      % get position of current cell
      p = [ matPos( c, 1 ) matPos( c, 2 ) matPos( c, 3 ) ];
      p = p - centerPosPerTimeStep{dataIndex}(curT,:);
      p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile, curI );
      
      % get position of next cell
      if curI ~= endI
        np = [ matNextPos( c, 1 ) matNextPos( c, 2 ) matNextPos( c, 3 ) ];
        np = np - centerPosPerTimeStep{dataIndex}(nextT,:);
        np = applyTransformations( np, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile, curI+1 );
      end
      
      centerEllipse(c, :) = p;
      
      % determine tileIndex of current nuclei position and add it to the
      % tile grid in order to average all positions later
      tileIndex = getTileIndex( p, [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns );
      tileGrid{ tileIndex } = [ tileGrid{ tileIndex }; [ p(1) p(2) p(3) ] ];
      if curI ~= endI
        nextTileIndex = getTileIndex( np, [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns );
        nextTileGrid{ nextTileIndex } = [ tileGrid{ nextTileIndex }; [ np(1) np(2) np(3) ] ];
      end
      
      if drawNuclei == 1
        if overlapping == 1
          zOffset = 0.5;
        else
          zOffset = 0.;
        end
        [ ELLIP(numTotalCells) ELLIPPATCH(numTotalCells) ] =...
          drawEllipse3d( p(1), p(2), p(3)+dataIndex*zOffset+0.2, radEllip, radEllip, 0, 0 );
        set( ELLIP(numTotalCells), 'color', colors( dataIndex, : ), 'LineWidth', lineWidth );
        set( ELLIPPATCH(numTotalCells), 'FaceColor', colors( dataIndex, : ), 'FaceLighting', 'none' );
        numTotalCells = numTotalCells+1;
      end
    end
    
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
  if drawAverageNuclei == 1
    if overlapping == 1
      zOffset = 0.5;
    else
      zOffset = 0.;
    end
    averagePos = [];
    nextAveragePos = [];
    for t=1:rows*columns
      hold on
      if size( tileGrid{ t }, 1 ) ~= 0
        % compute the average of all positions in the tile
        sum = 0;
        for i=1:size( tileGrid{ t }, 1 )
          sum = sum + [ tileGrid{ t }(i, 1) tileGrid{ t }(i, 2) tileGrid{ t }(i, 3)];
        end
        sum = sum/size( tileGrid{ t }, 1 );
        averagePos = [ averagePos ; sum(1), sum(2) 0 ];
        [ AV(t) AVPATCH(t) ] = drawEllipse3d( sum(1), sum(2), sum(3)+zOffset+0.2, radEllip, radEllip, 0, 0 );
        set( AV(t), 'color', colors( averageColorIndex, : ), 'LineWidth', lineWidth );
        set( AVPATCH(t), 'FaceColor', colors( averageColorIndex, : ), 'FaceLighting', 'none' );
      end
      
      % next step
      if curI ~= endI
        if size( nextTileGrid{ t }, 1 ) ~= 0
          % compute the average of all positions in the tile
          sum = 0;
          for i=1:size( nextTileGrid{ t }, 1 )
            sum = sum + [ nextTileGrid{ t }(i, 1) nextTileGrid{ t }(i, 2) nextTileGrid{ t }(i, 3)];
          end
          sum = sum/size( nextTileGrid{ t }, 1 );
          nextAveragePos = [ nextAveragePos ; sum(1), sum(2) 0 ];
        end
      end
    end
    
    % draw contour of data and the single marks
    if autoContour == 0
      factor = curI/maxI;
      cPoints = (1-factor) * cPointsFirst + factor * cPointsLast;
    else
      cPoints = generateAutomaticContourPoints( averagePos, cellDist, conOffset,...
        curI == startI, 'Average' );
    end
    
    % include and render contour points
    radii = 5;
    for cc=1:size( cPoints, 1 )
      hold on
      [ AV(rows*columns+cc) AVPATCH(rows*columns+cc) ] = drawEllipse3d( cPoints(cc,1), cPoints(cc,2), cPoints(cc,3)+zOffset+0.2,...
        radEllip, radEllip, 0, 0 );
      set( AV(rows*columns+cc), 'color', 'k', 'LineWidth', lineWidth );
      set( AVPATCH(rows*columns+cc), 'FaceColor', 'k', 'FaceLighting', 'none' );
    end
    
    % add the contour points to generate a complete triangulation
    for cc=1:size( cPoints, 1 )
      averagePos = [ averagePos; cPoints( cc, : ) ];
    end
    
    averagePos
    
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
    
    % export triangulation properties
    exportTriangulation( averageTri, averagePos, curI, 'Average', triangulationType, autoContour );
    
    if curI ~= endI
      if autoContour == 0
        fac = (curI+1)/maxI;
        nextInterPoints = (1-fac) * cPointsFirst + fac * cPointsLast;
        nextInterPoints
        nextAveragePos
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
    totalMinAxes(2) totalMaxAxes(2)...
    totalMinAxes(3)-viewOffset totalMaxAxes(3)+viewOffset ...
    totalMinAxes(3) totalMaxAxes(3) ] );
  axis on
  daspect( [ 1 1 1 ] );
  
  % legend
  leg = legend( '120830', '121204', '121211', '130508', '130607', 'Average' );
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
ElapsedTimeIndexLoop = cputime - cpuT
