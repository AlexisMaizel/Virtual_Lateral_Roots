geomPath = strcat( pwd, '/geom3d' );
addpath( geomPath );

%%%%% setting of properties %%%%%%
% color range for different data sets
% { 'm' 'k' 'r' 'g' 'b' 'c' };
colors = [ [ 1 0 1 ]; [ 0 0 0 ]; [ 1 0 0 ]; [ 0 1 0 ]; [ 0 0 1 ] ];
% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
cView = 2;
% startIndex
startI = 20;
% endIndex
endI = 20;
% min and max index
minI = 1;
maxI = 20;
% draw delaunay tri?
drawDelaunay = 0;
% draw point set as ellipses
drawNuclei = 0;
% draw contour of cell nuclei
drawContour = 1;
% if set to one then only a single time step
% is rendered given by startT
exportType = 2;
% vector of data strings
exportTypeStr = { 'SingleFigure' 'AsImages' };
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
  prepareData( dataStr, startData, endData, numData, visualizationType( 1, visType ), renderMasterFile, cView );

if drawContour == 1
  CONTOUR = [];
end
if drawNuclei == 1
  % ellipse instance
  ELLIP = [];
  ELLIPPATCH = [];
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

hold on;
[ rows, columns ] = generate2DGrid( [totalMinAxes(1) totalMinAxes(2)], [totalMaxAxes(1) totalMaxAxes(2)], resGrid );

cpuT = cputime;
% loop over all normalized steps
for curI=startI:endI
  
  if drawContour == 1
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
    
    % number of cells for current time step and data set
    numCells = numCellsPerTimeStep{dataIndex}(curT, 1);
    % matrix of positions for current time step
    matPos = zeros( numCells, 3 );
    % vector of object ids in order to access the cell
    % file information in the cell file map
    idToCF = zeros( numCells, 3 );
    
    % determine cell position and cell file information
    nc = 1;
    for j=1:dimData( dataIndex )
      if cellDatas{ dataIndex }{ j, 5 } == curT
        pos = [ cellDatas{ dataIndex }{ j, 2 } cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
        matPos(nc, :) = pos;
        idToCF(nc, :) = cellDatas{ dataIndex }{ j, 1 };
        nc = nc + 1;
      end
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
    
    numConsideredCells = 0;
    for c=1:numCells
      % only consider the master cell file if required
      if renderMasterFile == 1 &&...
          0 ~= cellFileMap{ dataIndex }( idToCF(c,1) )
        continue;
      end
      numConsideredCells = numConsideredCells + 1;
    end
    
    centerEllipse = zeros( numConsideredCells, 3 );
    nc = 1;
    % draw an ellipsoid for each cell
    for c=1:numCells
      % only consider the master cell file if required
      if renderMasterFile == 1 &&...
          0 ~= cellFileMap{ dataIndex }( idToCF(c,1) )
        continue;
      end
      
      % get position of current cell
      p = [ matPos( c, 1 ) matPos( c, 2 ) matPos( c, 3 ) ];
      p = p - centerPosPerTimeStep{dataIndex}(curT,:);
      p = applyTransformations( p, planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile, curI );
      centerEllipse(nc, :) = p;
      
      if drawNuclei == 1
        if overlapping == 1
          zOffset = 0.5;
        else
          zOffset = 0.;
        end
        hold on;
        [ ELLIP(dataIndex*numConsideredCells + nc) ELLIPPATCH(dataIndex*numConsideredCells + nc) ] =...
          drawEllipse3d( p(1), p(2), p(3)+dataIndex*zOffset+0.2, radEllip, radEllip, 0, 0 );
        set( ELLIP(dataIndex*numConsideredCells + nc), 'color', colors( dataIndex, : ), 'LineWidth', lineWidth );
        set( ELLIPPATCH(dataIndex*numConsideredCells + nc), 'FaceColor', colors( dataIndex, : ), 'FaceLighting', 'none' );
      end
      
      nc = nc + 1;
    end
    
    % if at least three cells exists
    if numCells > 3 && drawDelaunay == 1
      hold on;
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
        endP = applyTransformations( coeff(:, a), planePos, u, v, TF, dataStr( 1, dataIndex ), renderMasterFile, curI );
        hold on;
        P(a) = quiver3( start(1), start(2), start(3),...
          endP(1), endP(2), endP(3),...
          arrowLength, 'LineWidth', 3,...
          'Color', cm( a, : ) );
      end
    end
    
    if drawContour == 1
      hold on;
      dimL = size( centerEllipse, 1 );
      if dimL > 2
        if triangulationType == 1
          K = convhull( centerEllipse( :, 1 ), centerEllipse( :, 2 ) );
          CONTOUR(dataIndex) = line( centerEllipse( K, 1 ), centerEllipse( K, 2 ), centerEllipse( K, 3 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
        else
          [VolC,ShapeC] = alphavol( [ centerEllipse( :, 1 ), centerEllipse( :, 2 ) ], sqrt( alphaRadiiVector( curT, 1 )) );
          K = ShapeC.bnd(:,1);
          dimK = size( K, 1 );
          if dimK > 1
            K(dimK+1,:) = K(1,:);
            CONTOUR(dataIndex) = line( centerEllipse( K, 1 ), centerEllipse( K, 2 ), centerEllipse( K, 3 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
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
  
  % legend
  hold on;
  leg = legend( '120830', '121204', '121211', '130508', '130607' );
  set(leg, 'Location', 'NorthWestOutside');
  linecolors = { [ 1 0 1 ], [ 0 0 0 ], [ 1 0 0 ], [ 0 1 0 ], [ 0 0 1 ] };
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
