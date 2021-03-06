geomPath = strcat( pwd, '/geom3d' );
addpath( geomPath );

setenv('LC_ALL','C')

%%%%% setting of properties %%%%%%
% data Index:
% 1 -> 120830
% 2 -> 121204
% 3 -> 121211
% 4 -> 130508
% 5 -> 130607
% 6 -> 131203
% 7 -> all
dataId = 1;
% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
cView = 2;
% start with the current time step
startT = 1;
maxT = 1;
maxDataT = 300;
% draw delaunay tri?
drawDelaunay = 1;
% render only master file?
renderSingleCellFile = 1;
% render contour (=convexHull) instead of ellipses
renderContour = 0;
% include contour points
includeContourPoints = 1;
% use triangulation based on delaunay or alpha shape
% 1 -> delaunay
% 2 -> alpha shape
triangulationType = 1;
% vector of data strings
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' '131203_raw' };
% vector of view strings
viewStr = { 'Top' 'Side' 'Radial' '3D' };
% either loop over all data sets creating five figures or only
% show one specific one
startD = 1;
endD = 6;
if dataId ~= 7
  startD = dataId;
  endD = dataId;
end

% define offset for boundary points which is the distance between the
% acutal location of the contour points and their initial position within
% the model
conOffset = 5;

% distance from contour point to nearest nuclei position
cellDist = 25;

% generate contour points automatically or use the points set manually
autoContour = 0;

color = [ 0 1 0 ];

% path to image output
imageDir = strcat( 'images/AllData/' );

% create directory if required
mkdir( char(imageDir) );

% loop over all data sets
for dataIndex=startD:endD
  % output format of values
  format longG
  % reading raw data
  path = strcat( '../FinalVLRForMatlab/', dataStr( 1, dataIndex ), '.csv' );
  fileID = fopen( char(path) );
  % format of data sets:
  % ObjectID X Y Z Timepoint Radius Precursors Color Lineage TrackID TrackColor TrackGroup Layer DivisionType
  formatSpec = '%d %f %f %f %d %d %q %q %d %q %q %d %d %q';
  if strcmp( dataStr( 1, dataIndex ), '131203_raw' )
    formatSpec = '%d %f %f %f %d %d %q %q %d %d %d %q';
  end
  % read data and ignore the first four header lines
  data = textscan( fileID, formatSpec, 'HeaderLines', 4, 'Delimiter', ';' );
  fclose(fileID);
  
  % get dimension aka number of lines
  col = size(data{1});
  numLines = col(1,1);
  % store relevant columns in corresponding data structures
  % and ignore the others
  % Object Id
  IdCol = data{1};
  % X coord
  XCol = data{2};
  % Y coord
  YCol = data{3};
  % z coord; Note that due to resampling, the z value is multiplied with 2
  if strcmp( dataStr( 1, dataIndex ), '121204_raw_2014' )
    ZCol = -2 * data{4};
  else
    ZCol = 2 * data{4};
  end
  % Time Step
  TCol = data{5};
  % precursors of cells
  PCol = data{7};
  % Lineage Tree Id
  LCol = data{9};
  % temporal for data 131203
  if strcmp( dataStr( 1, dataIndex ), '131203_raw' )
    % Cell File Id
    CFCol = data{10};
    % Cell Layer Id
    CLCol = data{11};
    % Divison Type
    DCol = data{12};
  else
    % Cell File Id
    CFCol = data{12};
    % Cell Layer Id
    CLCol = data{13};
    % Divison Type
    %DCol = data{14};
  end
  
  l = 1;
  % first determine dimension of cellData
  % in order to initialize the size of cellData
  % -> huge performance speed up
  dim = 0;
  while (l < numLines+1)
    firstCellId = IdCol(l);
    secondCellId = IdCol(l+1);
    if firstCellId == secondCellId
      firstTS = TCol(l);
      secondTS = TCol(l+1);
      deltaTS = secondTS - firstTS;
      dim = dim + deltaTS + 1;
      l = l+2;
    else
      dim = dim+1;
      l = l+1;
    end
  end
  
  % cellData is the main array with all relevant information
  % for the further analysis:
  % ObjectID X Y Z Timepoint LineageID TrackGroup Last precursor
  cellData = cell( dim, 8 );
  cellFileMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
  
  l = 1;
  nl = 1;
  % interpolate the missing positions in between
  while (l < numLines+1)
    firstCellId = IdCol(l);
    secondCellId = IdCol(l+1);
    
    % check last precursor to detect daughter cells and the id of
    % their last parent
    % if the id is -1 than the cell starts at the root which we use as an
    % identification of non daughter cells
    % we set this value for all points interpolated in between
    lastPrecur = getLastPrecursorID( PCol(l) );
    
    % insert first line
    cellData( nl, : ) = {firstCellId XCol(l) YCol(l) ZCol(l) TCol(l) LCol(l) CFCol(l) lastPrecur};
    nl = nl+1;
    
    % interpolate between cell positions
    if firstCellId == secondCellId
      firstTS = TCol(l);
      secondTS = TCol(l+1);
      deltaTS = secondTS - firstTS;
      steps = double( 1./double(deltaTS) );
      
      for s=1:deltaTS-1
        t = firstTS + s;
        k = double( double(s)*steps );
        x = double(1-k) * XCol(l) + double(k) * XCol(l+1);
        y = double(1-k) * YCol(l) + double(k) * YCol(l+1);
        z = double(1-k) * ZCol(l) + double(k) * ZCol(l+1);
        % insert all relevant data into main data structure
        cellData( nl, : ) = {firstCellId x y z t LCol(l) CFCol(l) -1};
        nl = nl+1;
      end
      
      % insert last line
      cellData( nl, : ) = {firstCellId XCol(l+1) YCol(l+1) ZCol(l+1) TCol(l+1) LCol(l+1) CFCol(l+1) -1};
      nl = nl+1;
      % update cell file map
      cellFileMap( firstCellId ) = CFCol(l);
      
      % increment loop index
      l = l+2;
      % else cell exists only for one time step
    else
      % update cell file map
      cellFileMap( firstCellId ) = CFCol(l);
      
      l = l+1;
    end
  end
  
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
  C = strsplit( char( dataStr( 1, dataIndex ) ), '_' );
  title( strcat( C( 1, 1 ), ' Time Step ', num2str(startT) ) );
  camproj( 'orthographic' )
  
  % gca is the current axes handle
  set( gca,'nextplot','replacechildren' );
  % gcf is the current figure handle
  %lighting phong
  set( gcf, 'Renderer', 'zbuffer' );
  lighting gouraud
  set( gcf, 'Renderer', 'OpenGL' );
  set( gcf,'nextplot','replacechildren' );
  set( gcf, 'color', [ 1 1 1 ] );
  
  % surface instance
  S = [];
  % PC instance
  P = [];
  DATA = [];
  DATAF = [];
  DATAL = [];
  CSF = [];
  CSL = [];
  TEXT = [];
  
  % get stored eigenvectors for the last time step to set the same
  % direction view for each time step
  coeff = getPrincipalComponents( dataStr( 1, dataIndex ), renderSingleCellFile );
  
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
      v = -coeff(:,2);
    end
  elseif cView == 3
    dir = -coeff(:,1);
    u = -coeff(:,3);
    v = coeff(:,2);
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
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Traversal over all time steps %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if autoContour == 1
    fileName = strcat( '/tmp/triangulation-', dataStr( 1, dataIndex ), '_auto.txt' );
  else
    fileName = strcat( '/tmp/triangulation-', dataStr( 1, dataIndex ), '.txt' );
  end
  fileId = fopen( char(fileName), 'w' );
  % first write the maximum number of time steps
  %fprintf( fileId, '%1d\n', startT );
  %fprintf( fileId, '%1d\n', maxT );

  % contour points for the first and last time step
  if autoContour == 0
    cPointsFirst = generateContourPoints( dataStr( 1, dataIndex ), true, conOffset );
    cPointsLast = generateContourPoints( dataStr( 1, dataIndex ), false, conOffset );
  end
  
  % loop over all time steps
  for curT=startT:maxT
    % draw contour of data and the single marks
    if autoContour == 0
      factor = curT/maxDataT;
      cPoints = (1-factor) * cPointsFirst + factor * cPointsLast;
    end
    
    % number of cells for the current time step
    numCells = 0;
    
    % matrix of positions for current time step
    curPos = [];
    
    % cell ids
    cellIDs = [];
    
    for j=1:dim
      if cellData{ j, 5 } == curT
        
        if renderSingleCellFile == 1 && cellData{ j, 7 } ~= 0
          continue;
        end
        
        % this is a special case for the data set 130508 for which we
        % ignore the two cells that arise in the master cell file with
        % lineage ID 5
        if strcmp( dataStr( 1, dataIndex ), '130508_raw' ) && cellData{ j, 6 } == 5
          continue;
        end
        
        pos = [ cellData{ j, 2 } cellData{ j, 3 } cellData{ j, 4 } ];
        pos = projectOnPlane( pos, planePos, u, v );
        pos = transformPoint3d( pos, TF );
        
        curPos = [curPos; pos];
        cellIDs = [cellIDs; cellData{ j, 1 }];
        numCells = numCells +1;
      end
    end
    % not that elegenat but necessary to go over the whole loop again
    nextPos = zeros( numCells, 3 );
    daughterPos = [];
    for j=1:dim
      % also check next time step
      if cellData{ j, 5 } == curT+1
        
        if renderSingleCellFile == 1 && cellData{ j, 7 } ~= 0
          continue;
        end
        
        % this is a special case for the data set 130508 for which we
        % ignore the two cells that arise in the master cell file with
        % lineage ID 5
        if strcmp( dataStr( 1, dataIndex ), '130508_raw' ) && cellData{ j, 6 } == 5
          continue;
        end
        
        pos = [ cellData{ j, 2 } cellData{ j, 3 } cellData{ j, 4 } ];
        pos = projectOnPlane( pos, planePos, u, v );
        pos = transformPoint3d( pos, TF );
        
        % ID of precursor
        cellID = cellData{ j, 1 };
        precurID = cellData{ j, 8 };
        
        % if the cell is a daughter cell of a non-division cell
        % then the pos of the next time step is located at the same index
        % position as in the previous time step
        if precurID == -1
          [row,col] = find( cellIDs == cellID );
          nextPos( row, :) = pos;
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
      nextPos( index, :) = newPos;
      d = d + 2;
    end
    
    % clean figure content by removing the
    % last ellipsoids and lines in the previous time step
    if ishandle( S )
      set( S, 'Visible', 'off' );
    end
    if ishandle( P )
      set( P, 'Visible', 'off' );
    end
    if ishandle( DATA )
      set( DATA, 'Visible', 'off' );
    end
    if ishandle( CSF )
      set( CSF, 'Visible', 'off' );
    end
    if ishandle( CSL )
      set( CSL, 'Visible', 'off' );
    end
    if ishandle( TEXT )
      set( TEXT, 'Visible', 'off' );
    end
    
    % generate the contour points automatically
    if autoContour == 1
      cPoints = generateAutomaticContourPoints( curPos, cellDist, conOffset,...
        curT == startT, dataStr( 1, dataIndex ) );
    end
    
    % add the contour points to generate a complete triangulation
    if includeContourPoints == 1
      for cc=1:size( cPoints, 1 )
        curPos = [ curPos; cPoints( cc, : ) ];
      end
    end
    
    % determine min and max values of x,y for view options
    minX = min( curPos(:,1) )
    minY = min( curPos(:,2) )
    maxX = max( curPos(:,1) )
    maxY = max( curPos(:,2) )
    
    if triangulationType == 1
      % delaunay triangulation
      curTri = delaunayTriangulation( curPos(:,1), curPos(:,2) );
    else
      % alpha shape triangulation
      [Vol,Shape] = alphavol( [ curPos(:,1) curPos(:,2) ], 85 );
      curTri = Shape.tri;
    end
    
    % export triangulation properties
    exportTriangulation( curTri, curPos, curT, dataStr( 1, dataIndex ), triangulationType, autoContour );
    
    if curT ~= maxT
      if includeContourPoints == 1
        if autoContour == 0
          fac = (curT+1)/maxDataT;
          nextInterPoints = (1-fac) * cPointsFirst + fac * cPointsLast;
          numConMarks = size( nextInterPoints, 1 );
          exportNewPosOfTriangulation( nextInterPoints, nextPos, dataStr( 1, dataIndex ), autoContour );
        else
          newCPoints = generateAutomaticContourPoints( nextPos, cellDist, conOffset,...
          false, dataStr( 1, dataIndex ) );
          numCPoints = size( newCPoints, 1 );
          exportNewPosOfTriangulation( newCPoints, nextPos, dataStr( 1, dataIndex ), autoContour );
        end
      else
        fileName = strcat( '/tmp/triangulation-', dataStr( 1, dataIndex ), '.txt' );
        fileId = fopen( char(fileName), 'a' );
        nextPos = nextPos';
        fprintf( fileId, '%4f %4f\n', nextPos( 1:2, : ) );
        fprintf( fileId, '\n' );
        fclose( fileId );
      end
    end
    
    % if at least three cells exists
    if drawDelaunay == 1 && size( curPos, 1 ) > 2
      if triangulationType == 1
        triplot( curTri, 'b' );
        % draw triangle labels in the center of each triangle
        %TEXT = drawTriangleLabels( curTri );
      else
        trisurf( curTri, curPos(:,1), curPos(:,2), curPos(:,3),...
                 'FaceColor', 'blue', 'FaceAlpha', 0. );
      end
    end
    
    points = zeros( numCells, 3 );
    % draw an ellipsoid for each cell
    for c=1:numCells
      % get position of current cell
      p1 = [ curPos( c, 1 ) curPos( c, 2 ) curPos( c, 3 ) ];
      
      points( c, : ) = p1;
      
      % draw the single cell as sphere
      radii = 8;
      [ X, Y, Z ] = ellipsoid( p1(1), p1(2), p1(3), radii/2., radii/2., radii/2., 20 );
      
      % render sphere surfaces
      S(c) = surface( X, Y, Z, 'FaceColor', [ 1 0 0 ], 'EdgeColor', 'none', 'FaceLighting', 'gouraud' );
    end
    
    if size( points, 1 ) > 2 && renderContour == 1
      K = convhull( points( :, 1 ), points( :, 2 ) );
      DATA(curT) = line( points( K, 1 ), points( K, 2 ), points( K, 3 ), 'Color', color, 'LineWidth', 1.2 );
    end
    
    if includeContourPoints == 1
      radii = 5;
      for cc=1:size( cPoints, 1 )
        [ Xc, Yc, Zc ] = ellipsoid( cPoints(cc,1), cPoints(cc,2), cPoints(cc,3), radii/2., radii/2., radii/2., 20 );
        CSF(cc) = surface( Xc, Yc, Zc, 'FaceColor', [ 0 0 0 ], 'EdgeColor', 'none', 'FaceLighting', 'gouraud' );
      end
    end

    hold off;
    set( f,'nextplot','replacechildren' );
    viewOffset = 10;
    minMax = getTotalMinMax( dataStr( 1, dataIndex ), autoContour );
    axis( [ minMax(1,1)-viewOffset minMax(1,2)+viewOffset minMax(1,3)-viewOffset minMax(1,4)+viewOffset ] );
    axis on
    daspect( [ 1 1 1 ] );
    
    grid off;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    C = strsplit( char( dataStr( 1, dataIndex ) ), '_' );
    title( strcat( C( 1, 1 ), ' Time Step ', num2str(curT) ) );
    
    % name and type of contour points
    if autoContour == 0
      name = strcat( C( 1, 1 ), '_' );
    else
      name = strcat( C( 1, 1 ), '_auto_' );
    end
    
    if curT < 10
      digit = strcat( viewStr( 1, cView ), '_00' );
    elseif curT < 100
      digit = strcat( viewStr( 1, cView ), '_0' );
    else
      digit = strcat( viewStr( 1, cView ), '_' );
    end
    
    % output with number of cells
    filePath = strcat( imageDir, name, digit, num2str(curT), '.png' );
    
    saveas( gcf, char(filePath) );
    
  end
  % close figure after processing
  %close(f);
end
