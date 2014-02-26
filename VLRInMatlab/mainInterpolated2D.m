%%%% The mesh deformation is based on the paper of Graner et al,
%%%% Discrete rearranging disordered patterns, part I: Robust
%%%% statistical tools in two or three dimensions, Eur. Phys. J.,
%%%% pp 349 -- 369, 2008

addpath( '/home/necrolyte/Data/VLR/Virtual_Lateral_Roots/VLRInMatlab/geom3d' );

%%%%% setting of properties %%%%%%
% data Index:
% 1 -> 120830
% 2 -> 121204
% 3 -> 121211
% 4 -> 130508
% 5 -> 130607
% 6 -> 131203
% 7 -> all
dataId = 2;
% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
cView = 2;
% start with the current time step
startT = 300;
% draw delaunay tri?
drawDelaunay = 0;
% if set to one then only a single time step
% is rendered given by startT
exportType = 2;
% vector of data strings
exportTypeStr = { 'SingleFigure' 'AsImages' 'AsVideo' };
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
renderCellRanges = 0;
% only generate output images if the number of cells are within a desired
% range
cellRange = [ 20 40 60 80 100 120 140 ];
% offset for cell ranges
epsilon = 2;
% choose which major lines of the ellipsoids should be rendered
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

% vector of data strings
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' '131203_raw' };

% vector of view strings
viewStr = { 'Top' 'Side' 'Radial' '3D' };

% master cell file information taken by the picture made of Daniel located in dropbox
% and the trackGroup information of the raw data sets
%masterCellFile = [ 4 3 4 2 3 0 ];
masterCellFile = [ 4 4 4 3 3 0 ];

% either loop over all data sets creating five figures or only
% show one specific one
startD = 1;
endD = 6;
if dataId ~= 7
  startD = dataId;
  endD = dataId;
end

% loop over all data sets
for dataIndex=startD:endD
  
  % render the corresponding master cell file
  singleCellFile = masterCellFile( 1, dataIndex );
  
  % path to movie output
  movieDir = strcat( 'videos/', dataStr( 1, dataIndex ), '_',...
    viewStr( 1, cView ), '_test_movie.avi');
  
  % create directory if required
  if strcmp( exportTypeStr( 1, exportType ), 'AsVideos' )
    mkdir( 'videos/' );
  end
  
  % path to image output
  imageDir = strcat( 'images/', dataStr( 1, dataIndex ), '/' );
  
  % create directory if required
  if strcmp( exportTypeStr( 1, exportType ), 'AsImages' )
    mkdir( char(imageDir) );
  end
  
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
  
  % get maximum of time steps
  % if a movie should be created
  if strcmp( exportTypeStr( 1, exportType ), 'AsVideo' ) ||...
     strcmp( exportTypeStr( 1, exportType ), 'AsImages' )
    maxT = max( TCol );
    % else only render current time step
  else
    maxT = startT;
  end
  
  % get bounding box are by determining
  % min and max of x/y/z values
  % loop over all time steps
  minX = min( XCol );
  minY = min( YCol );
  minZ = min( ZCol );
  maxX = max( XCol );
  maxY = max( YCol );
  maxZ = max( ZCol );
  
  % store min and max of cell files
  minCF = min( CFCol );
  maxCF = max( CFCol );
  
  % number of subdivisions for ellipsoids
  nEllip = 10;
  
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
  % ObjectID X Y Z Timepoint LineageID TrackGroup
  %cellData = double.empty( 0, 7 );
  cellData = cell( dim, 7 );
  cellFileMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
  
  l = 1;
  nl = 1;
  % interpolate the missing positions in between
  while (l < numLines+1)
    firstCellId = IdCol(l);
    secondCellId = IdCol(l+1);
    % insert first line
    cellData( nl, : ) = {firstCellId XCol(l) YCol(l) ZCol(l) TCol(l) LCol(l) CFCol(l)};
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
        cellData( nl, : ) = {firstCellId x y z t LCol(l) CFCol(l)};
        nl = nl+1;
      end
      
      % insert last line
      cellData( nl, : ) = {firstCellId XCol(l+1) YCol(l+1) ZCol(l+1) TCol(l+1) LCol(l+1) CFCol(l+1)};
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
  
  % number of cell files
  numCellFiles = double( maxCF-minCF+1 );
  ticks = [ minCF:maxCF ];
  % set colormap and colorbar depending on the number of cell files
  cm = hsv( numCellFiles );
  colormap( cm );
  %h = colorbar;
  %set( h, 'YTickMode', 'manual' );
  %set( h, 'YTickLabel', ticks );
  
  % gca is the current axes handle
  set( gca,'nextplot','replacechildren' );
  % gcf is the current figure handle
  %lighting phong
  set( gcf, 'Renderer', 'zbuffer' );
  lighting gouraud
  set( gcf, 'Renderer', 'OpenGL' );
  set( gcf,'nextplot','replacechildren' );
  set( gcf, 'color', [ 1 1 1 ] );
  
  % video output options
  if strcmp( exportTypeStr( 1, exportType ), 'AsVideo' )
    writerObj = VideoWriter( char(movieDir) );
    writerObj.FrameRate = 3;
    writerObj.Quality = 100;
    open( writerObj );
  end
  
  maxTimeStep = max( TCol );
  numCellsAtStart = 0;
  numCellsAtEnd = 0;
  % get maximum number of cells at the beginning and at the end
  for j=1:dim
    % next time step
    if cellData{ j, 5 } == 1
      numCellsAtStart = numCellsAtStart +1;
    elseif cellData{ j, 5 } == maxTimeStep
      numCellsAtEnd = numCellsAtEnd +1;
    end
  end
  
  % surface instance
  S = [];
  % semi axes instances
  MIN = [];
  MID = [];
  MAX = [];
  % PC instance
  P = [];
  % ellipse instance
  ELLIP = [];
  ELLIPPATCH = [];
  
  % get stored eigenvectors for the last time step to set the same
  % direction view for each time step
  coeff = getPrincipalComponents( dataStr( 1, dataIndex ) );
  
  % set PC depending on the viewing direction
  if cView == 1
    dir = coeff(:,2);
    u = coeff(:,1);
    v = coeff(:,3);
  elseif cView == 2
    dir = coeff(:,3);
    u = coeff(:,1);
    v = coeff(:,2);
  elseif cView == 3
    dir = coeff(:,1);
    u = coeff(:,2);
    v = coeff(:,3);
  end
  
  % set plane position
  planePos = dir * 1;
  
  plane = [ planePos(1) planePos(2) planePos(3)...
    u(1) u(2) u(3)...
    v(1) v(2) v(3) ];
  TF = createBasisTransform3d( 'g', plane );
  
  % set axes properties
  minAxes = projectOnPlane( [ minX minY minZ ], planePos, u, v );
  minAxes = transformPoint3d( minAxes, TF );
  maxAxes = projectOnPlane( [ maxX maxY maxZ ], planePos, u, v );
  maxAxes = transformPoint3d( maxAxes, TF );
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Traversal over all time steps and texture field generation %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  imgStart = 1;
  % loop over all time steps
  for curT=startT:maxT
    % number of cells for the current time step
    numCells = 0;
    
    % matrix of positions for current time step
    matPos = [];
    
    % vector of object ids in order to access the cell
    % file information in the cell file map
    idToCF = [];
    
    for j=1:dim
      if cellData{ j, 5 } == curT
        pos = [ cellData{ j, 2 } cellData{ j, 3 } cellData{ j, 4 } ];
        matPos = [matPos; pos];
        numCells = numCells +1;
        idToCF = [ idToCF ; cellData{ j, 1 } ];
      end
    end
    
    % clean figure content by removing the
    % last ellipsoids and lines in the previous time step
    set( S, 'Visible', 'off' );
    set( MAX, 'Visible', 'off' );
    set( MID, 'Visible', 'off' );
    set( MIN, 'Visible', 'off' );
    set( ELLIP, 'Visible', 'off' );
    set( ELLIPPATCH, 'Visible', 'off' );
    set( P, 'Visible', 'off' );
    
    % if at least three cells exists
    if numCells > 3
      % delaunay triangulation visualzation
      tri = delaunayTriangulation( matPos(:,1), matPos(:,2), matPos(:,3) );
      
      % number of total links in current time step
      %numTotalLinks = size( edges(tri), 1 );
      
      if drawDelaunay == 1
        tetramesh(tri, 'FaceColor', cm( 1, : ), 'FaceAlpha', 0.9 );
      end
      
      cellFileMat = [];    
      linePos = [];
      minMaxSemiAxisVector = [];
      centerEllipse = [];
      
      % check if the current time step should be skipped depending
      % on the number of cells
      if renderCellRanges == 1
        skipTimeStep = 1;
        for nc=1:size( cellRange, 2 )
          if numCells - epsilon < cellRange(nc)...
              && numCells + epsilon > cellRange(nc)
            skipTimeStep = 0;
            break;
          end
        end
        if skipTimeStep == 1
          continue;
        end
      end
      
      % draw an ellipsoid for each cell
      for c=1:numCells
        % get position of current cell
        p1 = [ tri.Points( c, 1 ) tri.Points( c, 2 ) tri.Points( c, 3 ) ];
        
        % get color for ellipsoid
        if renderSingleCellFile == 0
          color = cm( cellFileMap( idToCF(c,1) ), : );
        else
          if singleCellFile == cellFileMap( idToCF(c,1) )
            color = cm( 1, : );
          else
            continue;
          end
        end
        
        cellFileMat = [ cellFileMat ; p1(1) p1(2) p1(3) ];
        
        % get neighbor for specific vertex ID = c
        nVec = getNeighbors( c, tri, numCells );
        
        % number of links of current cell
        numLinks = size( nVec, 2 );
        
        % initialize texture with zeros for each cell
        M = zeros(3);
        
        % loop over all linked neighbors
        for n=1:numLinks
          % vertex ID
          verID = nVec(1,n);
          
          % get the corresponding neighbor position
          p2 = [ tri.Points( verID, 1 ) tri.Points( verID, 2 ) tri.Points( verID, 3 ) ];
          
          % compute link matrix
          lMat = getLinkMatrix( p2-p1, p2-p1 );
          
          % update texture matrix
          M = M + lMat;
        end
        
        % after processing each neighbor, divide each entry by number
        % of neighbors -> averaging
        M = M./numLinks;
        
        % compute the eigenvectors and eigenvalues of matrix M
        % The columns of Q are the eigenvectors and the diagonal
        % elements of D are the eigenvalues
        [Q,D] = eig(M);
        
        % check if the eigenvalues are smaller then zero; if so, then do
        % not draw a line and consider the absolute value of it -> TODO
        positiveEigenvalue = [ 1 ; 1 ; 1 ];
        radii = diag(D);
        for e=1:3
          if radii(e) < 0.
            radii(e) = -radii(e);
            positiveEigenvalue( e, 1 ) = 0;
          end
        end
        
        % radii of the ellipsoid
        radii = sqrt( radii );
        
        xEigVec = Q(:, 1);
        yEigVec = Q(:, 2);
        zEigVec = Q(:, 3);
        
        % draw the single cell as ellipsoid
        [ x, y, z ] = ellipsoid( 0, 0, 0, radii(1)/2., radii(2)/2., radii(3)/2., nEllip );
        X = p1(1) + x*xEigVec(1) + y*yEigVec(1) + z*zEigVec(1);
        Y = p1(2) + x*xEigVec(2) + y*yEigVec(2) + z*zEigVec(2);
        Z = p1(3) + x*xEigVec(3) + y*yEigVec(3) + z*zEigVec(3);
        
        semiLines = [];
        % draw the three major axes in the origin which are then
        % rotated according to the eigenvectors
        for l=1:3
          if l == 1
            sX = [ -radii(1)/2., radii(1)/2. ];
            sY = [ 0, 0 ];
            sZ = [ 0, 0 ];
          elseif l == 2
            sX = [ 0, 0 ];
            sY = [ -radii(2)/2., radii(2)/2. ];
            sZ = [ 0, 0 ];
          else
            sX = [ 0, 0 ];
            sY = [ 0, 0 ];
            sZ = [ -radii(3)/2., radii(3)/2. ];
          end
          lineX = p1(1) + sX*xEigVec(1) + sY*yEigVec(1) + sZ*zEigVec(1);
          lineY = p1(2) + sX*xEigVec(2) + sY*yEigVec(2) + sZ*zEigVec(2);
          lineZ = p1(3) + sX*xEigVec(3) + sY*yEigVec(3) + sZ*zEigVec(3);
          
          projLine1 = projectOnPlane( [ lineX(1) lineY(1) lineZ(1) ], planePos, u, v );
          projLine2 = projectOnPlane( [ lineX(2) lineY(2) lineZ(2) ], planePos, u, v );
          projLine1 = transformPoint3d( projLine1, TF );
          projLine2 = transformPoint3d( projLine2, TF );
          
          % and store the start/end points of the lines in linePos
          semiLines = [ semiLines projLine1 projLine2 ];
        end
        
        % update semi axes in 3D
        linePos = [ linePos ; semiLines ];
        
        % project each vertex of the ellipsoid onto the plane
        dimP = size( X, 1 );
        for q=1:dimP
          for p=1:dimP
            curPos = projectOnPlane( [ X(p,q) Y(p,q) Z(p,q) ], planePos, u, v );
            curPos = transformPoint3d( [curPos(1) curPos(2) curPos(3)], TF );
            X(p,q) = curPos(1);
            Y(p,q) = curPos(2);
            Z(p,q) = curPos(3);
          end
        end
        
        p1 = projectOnPlane( p1, planePos, u, v );
        p1 = transformPoint3d( p1, TF );
        centerEllipse = [ centerEllipse ; p1 ];
        minMaxS = determineAxes( X, Y, Z, p1, dir );
        minMaxSemiAxisVector = [ minMaxSemiAxisVector ; minMaxS ];  
      end
      
      % draw principal components
      if renderPrincipalComponents == 1
        start = mean( cellFileMat );
        %coeff = pca( cellFileMat );
        arrowLength = 150;
        for a=1:3
          P(a) = quiver3( start(1), start(2), start(3),...
            coeff(1,a), coeff(2,a), coeff(3,a),...
            arrowLength, 'LineWidth', 3,...
            'Color', cm( a, : ) );
          hold on;
        end
      end
      
      if overlapping == 1
        zOffset = 1.;
      else
        zOffset = 0.;
      end
      
      dimL = size( centerEllipse, 1 );
      % draw ellipses and lines
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
        minorSemiAxes = cross( dir, maxSemiPoint-c );
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
        [ ELLIP(l) ELLIPPATCH(l) ] = drawEllipse3d( c(1), c(2), c(3)+l*zOffset+0.2, maxLength, minLength, 0, theta );
        set( ELLIP(l), 'color', [ 0 0 0 ], 'LineWidth', lineWidth );
        set( ELLIPPATCH(l), 'FaceColor', [ 1 1 1 ], 'FaceLighting', 'gouraud' );
        
        if strcmp( lineStr( 1, lineRenderType ), 'renderLargest2DElongation' )
          MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', 'r', 'LineWidth', lineWidth );
          hold on;
        elseif strcmp( lineStr( 1, lineRenderType ), 'renderAll2DElongation' )
          MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', 'r', 'LineWidth', lineWidth );
          hold on;
          MIN(l) = line( lineMinX, lineMinY, lineMinZ+l*zOffset+0.2, 'Color', 'r', 'LineWidth', lineWidth );
          hold on;
        elseif strcmp( lineStr( 1, lineRenderType ), 'renderLargest3DElongation' )
          lineMaxX = [ linePos( l, 13 ), linePos( l, 16 ) ];
          lineMaxY = [ linePos( l, 14 ), linePos( l, 17 ) ];
          lineMaxZ = [ linePos( l, 15 ), linePos( l, 18 ) ];
          MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', 'r', 'LineWidth', lineWidth );
          hold on;
        elseif strcmp( lineStr( 1, lineRenderType ), 'renderAll3DElongation' )
          lineMinX = [ linePos( l, 1 ), linePos( l, 4 ) ];
          lineMinY = [ linePos( l, 2 ), linePos( l, 5 ) ];
          lineMinZ = [ linePos( l, 3 ), linePos( l, 6 ) ];
          lineMidX = [ linePos( l, 7 ), linePos( l, 10 ) ];
          lineMidY = [ linePos( l, 8 ), linePos( l, 11 ) ];
          lineMidZ = [ linePos( l, 9 ), linePos( l, 12 ) ];
          lineMaxX = [ linePos( l, 13 ), linePos( l, 16 ) ];
          lineMaxY = [ linePos( l, 14 ), linePos( l, 17 ) ];
          lineMaxZ = [ linePos( l, 15 ), linePos( l, 18 ) ];
          MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', 'r', 'LineWidth', lineWidth );
          hold on;
          MID(l) = line( lineMidX, lineMidY, lineMidZ+l*zOffset+0.2, 'Color', 'k', 'LineWidth', lineWidth );
          hold on;
          MIN(l) = line( lineMinX, lineMinY, lineMinZ+l*zOffset+0.2, 'Color', 'b', 'LineWidth', lineWidth );
          hold on;
        end
      end
      
      hold off;
      set( f,'nextplot','replacechildren' );
      viewOffset = viewOffsets( dataIndex );
      axis( [ minAxes(1)-viewOffset maxAxes(1)+viewOffset...
        minAxes(2)-viewOffset maxAxes(2)+viewOffset...
        minAxes(3)-viewOffset maxAxes(3)+viewOffset...
        minAxes(3)-viewOffset maxAxes(3)+viewOffset ] );
      axis off
      daspect( [ 1 1 1 ] );
      
      grid off;
      xlabel('X');
      ylabel('Y');
      zlabel('Z');
      C = strsplit( char( dataStr( 1, dataIndex ) ), '_' );
      title( strcat( C( 1, 1 ), ' Time Step ', num2str(curT) ) );
      
      if strcmp( exportTypeStr( 1, exportType ), 'AsVideo' )
        writeVideo( writerObj, getframe(gcf) );
      end
      
      % if images instead of a video should be exported
      if strcmp( exportTypeStr( 1, exportType ), 'AsImages' )
        if imgStart < 10
          digit = strcat( dataStr( 1, dataIndex ), viewStr( 1, cView ), '_00' );
        elseif imgStart < 100
          digit = strcat( dataStr( 1, dataIndex ), viewStr( 1, cView ), '_0' );
        else
          digit = strcat( dataStr( 1, dataIndex ), viewStr( 1, cView ), '_' );
        end
        
        % output with number of cells
        filePath = strcat( imageDir, digit, num2str(imgStart), '_', num2str(numCells), '.png' );
        % output without number of cells
        %filePath = strcat( imageDir, digit, num2str(imgStart), '.png' );
        
        saveas( gcf, char(filePath) );
        
        imgStart = imgStart + 1;
      end
      
      % at last delete lightsource
      delete(findall(gcf, 'Type', 'light'))
      camlight headlight;
      
    end
  end
  
  if strcmp( exportTypeStr( 1, exportType ), 'AsVideo' )
    close(writerObj);
  end
  
  % close figure after processing
  %close(f);
  
end
