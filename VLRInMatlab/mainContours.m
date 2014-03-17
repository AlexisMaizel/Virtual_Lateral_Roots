%%%% The mesh deformation is based on the paper of Graner et al,
%%%% Discrete rearranging disordered patterns, part I: Robust
%%%% statistical tools in two or three dimensions, Eur. Phys. J.,
%%%% pp 349 -- 369, 2008

addpath( '/home/necrolyte/Data/VLR/Virtual_Lateral_Roots/VLRInMatlab/geom3d' );

%%%%% setting of properties %%%%%%
% color range for different data sets
% { 'm' 'c' 'r' 'g' 'b' 'k' };
colors = [ [ 1 0 1 ]; [ 0 1 1 ]; [ 1 0 0 ]; [ 0 1 0 ]; [ 0 0 1 ]; [ 0 0 0 ] ];
% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
cView = 2;
% start with the current time step
startT = 300;
% startIndex
startI = 1;
% endIndex
endI = 20;
% draw delaunay tri?
drawDelaunay = 0;
% if set to one then only a single time step
% is rendered given by startT
exportType = 2;
% vector of data strings
exportTypeStr = { 'SingleFigure' 'AsImages' };
% render only master file?
renderSingleCellFile = 0;
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
visType = 3;
renderCellRanges = 0;
% only generate output images if the number of cells are within a desired
% range
cellRange = [ 20 40 60 80 100 120 140 ];
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
triangulationType = 1;
% vector of data strings
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' '131203_raw' };
% vector of view strings
viewStr = { 'Top' 'Side' 'Radial' '3D' };
% master cell file information taken by the picture made of Daniel located in dropbox
% and the trackGroup information of the raw data sets
%masterCellFile = [ 4 3 4 2 3 0 ];
masterCellFile = [ 4 4 4 3 3 0 ];
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
endData = 6;
% num of data
numData = 6;
% number of cell files
numCellFiles = 6;%double( maxCF-minCF+1 );

% cellData is the main array with all relevant information
% for the further analysis:
% ObjectID X Y Z Timepoint LineageID TrackGroup
cellDatas = cell( numData );
cellFileMap = cell( numData );
maxT = zeros( numData, 1 );
minAxes = zeros( numData, 3 );
maxAxes = zeros( numData, 3 );
centerData = zeros( numData, 3 );
totalMinAxes = [ 5000 5000 5000 ];
totalMaxAxes = [ -5000 -5000 -5000 ];
numCellsPerTimeStep = cell( numData, 1 );
dimData = zeros( numData, 1 );

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

% preprocessing of data sets to store and determine different properties
% loop over all data sets
%for dataIndex=1:numData
for dataIndex=startData:endData
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
    
%     % get maximum of time steps
%     % if a movie should be created
%     if strcmp( exportTypeStr( 1, exportType ), 'AsImages' )
%       maxT( dataIndex ) = max( TCol );
%       % else only render current time step
%     else
%       maxT( dataIndex ) = startT;
%     end
    
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
    
    l = 1;
    % first determine dimension of cellData
    % in order to initialize the size of cellData
    % -> huge performance speed up
    while (l < numLines+1)
      firstCellId = IdCol(l);
      secondCellId = IdCol(l+1);
      if firstCellId == secondCellId
        firstTS = TCol(l);
        secondTS = TCol(l+1);
        deltaTS = secondTS - firstTS;
        dimData( dataIndex ) = dimData( dataIndex ) + deltaTS + 1;
        l = l+2;
      else
        dimData( dataIndex ) = dimData( dataIndex )+1;
        l = l+1;
      end
    end
    
    cellData = cell( dimData( dataIndex ), 7 );
    cellFileMap{ dataIndex } = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
    
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
        cellFileMap{ dataIndex }( firstCellId ) = CFCol(l);
        
        % increment loop index
        l = l+2;
        % else cell exists only for one time step
      else
        % update cell file map
        cellFileMap{ dataIndex }( firstCellId ) = CFCol(l);
        
        l = l+1;
      end
    end
    
    cellDatas{ dataIndex } = cellData;
    
    % initialize number of rows depending on the number of time steps of
    % current data set
    maxT( dataIndex ) = max( TCol );
    numCellsPerTimeStep{dataIndex} = zeros(maxT(dataIndex), 1);
    
    % get maximum number of cells for each data set and time step
    for j=1:dimData( dataIndex )
      timeStep = cellDatas{dataIndex}{j, 5};
      numCellsPerTimeStep{dataIndex}(timeStep, 1) = numCellsPerTimeStep{dataIndex}(timeStep, 1) + 1;
    end
    
    if strcmp( visualizationType( 1, visType ), 'Ellipses' ) ||...
        strcmp( visualizationType( 1, visType ), 'Contour' )
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

      % set axes properties
      minAxes(dataIndex, :) = projectOnPlane( [ minX minY minZ ], planePos, u, v );
      minAxes(dataIndex, :) = transformPoint3d( minAxes(dataIndex, :), TF );
      maxAxes(dataIndex, :) = projectOnPlane( [ maxX maxY maxZ ], planePos, u, v );
      maxAxes(dataIndex, :) = transformPoint3d( maxAxes(dataIndex, :), TF );
      
      % after transformation the individual coords of min and max
      % may be switched
      for mm=1:3
        if minAxes(dataIndex, mm) > maxAxes(dataIndex, mm)
          tmp = minAxes(dataIndex, mm);
          minAxes(dataIndex, mm) = maxAxes(dataIndex, mm);
          maxAxes(dataIndex, mm) = tmp;
        end
      end      
    elseif strcmp( visualizationType( 1, visType ), 'Ellipsoids' )
      minAxes(dataIndex, :) = [ minX minY minZ ];
      maxAxes(dataIndex, :) = [ maxX maxY maxZ ];
    end
    
    for mm=1:3
      if minAxes(dataIndex, mm) < totalMinAxes(mm)
        totalMinAxes(mm) = minAxes(dataIndex, mm);
      end
      if maxAxes(dataIndex, mm) >= totalMaxAxes(mm)
        totalMaxAxes(mm) = maxAxes(dataIndex, mm);
      end
    end
end

% get center of total min and max values of data sets
centerAllData = [ (totalMaxAxes(1) + totalMinAxes(1))/2. ...
  (totalMaxAxes(2) + totalMinAxes(2))/2. ...
  (totalMaxAxes(3) + totalMinAxes(3))/2. ];

% and translate total min and max values according to center
totalMinAxes = totalMinAxes - centerAllData;
totalMaxAxes = totalMaxAxes - centerAllData;

if strcmp( visualizationType( 1, visType ), 'Contour' )
  % contour instance
  CONTOUR = [];
elseif strcmp( visualizationType( 1, visType ), 'Ellipsoids' )
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

% loop over all normalized steps
for curI=startI:endI
  
  if strcmp( visualizationType( 1, visType ), 'Contour' )
    if ishandle( CONTOUR )
      set( CONTOUR, 'Visible', 'off' );
    end
  elseif strcmp( visualizationType( 1, visType ), 'Ellipsoids' )
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
  
  if renderPrincipalComponents == 1
    if ishandle( P )
      set( P, 'Visible', 'off' );
    end
  end
  
  leg = zeros( numData, 1 );
  
  % loop over all data sets
  for dataIndex=startData:endData
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
    
    % number of cells for the current time step
    numCells = 0;
    
    % matrix of positions for current time step
    matPos = [];
    
    % vector of object ids in order to access the cell
    % file information in the cell file map
    idToCF = [];
    
    % NEW: only continue with this time step that is synchronized in
    % number of cells for each data set
    numNormCells = getNormalizedCellNumber( curI, 18, 143, 1, 20 );
    curT = 1;
    for j=1:maxT(dataIndex)
      if numNormCells - epsilon < numCellsPerTimeStep{dataIndex}(j,1)...
          && numNormCells + epsilon > numCellsPerTimeStep{dataIndex}(j,1)
        curT = j;
        break;
      end
    end
    
    for j=1:dimData( dataIndex )
      if cellDatas{ dataIndex }{ j, 5 } == curT
        pos = [ cellDatas{ dataIndex }{ j, 2 } cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
        matPos = [matPos; pos];
        numCells = numCells +1;
        idToCF = [ idToCF ; cellDatas{ dataIndex }{ j, 1 } ];
      end
    end
    
    % if at least three cells exists
    if numCells > 3
      
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
      
      if triangulationType == 1
        % delaunay triangulation
        tri = delaunayTriangulation( matPos(:,1), matPos(:,2), matPos(:,3) );
      else
        % alpha shape triangulation
        [Vol,Shape] = alphavol( [ matPos(:,1) matPos(:,2) matPos(:,3) ], sqrt( alphaRadiiVector( curT, 1 )) );
        tri = Shape.tri;
      end
      
      % number of total links in current time step
      %numTotalLinks = size( edges(tri), 1 );
      
      if drawDelaunay == 1 && triangulationType == 1
        tetramesh( tri, 'FaceColor', cm( 1, : ), 'FaceAlpha', 0.9 );
      end
      
      cellFileMat = [];
      linePos = [];
      minMaxSemiAxisVector = [];
      centerEllipse = [];
      
      % compute center of data
      centerData = [ 0 0 0 ];
      numConsideredCells = 0;
      for c=1:numCells
        % get position of current cell
        if triangulationType == 1
          p1 = [ tri.Points( c, 1 ) tri.Points( c, 2 ) tri.Points( c, 3 ) ];
        else
          p1 = [ matPos( c, 1 ) matPos( c, 2 ) matPos( c, 3 ) ];
        end
        
        % get color for ellipsoid
        if renderSingleCellFile == 1
          if singleCellFile ~= cellFileMap{ dataIndex }( idToCF(c,1) )
            continue;
          else
            numConsideredCells = numConsideredCells + 1;
          end
        else
          numConsideredCells = numConsideredCells + 1;
        end
        
         centerData = centerData + p1;
      end
      
      if numConsideredCells > 0
        centerData = centerData./numConsideredCells;
      end
      
      % draw an ellipsoid for each cell
      for c=1:numCells
        % get position of current cell
        if triangulationType == 1
          p1 = [ tri.Points( c, 1 ) tri.Points( c, 2 ) tri.Points( c, 3 ) ];
        else
          p1 = [ matPos( c, 1 ) matPos( c, 2 ) matPos( c, 3 ) ];
        end
        
        p1 = p1 - centerData;
        
        % get color for ellipsoid
        if renderSingleCellFile == 1
          if singleCellFile ~= cellFileMap{ dataIndex }( idToCF(c,1) )
            continue;
          end
        end
        
        cellFileMat = [ cellFileMat ; p1(1) p1(2) p1(3) ];
        
        % get neighbor for specific vertex ID = c
        if triangulationType == 1
          nVec = getNeighbors( c, tri, numCells );
        else
          nVec = getAlphaShapeNeighbors( c, tri );
        end
        
        % number of links of current cell
        numLinks = size( nVec, 2 );
        
        if size( nVec, 1 ) == 0
          continue;
        end
        
        % initialize texture with zeros for each cell
        M = zeros(3);
        
        % loop over all linked neighbors
        for n=1:numLinks
          % vertex ID
          verID = nVec(1,n);
          
          % get the corresponding neighbor position
          if triangulationType == 1
            p2 = [ tri.Points( verID, 1 ) tri.Points( verID, 2 ) tri.Points( verID, 3 ) ];
          else
            p2 = [ matPos( verID, 1 ) matPos( verID, 2 ) matPos( verID, 3 ) ];
          end
          
          p2 = p2 - centerData;
          
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
        
        if strcmp( visualizationType( 1, visType ), 'Ellipsoids' )
          S(c) = surface( X, Y, Z, 'FaceColor', colors( dataIndex, : ),...
            'EdgeColor', 'none', 'EdgeAlpha', 0,...
            'FaceLighting', 'gouraud' );
        end
        
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
        
%         PS(c) = surface( X, Y, Z, 'FaceColor', colors( dataIndex, : ),...
%           'EdgeColor', 'none', 'EdgeAlpha', 0,...
%           'FaceLighting', 'none' );
        
        p1 = projectOnPlane( p1, planePos, u, v );
        p1 = transformPoint3d( p1, TF );
        centerEllipse = [ centerEllipse ; p1 ];
        % the direction is now the normal of the x-y plane
        minMaxS = determineAxes( X, Y, Z, p1, [ 0 0 1 ] );
        minMaxSemiAxisVector = [ minMaxSemiAxisVector ; minMaxS ];
      end
      
      % draw principal components
      if renderPrincipalComponents == 1
        start = mean( cellFileMat );
        %coeff = pca( cellFileMat )
        if strcmp( visualizationType( 1, visType ), 'Ellipses' ) ||...
            strcmp( visualizationType( 1, visType ), 'Contour' )
          start = projectOnPlane( start, planePos, u, v );
          start = transformPoint3d( start, TF );
        end
        arrowLength = 150;
        for a=1:3
          if strcmp( visualizationType( 1, visType ), 'Ellipses' ) ||...
              strcmp( visualizationType( 1, visType ), 'Contour' )
            endP = projectOnPlane( coeff(:, a), planePos, u, v );
            endP = transformPoint3d( endP, TF );
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
      
      if strcmp( visualizationType( 1, visType ), 'Ellipses' )
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
          [ ELLIP(l) ELLIPPATCH(l) ] = drawEllipse3d( c(1), c(2), c(3)+l*zOffset+0.2, maxLength, minLength, 0, theta );
          set( ELLIP(l), 'color', [ 0 0 0 ], 'LineWidth', lineWidth );
          set( ELLIPPATCH(l), 'FaceColor', [ 1 1 1 ], 'FaceLighting', 'none' );
          
          if strcmp( lineStr( 1, lineRenderType ), 'renderLargest2DElongation' )
            hold on;
            MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', 'k', 'LineWidth', lineWidth );
          elseif strcmp( lineStr( 1, lineRenderType ), 'renderAll2DElongation' )
            hold on;
            MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', 'k', 'LineWidth', lineWidth );
            hold on;
            MIN(l) = line( lineMinX, lineMinY, lineMinZ+l*zOffset+0.2, 'Color', 'k', 'LineWidth', lineWidth );
          elseif strcmp( lineStr( 1, lineRenderType ), 'renderLargest3DElongation' )
            lineMaxX = [ linePos( l, 13 ), linePos( l, 16 ) ];
            lineMaxY = [ linePos( l, 14 ), linePos( l, 17 ) ];
            lineMaxZ = [ linePos( l, 15 ), linePos( l, 18 ) ];
            hold on;
            MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', 'r', 'LineWidth', lineWidth );
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
            hold on;
            MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', 'r', 'LineWidth', lineWidth );
            hold on;
            MID(l) = line( lineMidX, lineMidY, lineMidZ+l*zOffset+0.2, 'Color', 'k', 'LineWidth', lineWidth );
            hold on;
            MIN(l) = line( lineMinX, lineMinY, lineMinZ+l*zOffset+0.2, 'Color', 'b', 'LineWidth', lineWidth );
          end
        end
      elseif strcmp( visualizationType( 1, visType ), 'Contour' )
        if size( centerEllipse, 1 ) > 2
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
  end
  
  hold off;
  set( f,'nextplot','replacechildren' );
  viewOffset = 100;%viewOffsets( dataIndex );
%   axis( [ totalMinAxes(1)-viewOffset totalMaxAxes(1)+viewOffset...
%     totalMinAxes(2)-viewOffset totalMaxAxes(2)+viewOffset...
%     totalMinAxes(3)-viewOffset*3 totalMaxAxes(3)+viewOffset*3 ...
%     totalMinAxes(3)-viewOffset totalMaxAxes(3)+viewOffset ] );
    axis( [ totalMinAxes(1)-viewOffset totalMaxAxes(1)+viewOffset...
    -120 120 ...
    totalMinAxes(3)-viewOffset*3 totalMaxAxes(3)+viewOffset*3 ...
    totalMinAxes(3)-viewOffset totalMaxAxes(3)+viewOffset ] );
  %axis equal;
  axis off
  daspect( [ 1 1 1 ] );
  
  % legend
  hold on;
  leg = legend( '120830', '121204', '121211', '130508', '130607', '131203' );
  set(leg, 'Location', 'NorthWestOutside');
  linecolors = { [ 1 0 1 ], [ 0 1 1 ], [ 1 0 0 ], [ 0 1 0 ], [ 0 0 1 ], [ 0 0 0 ] };
  legendlinestyles( leg, {}, {}, linecolors );
  
  grid off;
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  title( strcat( 'Normalized Step ', num2str(curI), '-Cells', num2str(numNormCells) ) );
  
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
    filePath = strcat( imageDir, digit, num2str(curI), '-Cells', num2str(numNormCells), '.png' );
    
    saveas( gcf, char(filePath) );
  end
end
