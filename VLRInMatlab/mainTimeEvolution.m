%%%% The mesh deformation is based on the paper of Graner et al,
%%%% Discrete rearranging disordered patterns, part I: Robust
%%%% statistical tools in two or three dimensions, Eur. Phys. J.,
%%%% pp 349 -- 369, 2008

geomPath = strcat( pwd, '/geom3d' );
addpath( geomPath );

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
% 4 -> 3D
cView = 2;
% start with the current time step
startT = 1;
% deltaT value based on the paper mentioned above
% have to be a divider of the max time step value!
deltaT = 10;
% decide which term should be included in the time evolution
renderTermType = 3;
termTypeStr = { 'B' 'T' 'All' };
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
% choose which major lines of the ellipsoids should be rendered
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
% 1 -> convex hull
% 2 -> alpha shape
triangulationType = 2;
% vector of data strings
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' '131203_raw' };

% vector of view strings
viewStr = { 'Top' 'Side' 'Radial' '3D' };

% master cell file information taken by the picture made of Daniel located in dropbox
% and the trackGroup information of the raw data sets
%masterCellFile = [ 4 3 4 2 3 ];
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
  if triangulationType == 1
    imageDir = strcat( 'images/', dataStr( 1, dataIndex ),...
      'Deformation_CH', '/' );
  else
    imageDir = strcat( 'images/', dataStr( 1, dataIndex ),...
      'Deformation_AS' , '/' );
  end
  
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
  % string of precursors
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
  
  % get maximum of time steps
  % if a movie should be created
  if strcmp( exportTypeStr( 1, exportType ), 'AsVideo' ) ||...
      strcmp( exportTypeStr( 1, exportType ), 'AsImages' )
    maxT = max( TCol );
    % else only render between two time steps
  else
    maxT = startT + deltaT;
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
  % ObjectID X Y Z Timepoint LineageID TrackGroup Precursors
  %cellData = double.empty( 0, 8 );
  cellData = cell( dim, 8 );
  cellFileMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
  
  l = 1;
  nl = 1;
  % interpolate the missing positions in between
  while (l < numLines+1)
    firstCellId = IdCol(l);
    secondCellId = IdCol(l+1);
    % insert first line
    cellData( nl, : ) = {firstCellId XCol(l) YCol(l) ZCol(l) TCol(l) LCol(l) CFCol(l) PCol(l)};
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
        cellData( nl, : ) = {firstCellId x y z t LCol(l) CFCol(l) PCol(l)};
        nl = nl+1;
      end
      
      % insert last line
      cellData( nl, : ) = {firstCellId XCol(l+1) YCol(l+1) ZCol(l+1) TCol(l+1) LCol(l+1) CFCol(l+1) PCol(l+1)};
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
  %set( gcf, 'Renderer', 'zbuffer' );
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
  % projected surface instance
  PS = [];
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
  minAxes = applyTransformations( [ minX minY minZ ], planePos, u, v, TF, dataStr( 1, dataIndex ) );
  maxAxes = applyTransformations( [ maxX maxY maxZ ], planePos, u, v, TF, dataStr( 1, dataIndex ) );
  
  % after transformation the individual coords of min and max
  % may be switched
  for mm=1:3
    if minAxes(mm) > maxAxes(mm)
      tmp = minAxes(mm);
      minAxes(mm) = maxAxes(mm);
      maxAxes(mm) = tmp;
    end
  end
  
  % get the alpha shape radii for all time steps
  if triangulationType == 2
    alphaRadiiVector = getAlphaRadius( dataStr( 1, dataIndex ) );
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Traversal over all time steps and time evolution generation %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % this variable is set to zero after the first traversal of the loop
  % since it indicated the begin of the loop; this is required in order to
  % set the information of the next time step to the infos of the current
  % one in each time step which prevent redundant computations, but this is
  % valid for all loop traversals except the first one
  begin = 1;
  
  % the variables of two successive time steps are added a 'C' (Current)
  % or 'N' (Next) at the end
  % number of cells for the current and next time step
  numCellsC = 0;
  numCellsN = 0;
  
  % matrix of positions for current and next time step
  matPosC = [];
  matPosN = [];
  
  % vector of object ids in order to access the cell
  % file information in the cell file map
  cellIdsC = [];
  cellIdsN = [];
  
  % vector of strings storing the precursors
  cellPrecursorsN = [];
  
  % delaunay triangulation
  triC = [];
  triN = [];
  
  % unique edges in the triangulation
  uniqueEdgesC = [];
  uniqueEdgesN = [];
  
  % number of total links in current and next time step
  numTotalLinksC = 0;
  numTotalLinksN = 0;
  numTotalAveragedLinks = 0;
  
  % loop over all time steps
  imgStart = 1;
  curT = startT;
  while curT < maxT-deltaT+1
    % update cell information
    if begin ~= 1
      numCellsC = numCellsN;
      numCellsN = 0;
      matPosC = matPosN;
      matPosN = [];
      cellIdsC = cellIdsN;
      cellIdsN = [];
      cellPrecursorsN = [];
      
      for j=1:dim
        % next time step
        if cellData{ j, 5 } == curT+deltaT
          pos = [ cellData{ j, 2 } cellData{ j, 3 } cellData{ j, 4 } ];
          matPosN = [matPosN; pos];
          numCellsN = numCellsN +1;
          cellIdsN = [ cellIdsN ; cellData{ j, 1 } ];
          cellPrecursorsN = [ cellPrecursorsN ; cellData{ j, 8 } ];
        end
      end
    else
      for j=1:dim
        % current time step
        if cellData{ j, 5 } == curT
          pos = [ cellData{ j, 2 } cellData{ j, 3 } cellData{ j, 4 } ];
          matPosC = [matPosC; pos];
          numCellsC = numCellsC +1;
          cellIdsC = [ cellIdsC ; cellData{ j, 1 } ];
          % next time step
        elseif cellData{ j, 5 } == curT+deltaT%curT+deltaT-1
          pos = [ cellData{ j, 2 } cellData{ j, 3 } cellData{ j, 4 } ];
          matPosN = [matPosN; pos];
          numCellsN = numCellsN +1;
          cellIdsN = [ cellIdsN ; cellData{ j, 1 } ];
          cellPrecursorsN = [ cellPrecursorsN ; cellData{ j, 8 } ];
        end
      end
    end
    
    % clean figure content by removing the
    % last ellipsoids and lines in the previous time step
    if ishandle( S )
      set( S, 'Visible', 'off' );
    end
    if ishandle( PS )
      set( PS, 'Visible', 'off' );
    end
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
    if ishandle( P )
      set( P, 'Visible', 'off' );
    end
    
    % if at least three cells exists in both time steps
    if numCellsC > 3 && numCellsN > 3
      if begin ~= 1
        uniqueEdgesC = uniqueEdgesN;
        triC = triN;
        if triangulationType == 1
          % delaunay triangulation
          triN = delaunayTriangulation( matPosN(:,1), matPosN(:,2), matPosN(:,3) );
          uniqueEdgesN = edges(triN);
        else
          % alpha shape triangulation
          [VolN,ShapeN] = alphavol( [ matPosN(:,1) matPosN(:,2) matPosN(:,3) ], sqrt( alphaRadiiVector( curT, 1 )) );
          triN = ShapeN.tri;
          uniqueEdgesN = getUniqueEdges( triN );
        end
        numTotalLinksC = numTotalLinksN;
        numTotalLinksN = size( uniqueEdgesN, 1 );
        numTotalAveragedLinks = ( numTotalLinksC + numTotalLinksN )/2.;
      else
        if triangulationType == 1
          % delaunay triangulation
          triC = delaunayTriangulation( matPosC(:,1), matPosC(:,2), matPosC(:,3) );
          uniqueEdgesC = edges( triC );
          triN = delaunayTriangulation( matPosN(:,1), matPosN(:,2), matPosN(:,3) );
          uniqueEdgesN = edges( triN );
        else
          % alpha shape triangulation
          [VolC,ShapeC] = alphavol( [ matPosC(:,1) matPosC(:,2) matPosC(:,3) ], sqrt( alphaRadiiVector( curT, 1 )) );
          triC = ShapeC.tri;
          uniqueEdgesC = getUniqueEdges( triC );
          [VolN,ShapeN] = alphavol( [ matPosN(:,1) matPosN(:,2) matPosN(:,3) ], sqrt( alphaRadiiVector( curT, 1 )) );
          triN = ShapeN.tri;
          uniqueEdgesN = getUniqueEdges( triN );
        end
        numTotalLinksC = size( uniqueEdgesC, 1 );
        numTotalLinksN = size( uniqueEdgesN, 1 );
        numTotalAveragedLinks = ( numTotalLinksC + numTotalLinksN )/2.;
      end
      
      % compute time evolution for current deltaT and time step
      [cellsInFileMat, lineColorIndex, linePos, minMaxEigenValueIndex,...
        positiveEigenvalueVector, minMaxSemiAxisVector, centerEllipse, indexColorSet ]...
        = computeTimeEvolution( uniqueEdgesC, uniqueEdgesN, cellIdsC, cellIdsN,...
        numCellsN, triC, triN, matPosC, matPosN, cellPrecursorsN, triangulationType,...
        termTypeStr( 1, renderTermType ), dataStr( 1, dataIndex ), planePos, u, v, TF, deltaT,...
        renderSingleCellFile, singleCellFile, cellFileMap );
      
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
      
      % draw ellipses and lines
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
        [ ELLIP(l) ELLIPPATCH(l) ] = drawEllipse3d( c(1), c(2), c(3)+l*zOffset+0.2, maxLength, minLength, 0, theta );
        set( ELLIP(l), 'color', [ 0 0 0 ], 'LineWidth', lineWidth );
        set( ELLIPPATCH(l), 'FaceColor', [ 1 1 1 ], 'FaceLighting', 'none' );
        
        if strcmp( lineStr( 1, lineRenderType ), 'renderLargest2DElongation' )
          colorIndex = indexColorSet( l, 2 );
          if colorIndex == 0
            color = [ 0 0 1 ];
          else
            color = [ 1 0 0 ];
          end
          
          if strcmp( termTypeStr( 1, renderTermType ), 'T' ) ||...
              strcmp( termTypeStr( 1, renderTermType ), 'All' )
            color = [ 0 0 0 ];
          end
          MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
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
            MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
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
            MIN(l) = line( lineMinX, lineMinY, lineMinZ+l*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
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
          MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
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
              MIN(l) = line( lineX, lineY, lineZ+l*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
            elseif el == 2
              MID(l) = line( lineX, lineY, lineZ+l*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
            else
              MAX(l) = line( lineX, lineY, lineZ+l*zOffset+0.2, 'Color', color, 'LineWidth', lineWidth );
            end
            hold on;
          end
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
      if begin == 1
        %title( strcat( C( 1, 1 ), ' Time Step ', num2str(curT-1) ) );
        title( strcat( C( 1, 1 ), ' Time Step ', num2str(curT+deltaT) ) );
      else
        title( strcat( C( 1, 1 ), ' Time Step ', num2str(curT+deltaT) ) );
      end
    else
      % if too few cells are available for PCA then just continue
      curT = curT + deltaT;
      continue;
    end
    
    if strcmp( exportTypeStr( 1, exportType ), 'AsVideo' )
      writeVideo( writerObj, getframe(f) );
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
      
      filePath = strcat( imageDir, digit, num2str(imgStart), '.png' );
      
      saveas( gcf, char(filePath) );
      
      imgStart = imgStart + 1;
    end
    
    % at last set begin to false
    if begin == 1
      begin = 0;
      %curT = curT + deltaT - 1;
      curT = curT + deltaT;
    else
      curT = curT + deltaT;
    end
  end
  
  if strcmp( exportTypeStr( 1, exportType ), 'AsVideo' )
    close(writerObj);
  end
  
  % close figure after processing
  %close(f);
  
end
