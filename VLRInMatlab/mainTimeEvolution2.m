%%%% The mesh deformation is based on the paper of Graner et al,
%%%% Discrete rearranging disordered patterns, part I: Robust
%%%% statistical tools in two or three dimensions, Eur. Phys. J.,
%%%% pp 349 -- 369, 2008

%%%%% setting of properties %%%%%%
% data Index:
% 1 -> 120830
% 2 -> 121204
% 3 -> 121211
% 4 -> 130508
% 5 -> 130607
% 6 -> 131203
% 7 -> all
data = 7;
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
renderTermType = 2;
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
% equalAspectRation
% This is required to use equal axes in order to project the ellipsoids
% correctly, else the result looks blurred (but only the vis, the computation
% is still correct)
equalAspectRatio = 0;
% render elongation lines
renderLines = 1;
% choose which major lines of the ellipsoids should be rendered
lineRenderType = 4;
% vector of line render types
% 1. draw only those lines with the largest elongation of the ellipoids in 3D
% 2. draw only those lines with the largest elongation of the ellipoids
%    projected onto the generated 2D plane
% 3: render all three major lines of elongation of the ellipsoids in 3D
lineStr = { 'renderOnlyLargest3DElongation'...
            'renderOnlyLargestElongation'...
            'all2DLines'...
            'allLines' };
% in fact the lambda factor how far the plane is translated
% along the viewing direction vector
planePositionFactor = 350;

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
if data ~= 7
  startD = data;
  endD = data;
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
  imageDir = strcat( 'images/', dataStr( 1, dataIndex ), 'Deformation' , '/' );
  
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
  ZCol = 2 * data{4};
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
    maxT = startT + 1;
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
  nEllip = 20;
  
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
  
  
  % boundary offset for rendering since the elongation of the
  % ellipsoids can be quite large; special case for data set 130508
  offset = 50;
  if strcmp( dataStr( 1, dataIndex ), '130508_raw' ) ||...
     strcmp( dataStr( 1, dataIndex ), '131203_raw' )
    offset = 150;
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
  
  % axes properties
  if equalAspectRatio == 0
    axis( [ minX-offset maxX+offset minY-offset maxY+offset minZ-offset maxZ+offset minZ-offset maxZ+offset ] );
    axis off
    daspect( [ 1 1 1 ] );
  else
    axis equal
    axis off
  end
  
  grid off;
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  C = strsplit( char( dataStr( 1, dataIndex ) ), '_' );
  title( strcat( C( 1, 1 ), ' Time Step ', num2str(startT) ) );
  camlight headlight;
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
  % line instance
  L = [];
  % PC instance
  P = [];
  
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
    set( S, 'Visible', 'off' );
    set( L, 'Visible', 'off' );
    set( P, 'Visible', 'off' );
    
    % if at least three cells exists in both time steps
    if numCellsC > 3 && numCellsN > 3
      % delaunay triangulation
      if begin ~= 1
        triC = triN;
        triN = delaunayTriangulation( matPosN(:,1), matPosN(:,2), matPosN(:,3) );
        numTotalLinksC = numTotalLinksN;
        numTotalLinksN = size( edges(triN), 1 );
        numTotalAveragedLinks = ( numTotalLinksC + numTotalLinksN )/2.;
      else
        triC = delaunayTriangulation( matPosC(:,1), matPosC(:,2), matPosC(:,3) );
        triN = delaunayTriangulation( matPosN(:,1), matPosN(:,2), matPosN(:,3) );
        numTotalLinksC = size( edges(triC), 1 );
        numTotalLinksN = size( edges(triN), 1 );
        numTotalAveragedLinks = ( numTotalLinksC + numTotalLinksN )/2.;
      end
      
      cellsInFileMat = [];
      linePos = [];
      lineColorIndex = [];
      
      % first get mapping of vertex ids of delaunay triangulation to object
      % ids of raw data set and store the results as a dimx2 matrix
      % including the link information between two cells
      objectLinksC = generateLinksBetweenObjectsId( edges( triC ), cellIdsC );
      objectLinksN = generateLinksBetweenObjectsId( edges( triN ), cellIdsN );
      
      % determine the number of conserved links for all cells
      conservedLinks = intersect( objectLinksC, objectLinksN, 'rows' );
      
      % disappeared links
      disappearedLinks = setxor( objectLinksC, conservedLinks, 'rows' );
      
      % appeared links
      appearedLinks = setxor( objectLinksN, conservedLinks, 'rows' );
      
      % consider each cell between two time steps
      % and render the time evolution as an ellipsoid located at the
      % position of the cell at time step t + deltaT
      % Note that we consider the TIME STEP t + deltaT and look back at
      % time step t which cell is related to the second time step
      for c=1:numCellsN
        % geometrical term
        B = zeros(3);
        
        % topological term
        T1 = zeros(3);
        T2 = zeros(3);
        
        % get position of current cell
        p1 = [ triN.Points( c, 1 ) triN.Points( c, 2 ) triN.Points( c, 3 ) ];
        
        % get objectId of current cell
        objectIdN = cellIdsN( c );
        
        % get objectId of previous cell
        % which could be in the simplest case the same
        % as the objectId of the current cell
        objectIdC = objectIdN;
        
        % get color for ellipsoid
        if renderSingleCellFile == 0
          color = cm( cellFileMap( objectIdN ), : );
        else
          if singleCellFile == cellFileMap( objectIdN )
            color = cm( 1, : );
          else
            continue;
          end
        end

        % check if the current cell already existed in the last time
        % step; if not then back traverse its precursors until the
        % corresponding object id is found
        if isConserved( objectIdN, objectLinksC ) == 0
          objectIdC = getPrecursorID( objectIdN, cellPrecursorsN( c ), cellIdsC );
        else
          objectIdC = objectIdN;
        end
        
        % get the position of the previous cell in the previous time step
        p2 = getCellPosition( objectIdC, triC, cellIdsC );
        
        cellsInFileMat = [ cellsInFileMat ; p2(1) p2(2) p2(3) ];
        
        % get all neighbors for the two found cells
        nVecC = getNeighborsOfObjectId( objectIdC, objectLinksC );
        nVecN = getNeighborsOfObjectId( objectIdN, objectLinksN );
        
        % vector of objectIds of cells which are conserved in the next time
        % step
        conservedLinksPerCell = getNeighborsOfObjectId( objectIdN, conservedLinks );
        % number of conserved links between two time steps
        numConservedLinksPerCell = size( conservedLinksPerCell, 2 );
        
        % get vector of objectsIds of cell which are added in the next time
        % step
        appearedLinksPerCell = getNeighborsOfObjectId( objectIdN, appearedLinks );
        % number of added links between two time steps
        numAppearedLinksPerCell = size( appearedLinksPerCell, 2 );
        
        % get vector of objectsIds of cell which disappeared in the next time
        % step
        disappearedLinksPerCell = getNeighborsOfObjectId( objectIdN, disappearedLinks );
        % number of added links between two time steps
        numDisappearedLinksPerCell = size( disappearedLinksPerCell, 2 );
        
        % number of links of current and next cell
        numTotalLinksPerCellC = size( nVecC, 2 );
        numTotalLinksPerCellN = size( nVecN, 2 );
        % number of averaged links
        numAveragedLinksPerCell = ( numTotalLinksPerCellC + numTotalLinksPerCellN )/2.;
        
        % loop over all linked and conserved neighbors
        % determining B
        for n=1:numConservedLinksPerCell
          % current object id of neighbor
          neighborId = conservedLinksPerCell( 1, n );
          
          % link at time step t
          l1 = getCellPosition( neighborId, triC, cellIdsC ) - p2;
          % link at time step t + deltaT
          l2 = getCellPosition( neighborId, triN, cellIdsN ) - p1;
          
          % compute average of the two links
          la = ( l1 + l2 )/2.;
          
          % compute difference between the two links
          ld = l2 - l1;
          
          % compute link matrix which is the matrix C in the paper (Appendix C1)
          C = getLinkMatrix( la, ld/deltaT );
          
          % sum of C and the transposed one
          B = B + C + C';
        end
        
        % after processing each neighbor, divide each entry by number
        % of conserved neighbors -> averaging
        if numConservedLinksPerCell > 0
          B = B./numConservedLinksPerCell;
        end
         
        % multiply the factor of N_c/N_tot (see paper in Appendix C1)
        if numConservedLinksPerCell > 0
          B = B.*( numConservedLinksPerCell / numAveragedLinksPerCell);
        end
        
        % loop over all appeared linked neighbors determining the first part of T
        for a=1:numAppearedLinksPerCell
          % current object id of neighbor
          neighborId = appearedLinksPerCell( 1, a );
          
          % link at time step t + deltaT
          l = getCellPosition( neighborId, triN, cellIdsN ) -...
              getCellPosition( objectIdN, triN, cellIdsN );
            
          % compute link matrix
          m = getLinkMatrix( l, l );
          
          % update matrix T
          T1 = T1 + m;
        end
        
        % NEW: also add the link between the cell in time step t + deltaT
        % and the precursor which does not have the same id
        if isConserved( objectIdN, objectLinksC ) == 0
          % link between position at time step t and t + deltaT
          l = p2 - p1;
            
          % compute link matrix
          m = getLinkMatrix( l, l );
          
          T1 = T1 + m;
          
          numAppearedLinksPerCell = numAppearedLinksPerCell + 1;
        end
        
        % after processing each neighbor, divide each entry by number
        % of appeared neighbors -> averaging
        if numAppearedLinksPerCell > 0
          T1 = T1./numAppearedLinksPerCell;
        end
        
        % multiply the factor of deltaN_a/N_tot (see paper in Appendix C1)
        if numAppearedLinksPerCell > 0
          T1 = T1.*( numAppearedLinksPerCell / numAveragedLinksPerCell);
        end
        
        % divide by deltaT
        T1 = T1./deltaT;
        
        % loop over all disappeared linked neighbors determining the second part of T
        for d=1:numDisappearedLinksPerCell
          % current object id of neighbor
          neighborId = disappearedLinksPerCell( 1, d );
          
          % link at time step t
          l = getCellPosition( neighborId, triC, cellIdsC ) -...
              getCellPosition( objectIdC, triC, cellIdsC );
            
          % compute link matrix
          m = getLinkMatrix( l, l );
          
          % update matrix T
          T2 = T2 + m;
        end
        
        % after processing each neighbor, divide each entry by number
        % of disappeared neighbors -> averaging
        if numDisappearedLinksPerCell > 0
          T2 = T2./numDisappearedLinksPerCell;
        end
        
        % multiply the factor of deltaN_a/N_tot (see paper in Appendix C1)
        if numDisappearedLinksPerCell > 0
          T2 = T2.*( numDisappearedLinksPerCell / numAveragedLinksPerCell);
        end
        
        % divide by deltaT
        T2 = T2./deltaT;
        
        % compute the final texture for the current ellipsoid
        if strcmp( termTypeStr( 1, renderTermType ), 'B' )
          M = B;
        elseif strcmp( termTypeStr( 1, renderTermType ), 'T' )
          M = T1 - T2;
        else
          M = B + T1 - T2;
        end
        
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
        
        % line color by default black for each line
        lineColor = [ 0 0 0 ; 0 0 0 ; 0 0 0 ];
        % set line colors depending on the computed term and
        % eigenvalue
        if strcmp( termTypeStr( 1, renderTermType ), 'B' )
          for e=1:3
            % if the eigenvalue is positive then use a red color
            if positiveEigenvalue( e, 1 ) == 1
              lineColor( e, 1 ) = 1;
              lineColor( e, 2 ) = 0;
              lineColor( e, 3 ) = 0;
            % negative eigenvalue
            elseif positiveEigenvalue( e, 1 ) == 0
              lineColor( e, 1 ) = 0;
              lineColor( e, 2 ) = 0;
              lineColor( e, 3 ) = 1;
            end
          end
        end
        
        % radii of the ellipsoid
        radii = sqrt( radii );
        
        % scaling such that even with small deformation changes
        % the ellipsoids can be identified
        radii = radii.*5;
        
        xEigVec = Q(:, 1);
        yEigVec = Q(:, 2);
        zEigVec = Q(:, 3);
        
        % draw the single cell as ellipsoid
        [ x, y, z ] = ellipsoid( 0, 0, 0, radii(1)/2., radii(2)/2., radii(3)/2., nEllip );
        %ellipPos = (p1+p2)/2.;
        ellipPos = p1;
        X = ellipPos(1) + x*xEigVec(1) + y*yEigVec(1) + z*zEigVec(1);
        Y = ellipPos(2) + x*xEigVec(2) + y*yEigVec(2) + z*zEigVec(2);
        Z = ellipPos(3) + x*xEigVec(3) + y*yEigVec(3) + z*zEigVec(3);
        
        if renderLines == 1
          % if only the lines with the largest elongation should be drawn
          % then only traverse the last entry of the loop which corresponds
          % to the largest eigenvalue
          lineStart = 1;
          if strcmp( lineStr( 1, lineRenderType ), 'renderOnlyLargest3DElongation' )
            lineStart = 3;
          end
          
          % draw the three major axes in the origin which are then
          % rotated according to the eigenvectors
          for l=lineStart:3
            
            if strcmp( termTypeStr( 1, renderTermType ), 'T' )
              % skip the lines that correspond to a negative eigenvalue
              if positiveEigenvalue( e, 1 ) == 0
                continue;
              end
            end
            
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
            lineX = ellipPos(1) + sX*xEigVec(1) + sY*yEigVec(1) + sZ*zEigVec(1);
            lineY = ellipPos(2) + sX*xEigVec(2) + sY*yEigVec(2) + sZ*zEigVec(2);
            lineZ = ellipPos(3) + sX*xEigVec(3) + sY*yEigVec(3) + sZ*zEigVec(3);
            
            % and store the start/end points of the lines in linePos
            linePos = [ linePos ; lineX(1) lineY(1) lineZ(1) ; lineX(2) lineY(2) lineZ(2) ];
            lineColorIndex = [ lineColorIndex ; lineColor( l, : ) ];
            
            % draw line in 3D
            %line( lineX, lineY, lineZ );
          end
        end
        
        if renderSingleCellFile == 1
          S(c) = surface( X, Y, Z, 'FaceColor', 'w',...
            'EdgeColor', 'none', 'EdgeAlpha', 0,...
            'FaceLighting', 'gouraud' );
        else
          S(c) = surface( X, Y, Z, 'FaceColor', color,...
            'EdgeColor', 'none', 'EdgeAlpha', 0,...
            'FaceLighting', 'gouraud' );
        end
        hold on;
      end
      
      % after drawing all cells perform a pca for getting the principal axes
      % if at least four points exist (such that three PCs are generated)
      numPCACells = size( cellsInFileMat, 1 );
      if numPCACells > 3
        start = mean( cellsInFileMat );
        %coeff = pca( cellsInFileMat );
        coeff = getPrincipalComponents( dataStr( 1, dataIndex ) );
        arrowLength = 150;
        
        % draw principal components
        if renderPrincipalComponents == 1
          for a=1:3
            P(a) = quiver3( start(1), start(2), start(3),...
              coeff(1,a), coeff(2,a), coeff(3,a),...
              arrowLength, 'LineWidth', 3,...
              'Color', cm( a, : ) );
            hold on;
          end
        end
        
        % set PC depending on the viewing direction
        if cView == 1
          dir = coeff(:,2);
          u = coeff(:,1);
          v = coeff(:,3);
        elseif cView == 2 || cView == 4
          dir = coeff(:,3);
          u = coeff(:,1);
          v = coeff(:,2);
        elseif cView == 3
          dir = coeff(:,1);
          u = coeff(:,2);
          v = coeff(:,3);
        end
        
        if renderLines == 1
          % set plane position
          planePos = dir * planePositionFactor;
          
          % draw lines of the major ellipsoid axes onto the projected plane
          l = 1;
          dimL = size( linePos, 1 ) + 1;
          while l < dimL
            if strcmp( lineStr( 1, lineRenderType ), 'renderOnlyLargestElongation' )
              % check the length of all three main major axes and choose the
              % largest one to render only
              lengthVec = zeros( 3, 1 );
              lineVectors = [];
              for ll=1:3
                linep1 = projectOnPlane( linePos(l, :), planePos, u, v );
                linep2 = projectOnPlane( linePos(l+1, :), planePos, u, v );
                lineVectors = [ lineVectors ; linep2(1)-linep1(1) linep2(1)-linep1(1) linep2(1)-linep1(1) ];
                length = norm( linep2 - linep1 );
                lengthVec( ll, 1 ) = length;
                l = l+2;
              end
              [ C, I ] = max( lengthVec );
              linep1 = projectOnPlane( linePos( l-2*(4-I), :), planePos, u, v );
              linep2 = projectOnPlane( linePos( l+1-2*(4-I), :), planePos, u, v );
              lineX = [ linep1(1), linep2(1) ];
              lineY = [ linep1(2), linep2(2) ];
              lineZ = [ linep1(3), linep2(3) ];
              index = int16( round( (l+1-2*(4-I))/2 ) );
              L(l) = line( lineX, lineY, lineZ, 'Color', lineColorIndex( index, : ), 'LineWidth', 1.5 );
              hold on;
            elseif strcmp( lineStr( 1, lineRenderType ), 'all2DLines' )
              % only render the line with the longest elongation and the one perpendicular
              % to it
              lengthVec = zeros( 3, 1 );
              lineVectors = [];
              for ll=1:3
                linep1 = projectOnPlane( linePos(l, :), planePos, u, v );
                linep2 = projectOnPlane( linePos(l+1, :), planePos, u, v );
                lineVectors = [ lineVectors ; linep2(1)-linep1(1) linep2(2)-linep1(2) linep2(3)-linep1(3) ];
                length = norm( linep2 - linep1 );
                lengthVec( ll, 1 ) = length;
                l = l+2;
              end
              
              % compute angle between
              a = dot( lineVectors( 1, : ), lineVectors( 2, : ) )
              a = dot( lineVectors( 1, : ), lineVectors( 3, : ) )
              a = dot( lineVectors( 2, : ), lineVectors( 3, : ) )
              
              [ C, I ] = max( lengthVec );
              linep1 = projectOnPlane( linePos( l-2*(4-I), :), planePos, u, v );
              linep2 = projectOnPlane( linePos( l+1-2*(4-I), :), planePos, u, v );
              lineX = [ linep1(1), linep2(1) ];
              lineY = [ linep1(2), linep2(2) ];
              lineZ = [ linep1(3), linep2(3) ];
              index = int16( round( (l+1-2*(4-I))/2 ) );
              L(l) = line( lineX, lineY, lineZ, 'Color', lineColorIndex( index, : ), 'LineWidth', 1.5 );
              hold on;
            else
              linep1 = projectOnPlane( linePos(l, :), planePos, u, v );
              linep2 = projectOnPlane( linePos(l+1, :), planePos, u, v );
              lineX = [ linep1(1), linep2(1) ];
              lineY = [ linep1(2), linep2(2) ];
              lineZ = [ linep1(3), linep2(3) ];
              index = int16( round( (l+1)/2 ) );
              L(l) = line( lineX, lineY, lineZ, 'Color', lineColorIndex( index, : ), 'LineWidth', 1.5 );
              hold on;
              l = l+2;
            end
          end
        end
        
        hold off;
        set( f,'nextplot','replacechildren' );
        if equalAspectRatio == 0
          axis( [ minX-offset maxX+offset minY-offset maxY+offset minZ-offset maxZ+offset minZ-offset maxZ+offset ] );
          %axis off
          daspect( [ 1 1 1 ] );
        else
          %axis( [ minX-offset maxX+offset minY-offset maxY+offset minZ-offset maxZ+offset minZ-offset maxZ+offset ] );
          axis equal
          axis off
        end
        
        %daspect('manual')
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
        % perhaps include optical rotation information into current view
        %rot = [ 2 -5 41 ];
        if cView == 1
          % top view
          if dataIndex == 2
            view( [ -coeff(1,2), -coeff(2,2), -coeff(3,2) ] );
          else
            view( [ coeff(1,2), coeff(2,2), coeff(3,2) ] );
          end
          camup( [ coeff(:,3) ] );
        elseif cView == 2
          % side view
          view( [ coeff(1,3), coeff(2,3), coeff(3,3) ] );
          if dataIndex == 2
            camup( [ -coeff(:,2) ] );
          else
            camup( [ coeff(:,2) ] );
          end
        elseif cView == 3
          % radial view
          view( [ coeff(1,1), coeff(2,1), coeff(3,1) ] );
          if dataIndex == 2
            camup( [ -coeff(:,2) ] );
          else
            camup( [ coeff(:,2) ] );
          end
        else
          % 3D view
          view(3);
        end
      else
        % if too few cells are available for PCA then just continue
        curT = curT + deltaT;
        continue;
      end
      
      % first delete last lightsource
      delete(findall(gcf, 'Type', 'light'))
      camlight headlight;
      
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
