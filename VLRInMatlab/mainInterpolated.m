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
% 6 -> all
data = 5;
% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
% 4 -> 3D
cView = 2;
% start with the current time step
startT = 1;
% deltaT value based on the paper mentioned above
deltaT = 1;
% draw ellipsoids?
drawEllipsoids = 1;
% draw delaunay tri?
drawDelaunay = 0;
% render movie?
% if set to zero then only a single time step
% is rendered given by startT
renderMovie = 1;
% create images for each time step
% TODO: does not work at the moment, let it be 0
createImages = 0;
% render only master file?
renderSingleCellFile = 1;
% render principal components
renderPrincipalComponents = 0;
% equalAspectRation
% This is required to use equal axes in order to project the ellipsoids
% correctly, else the result looks blurred (but only the vis, the computation
% is still correct)
equalAspectRatio = 0;
% choose which major lines of the ellipsoids should be rendered
lineRenderType = 2;
% vector of line render types
% 1. draw only those lines with the largest elongation of the ellipoids in 3D
% 2. draw only those lines with the largest elongation of the ellipoids
%    projected onto the generated 2D plane
% 3: render all three major lines of elongation of the ellipsoids in 3D
lineStr = { 'renderOnlyLargest3DElongation'...
            'renderOnlyLargestElongation'...
            'allLines' };
% in fact the lambda factor how far the plane is translated
% along the viewing direction vector
planePositionFactor = 350;

% vector of data strings
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };

% vector of view strings
viewStr = { 'Top' 'Side' 'Radial' '3D' };

% master cell file information taken by the picture made of Daniel located in dropbox
% and the trackGroup information of the raw data sets
%masterCellFile = [ 4 3 4 2 3 ];
masterCellFile = [ 4 4 4 3 3 ];

% either loop over all data sets creating five figures or only
% show one specific one
startD = 1;
endD = 5;
if data ~= 6
  startD = data;
  endD = data;
end

% loop over all data sets
for dataIndex=startD:endD
  
  % render the corresponding master cell file
  singleCellFile = masterCellFile( 1, dataIndex );
  
  % path to movie output
  movieDir = strcat( 'videos/', dataStr( 1, dataIndex ), '_',...
                     viewStr( 1, cView ),  '_movie.avi');
  
  % path to image output
  imageDir = strcat( 'videos/', dataStr( 1, dataIndex ),...
                     viewStr( 1, cView ), '/');
  
  % output format of values
  format longG
  
  % reading raw data
  path = strcat( '../FinalVLRForMatlab/', dataStr( 1, dataIndex ), '.csv' );
  fileID = fopen( char(path) );
  % format of data sets:
  % ObjectID X Y Z Timepoint Radius Precursors Color Lineage TrackID TrackColor TrackGroup Layer DivisionType
  formatSpec = '%d %f %f %f %d %d %q %q %d %q %q %d %d %q';
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
  % Lineage Tree Id
  LCol = data{9};
  % Cell File Id
  CFCol = data{12};
  % Cell Layer Id
  CLCol = data{13};
  % Divison Type
  %DCol = data{14};
  
  % get maximum of time steps
  % if a movie should be created
  if renderMovie == 1
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
  
  % boundary offset for rendering since the elongation of the
  % ellipsoids can be quite large
  offset = 150;
  
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
    %axis( [ minX-offset maxX+offset minY-offset maxY+offset minZ-offset maxZ+offset minZ-offset maxZ+offset ] );
    axis equal
    axis off
  end
  
  grid off;
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  C = strsplit( char( dataStr( 1, dataIndex ) ), '_' );
  title( strcat( C( 1, 1 ), ' Time Step ', num2str(curT) ) );
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
  if renderMovie == 1
    writerObj = VideoWriter( char(movieDir) );
    writerObj.FrameRate = 10;
    writerObj.Quality = 50;
    open( writerObj );
  end
  
  % surface instance
  S = [];
  % line instance
  L = [];
  % PC instance
  P = [];
  
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
    set( L, 'Visible', 'off' );
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
      
      if drawEllipsoids == 1
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
            lMat = getLinkMatrix( p1, p2 );
            
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
          
          % radii of the ellipsoid
          radii = sqrt( diag(D) );
          
          xEigVec = Q(:, 1);
          yEigVec = Q(:, 2);
          zEigVec = Q(:, 3);
          
          % draw the single cell as ellipsoid
          [ x, y, z ] = ellipsoid( 0, 0, 0, radii(1)/2., radii(2)/2., radii(3)/2., nEllip );
          X = p1(1) + x*xEigVec(1) + y*yEigVec(1) + z*zEigVec(1);
          Y = p1(2) + x*xEigVec(2) + y*yEigVec(2) + z*zEigVec(2);
          Z = p1(3) + x*xEigVec(3) + y*yEigVec(3) + z*zEigVec(3);
          
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
            
            % and store the start/end points of the lines in linePos
            linePos = [ linePos ; lineX(1) lineY(1) lineZ(1) ; lineX(2) lineY(2) lineZ(2) ];
            
            % draw line in 3D
            %line( lineX, lineY, lineZ );
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
      end
      
      % after drawing all cells perform a pca for getting the principal axes
      % if at least four points exist (such that three PCs are generated)
      numPCACells = size( cellFileMat( :, 1 ), 1 );
      if numPCACells > 3
        start = mean( cellFileMat );
        %coeff = pca( cellFileMat );
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
        
        % set plane position
        planePos = dir * planePositionFactor;
        
        % draw lines of the major ellipsoid axes onto the projected plane
        l = 1;
        dimL = size( linePos( :, 1 ), 1 ) + 1;
        while l < dimL
          if strcmp( lineStr( 1, lineRenderType ), 'renderOnlyLargestElongation' )
            % check the length of all three main major axes and choose the
            % largest one to render only
            lengthVec = zeros( 3, 1 );
            for ll=1:3
              linep1 = projectOnPlane( linePos(l, :), planePos, u, v );
              linep2 = projectOnPlane( linePos(l+1, :), planePos, u, v );
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
            L(l) = line( lineX, lineY, lineZ, 'Color', [ 0. 0. 0. ], 'LineWidth', 2 );
            hold on;
          else
            linep1 = projectOnPlane( linePos(l, :), planePos, u, v );
            linep2 = projectOnPlane( linePos(l+1, :), planePos, u, v );
            lineX = [ linep1(1), linep2(1) ];
            lineY = [ linep1(2), linep2(2) ];
            lineZ = [ linep1(3), linep2(3) ];
            L(l) = line( lineX, lineY, lineZ, 'Color', [ 0. 0. 0. ], 'LineWidth', 2 );
            hold on;
            l = l+2;
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
        title( strcat( C( 1, 1 ), ' Time Step ', num2str(curT) ) );
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
        continue;
        hold off;
        set( f,'nextplot','replacechildren' );
        if equalAspectRatio == 0
          axis( [ minX-offset maxX+offset minY-offset maxY+offset minZ-offset maxZ+offset ] );
        else
          axis 'vis3d'
        end
        grid off;
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        C = strsplit( char( dataStr( 1, dataIndex ) ), '_' );
        title( strcat( C( 1, 1 ), ' Time Step ', num2str(curT) ) );
        view(2);
      end
      
      % first delete last lightsource
      delete(findall(gcf, 'Type', 'light'))
      camlight headlight;
      
      if renderMovie == 1
        writeVideo( writerObj, getframe(f) );
      end
      
      % TODO: does not work at the moment
      if createImages == 1
        F = getframe(S);
        I = image( F.cdata );
        
        if curT < 10
          digit = strcat( dataStr( 1, dataIndex ), '_00' );
        elseif curT < 100
          digit = strcat( dataStr( 1, dataIndex ), '_0' );
        else
          digit = strcat( dataStr( 1, dataIndex ), '_' );
        end
        
        imwrite(I, strcat( imageDir, digit, num2str(curT), '.png' ), 'png');
      end
    end
  end
  
  if renderMovie == 1
    close(writerObj);
  end
  
  % close figure after processing
  %close(f);
  
end
