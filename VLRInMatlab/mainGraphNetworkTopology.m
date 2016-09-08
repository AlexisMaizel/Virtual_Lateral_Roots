setWorkingPathProperties()

chosenData = 5;
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };
rawDataStr = { '120830' '121204' '121211' '130508' '130607' };
startT = 1;
endT = 300;

radEllip = 5;
lineWidth = 1;

% vector of view offsets for x/Y/Z coordinates
viewOffsets = [ 100 100 100 250 250 250 ];

% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
cView = 2;
viewStr = { 'Top' 'Side' 'Radial' };

% use triangulation based on delaunay or alpha shape
% 1 -> convex hull
% 2 -> alpha shape
triangulationType = 1;

% show specific cell file in [-3,3]; only valid for 3D graph vis
cellFileType = 0;

% drawing parameters
draw3DGraph = 0;
drawDistributions = 1;
drawTriangulation = 0;
drawVoronoi = 0;

numBins = 20;

% number of generated graph properties each displayed in a different
% subplot
numGraphProperties = 4;
graphPropString = { 'Degree' 'Closeness' 'Betweenness' 'PageRank' };

% get the alpha shape radii for all time steps
if triangulationType == 2
  alphaRadiiVector = getAlphaRadius( dataStr( 1, chosenData ) );
end

% output format of values
format longG

[ cellData, dimData, centerPosPerTimeStep, numCellsPerTimeStep,...
  minX, minY, minZ, maxX, maxY, maxZ, cellFileMap, minCF, maxCF ] =...
  readRawData( dataStr( 1, chosenData ) );

% figure properties
f = figure( 'Name', 'Graph Topology Analysis', 'Position', [10 10 1400 800] );
% activate orbit rotation by default
cameratoolbar( 'SetMode', 'orbit' );
% activate none coord system by default for not resetting the camera up
% vector when rotating
cameratoolbar( 'SetCoordSys', 'none' );
% show camera toolbar by default
cameratoolbar( 'Show' );

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
% ellipse instance
ELLIP = [];
ELLIPPATCH = [];

% get stored eigenvectors for the last time step to set the same
% direction view for each time step
coeff = getPrincipalComponents( dataStr( 1, chosenData ), 1 );

% set PC depending on the viewing direction
if cView == 1
  dir = coeff(:,2);
  u = coeff(:,1);
  v = coeff(:,3);
elseif cView == 2
  dir = coeff(:,3);
  u = coeff(:,1);
  v = coeff(:,2);
  % invert pc to get the correct position of the dome
  if chosenData == 2 || chosenData == 3
    v = v .* -1;
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

% set axes properties
minAxes = projectOnPlane( [ minX minY minZ ], planePos, u, v );
minAxes = transformPoint3d( minAxes, TF );
maxAxes = projectOnPlane( [ maxX maxY maxZ ], planePos, u, v );
maxAxes = transformPoint3d( maxAxes, TF );

% after transformation the individual coords of min and max
% may be switched
for mm=1:3
  if minAxes(mm) > maxAxes(mm)
    tmp = minAxes(mm);
    minAxes(mm) = maxAxes(mm);
    maxAxes(mm) = tmp;
  end
end

% set the possible cell files such that the colors are fixed for each data
% set
minCF = -4;
maxCF = 3;

imageOutputPath = strcat( 'I:\GraphTopologyAnalysis\', rawDataStr( 1, chosenData ), '\' );
mkdir( char(imageOutputPath) );

nucleiCounter = 1;
% loop over all time steps
for curT=startT:49:endT
  if curT < 10
    digit = '00';
  elseif curT < 100
    digit = '0';
  else
    digit = '';
  end
  % clean figure content by removing the
  % last ellipsoids and lines in the previous time step
  hideHandle( S );
  hideHandle( ELLIP );
  hideHandle( ELLIPPATCH );
  
  for pl=1:numGraphProperties
    subplot( numGraphProperties/2, numGraphProperties/2, pl, 'replace' )
    viewOffset = viewOffsets( chosenData );
    if drawDistributions == 0 && draw3DGraph == 1
      axis( [ minAxes(1)-viewOffset maxAxes(1)+viewOffset...
        minAxes(2)-viewOffset maxAxes(2)+viewOffset...
        minAxes(3)-viewOffset maxAxes(3)+viewOffset...
        minAxes(3)-viewOffset maxAxes(3)+viewOffset ] );
    end
    axis on
    if draw3DGraph == 1 && drawDistributions == 0
      daspect( [ 1 1 1 ] );
    end
    
    grid off;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    C = strsplit( char( dataStr( 1, chosenData ) ), '_' );
    title( strcat( C( 1, 1 ), ' Time Step ', num2str(curT), '\_', char( graphPropString( 1, pl ) ) ) );
  end
  
  % number of cells in the current time step
  numCells = numCellsPerTimeStep( curT, 1 )
  
  % matrix of positions for current time step
  matPos = zeros( numCells, 3 );
  
  % cell file mapper
  idToCF = zeros( numCells, 1 );
  
  cellCounter = 1;
  for j=1:dimData
    if cellData{ j, 5 } == curT
      pos = [ cellData{ j, 2 } cellData{ j, 3 } cellData{ j, 4 } ];
      matPos( cellCounter, : ) = pos;
      idToCF( cellCounter, 1 ) = cellData{ j, 1 };
      cellCounter = cellCounter + 1;
    end
  end
  
  % create empty graph matrix representing the connectivity information
  connMat = zeros( numCells, numCells );
  
  % if at least three cells exists
  if numCells > 3
    if triangulationType == 1
      % delaunay triangulation
      tri = delaunayTriangulation( matPos(:,1), matPos(:,2), matPos(:,3) );
      % loop over all edges
      E = edges(tri);
      for e=1:size( E, 1 )
        connMat( E( e, 1 ), E( e, 2 ) ) = 1;
        connMat( E( e, 2 ), E( e, 1 ) ) = 1;
      end
    else
      % alpha shape triangulation
      [Vol,Shape] = alphavol( [ matPos(:,1) matPos(:,2) matPos(:,3) ], sqrt( alphaRadiiVector( curT, 1 )) );
      tri = Shape.tri;
      % loop over all boundary facets
      for f=1:size( Shape.bnd, 1 )
        % TODO: symmetric case
        connMat( tri( f, 1 ), tri( f, 2 ) ) = 1;
        connMat( tri( f, 2 ), tri( f, 3 ) ) = 1;
        connMat( tri( f, 3 ), tri( f, 1 ) ) = 1;
      end
    end
    
    % compute graph properties
    G = graph( connMat );
    %p = plot( G, 'MarkerSize',5 );
    
    for pl=1:numGraphProperties
      subplot( numGraphProperties/2, numGraphProperties/2, pl )
      
      if pl == 1
        gp = centrality( G, 'degree' );
      elseif pl == 2
        gp = centrality( G, 'closeness' );
      elseif pl == 3
        gp = centrality( G, 'betweenness' );
      elseif pl == 4
        gp = centrality( G, 'pagerank' );
      end
      
      cm = jet( numCells );
      colormap( cm );
      minNC = min(gp);
      maxNC = max(gp);
      %colorbar( 'Ytick', [0, numCells/2, numCells], 'Yticklabel', {minNC, maxNC} );
      y = floor( ((gp-minNC)/range( gp )) * (numCells-1) ) + 1;
      
      if drawTriangulation == 1 && triangulationType == 1
        tetramesh(tri, 'FaceColor', cm( 2, : ), 'FaceAlpha', 0.0, 'EdgeAlpha', 0.1 );
      end
      
      numCellFiles = maxCF-minCF+1;
      cFToGP = cell( numCellFiles, 1 );
      % draw an ellipsoid for each cell
      hold on
      for c=1:numCells
        % get position of current cell
        if triangulationType == 1
          p = [ tri.Points( c, 1 ) tri.Points( c, 2 ) tri.Points( c, 3 ) ];
        else
          p = [ matPos( c, 1 ) matPos( c, 2 ) matPos( c, 3 ) ];
        end
        
        if draw3DGraph == 1 && drawDistributions == 0
          % transform points to fit into viewing types
          p = projectOnPlane( p, planePos, u, v );
          p = transformPoint3d( p, TF );
        end
        
        % separate the gp values into each cell file
        pureCellFile = cellFileMap( idToCF( c, 1 ) );
        % shift to start
        cf = pureCellFile - minCF + 1;
        cFToGP{ cf, 1 } = [ cFToGP{ cf, 1 } gp( c, 1 ) ];
        
        if draw3DGraph == 1 && cellFileType == pureCellFile
          color = cm( y(c, 1), : );
          [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
            drawEllipse3d( p(1), p(2), p(3), radEllip, radEllip, 0, 0 );
          set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
          set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
          nucleiCounter = nucleiCounter+1;
        end
      end
      
      if drawDistributions == 0 && draw3DGraph == 1
        colorbar
        caxis( [minNC, maxNC] )
      end
      
      if drawDistributions == 1
        drawCellFileBar( minNC, maxNC, numCellFiles, numBins, cFToGP,...
          pl, minCF, maxCF );
      end
      
      % TODO
      %     if drawVoronoi == 1 && triangulationType == 1
      %       [ V, R ] = voronoiDiagram( tri );
      %       for i=1:size( R, 1 )
      %         indVert = R{i};
      %         vertices = V( indVert, : );
      %         vertices( 1, : ) = [];
      %         if size( vertices( :, 1 ), 1 ) > 3
      %           %K = convhulln( vertices );
      %           k = boundary( vertices );
      %           patch( vertices(k, 1), vertices(k, 2), vertices(k, 3) );
      %         end
      %       end
      %     end
      
      %hold off;
      %set( f,'nextplot','replacechildren' );
    end
  end
  
  if draw3DGraph == 1
    suffixStr = strcat( '_Data_', viewStr( 1, cView ), '_CF', num2str(cellFileType) );
  else
    suffixStr = strcat( '_Distribution_', viewStr( 1, cView ), '_CF', num2str(cellFileType) );
  end
  
  filePath = strcat( imageOutputPath, rawDataStr( 1, chosenData ), '_T', digit, num2str(curT), suffixStr, '.png' );
  export_fig( gcf, char(filePath), '-m2', '-png' );
end
