setWorkingPathProperties()

chosenData = 1;
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };
rawDataStr = { '120830' '121204' '121211' '130508' '130607' };
startT = 1;
endT = 1;

radEllip = 5;
lineWidth = 1;

% vector of view offsets for x/Y/Z coordinates
viewOffsets = [ 100 100 100 250 250 250 ];

% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
cView = 2;

% use triangulation based on delaunay or alpha shape
% 1 -> convex hull
% 2 -> alpha shape
triangulationType = 1;

drawDelaunay = 1;
drawVoronoi = 0;

% get the alpha shape radii for all time steps
if triangulationType == 2
  alphaRadiiVector = getAlphaRadius( dataStr( 1, chosenData ) );
end

% output format of values
format longG

[ cellData, dimData, centerPosPerTimeStep, numCellsPerTimeStep,...
  minX, minY, minZ, maxX, maxY, maxZ ] = readRawData( dataStr( 1, chosenData ) );
  
% figure properties
f = figure( 'Name', 'Graph Topology Analysis', 'Position', [100 100 800 800] );
% activate orbit rotation by default
cameratoolbar( 'SetMode', 'orbit' );
% activate none coord system by default for not resetting the camera up
% vector when rotating
cameratoolbar( 'SetCoordSys', 'none' );
% show camera toolbar by default
cameratoolbar( 'Show' );

% set colormap and colorbar depending on the number of cell files
cm = hsv( 10 );
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

% surface instance
S = [];
% ellipse instance
ELLIP = [];
ELLIPPATCH = [];

% get stored eigenvectors for the last time step to set the same
% direction view for each time step
%coeff = getPrincipalComponents( dataStr( 1, chosenData ), 1 );

% set PC depending on the viewing direction
% if cView == 1
%   dir = coeff(:,2);
%   u = coeff(:,1);
%   v = coeff(:,3);
% elseif cView == 2
%   dir = coeff(:,3);
%   u = coeff(:,1);
%   v = coeff(:,2);
% elseif cView == 3
%   dir = -coeff(:,1);
%   u = -coeff(:,3);
%   v = coeff(:,2);
% end
%   
% % set plane position
% planePos = dir * 1;
% 
% plane = [ planePos(1) planePos(2) planePos(3)...
%   u(1) u(2) u(3)...
%   v(1) v(2) v(3) ];
% TF = createBasisTransform3d( 'g', plane );

% set axes properties
% minAxes = projectOnPlane( [ minX minY minZ ], planePos, u, v );
% minAxes = transformPoint3d( minAxes, TF );
% maxAxes = projectOnPlane( [ maxX maxY maxZ ], planePos, u, v );
% maxAxes = transformPoint3d( maxAxes, TF );

% after transformation the individual coords of min and max
% may be switched
% for mm=1:3
%   if minAxes(mm) > maxAxes(mm)
%     tmp = minAxes(mm);
%     minAxes(mm) = maxAxes(mm);
%     maxAxes(mm) = tmp;
%   end
% end

nucleiCounter = 1;
% loop over all time steps
for curT=startT:endT
  % clean figure content by removing the
  % last ellipsoids and lines in the previous time step
  hideHandle( S );
  hideHandle( ELLIP );
  hideHandle( ELLIPPATCH );
  
  hold on;
%   viewOffset = viewOffsets( chosenData );
%   axis( [ minAxes(1)-viewOffset maxAxes(1)+viewOffset...
%     minAxes(2)-viewOffset maxAxes(2)+viewOffset...
%     minAxes(3)-viewOffset maxAxes(3)+viewOffset...
%     minAxes(3)-viewOffset maxAxes(3)+viewOffset ] );
  axis on
  daspect( [ 1 1 1 ] );
  
  grid off;
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  C = strsplit( char( dataStr( 1, chosenData ) ), '_' );
  title( strcat( C( 1, 1 ), ' Time Step ', num2str(curT) ) );
  
  % number of cells in the current time step
  numCells = numCellsPerTimeStep( curT, 1 );
  
  % matrix of positions for current time step
  matPos = zeros( numCells, 3 );
  
  cellCounter = 1;
  for j=1:dimData
    if cellData{ j, 5 } == curT
      pos = [ cellData{ j, 2 } cellData{ j, 3 } cellData{ j, 4 } ];
      matPos( cellCounter, : ) = pos;
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
    
    if drawDelaunay == 1 && triangulationType == 1
      tetramesh(tri, 'FaceColor', cm( 2, : ), 'FaceAlpha', 0.0 );
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
    
    % draw an ellipsoid for each cell
    for c=1:numCells
      % get position of current cell
      if triangulationType == 1
        p = [ tri.Points( c, 1 ) tri.Points( c, 2 ) tri.Points( c, 3 ) ];
      else
        p = [ matPos( c, 1 ) matPos( c, 2 ) matPos( c, 3 ) ];
      end
      
      color = cm( 1, : );
      [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
        drawEllipse3d( p(1), p(2), p(3), radEllip, radEllip, 0, 0 );
      set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
      set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
      nucleiCounter = nucleiCounter+1;
    end
    
    hold off;
    set( f,'nextplot','replacechildren' );
  end
end
  