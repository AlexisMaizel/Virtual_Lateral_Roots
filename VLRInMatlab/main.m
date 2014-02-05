%%%% The mesh deformation is based on the paper of Graner et al,
%%%% Discrete rearranging disordered patterns, part I: Robust
%%%% statistical tools in two or three dimensions, Eur. Phys. J.,
%%%% pp 349 -- 369, 2008

function main()

%%%%% setting of properties %%%%%%
% start with the current time step
startT = 1;
% draw ellipsoids?
drawEllipsoids = 1;
% draw only cell positions?
drawSpheres = 0;
% draw delaunay tri?
drawDelaunay = 0;
% render movie?
renderMovie = 1;
% create images for each time step
% TODO: does not work at the moment
createImages = 0;
% render only master file?
renderMasterFile = 1;
% default view in positive z direction
dir = [ 0, 90 ];
% camera view:
% 0 -> top
% 1 -> side
% 2 -> radial
% 3 -> 3D
cView = 3;

% path to movie output
movieDir = '/home/necrolyte/Data/VLR/VLRDataForMatlab/121204_3D_310114_master.avi';

% path to image output
imageDir = '/home/necrolyte/Data/VLR/VLRDataForMatlab/121204/TopView/';

% output format of values
format shortG

% reading raw data
fileID = fopen('/home/necrolyte/Uni/LateralRootGrowth/MatlabProject/data/121204_raw_2014.csv');
formatSpec = '%d %f %f %f %d %d %s';
data = textscan( fileID, formatSpec, 'HeaderLines', 1, 'Delimiter', ';' );
fclose(fileID);

% get dimension aka number of lines
col = size(data{1});
numLines = col(1,1);
% store relevant columns in corresponding data structures
IdCol = data{1};
XCol = data{2};
YCol = data{3};
ZCol = data{4};
TCol = data{5};
LCol = data{7};

% get maximum of time steps
maxT = max( TCol );

% master cell file
masterFile = 4;

% manual cell file information at the moment
k = { 219, 220, 191, 234, 7, 181, 1212, 2233, 4, 5566, 8, 9910, 11, 201, 204 };
v = { 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 8, 8 };
cellFileMap = containers.Map( k, v );

% get bounding box are by determining
% min and max of x/y/z values
% loop over all time steps
minX = min( XCol );
minY = min( YCol );
minZ = min( ZCol );
maxX = max( XCol );
maxY = max( YCol );
maxZ = max( ZCol );

% boundary offset
offset = 130;

% window properties
f = figure( 'Name', 'Mesh Deformation', 'Position', [100 100 3000 1000] );
hold on;
cameratoolbar( 'SetMode', 'orbit' );
%zoom on;

if drawSpheres == 1
  axis vis3d
else
  axis( [ minX-offset maxX+offset minY-offset maxY+offset minZ-offset maxZ+offset ] );
end
grid off;
xlabel('X');
ylabel('Y');
zlabel('Z');
title( strcat( 'Time Step ', num2str(startT) ) );
camlight headlight;
view(dir);

% colormap
cm = hsv(8);

% gca is the current axes handle
set( gca,'nextplot','replacechildren' );
% gcf is the current figure handle
%lighting phong
%set( gcf, 'Renderer', 'zbuffer' );
lighting gouraud
set( gcf, 'Renderer', 'OpenGL' );
set( gcf,'nextplot','replacechildren' );

% video output options
if renderMovie == 1
  writerObj = VideoWriter(movieDir);
  writerObj.FrameRate = 10;
  writerObj.Quality = 50;
  open( writerObj );
end

% surface instance
S = [];

% loop over all time steps
for curT=startT:maxT
  % number of cells for the current time step
  numCells = 0;
  
  % matrix of positions for current time step
  matPos = [];
  
  for j=1:numLines
    if TCol(j) == curT %&& ( getLineageId( IdCol(j), LCol(j) ) == 1212 || getLineageId( IdCol(j), LCol(j) ) == 2233 )
      pos = [ XCol(j) YCol(j) ZCol(j) ];
      matPos = [matPos; pos];
      numCells = numCells +1;
    end
  end
  
  % clean figure content by removing the
  % last ellipsoids in the previous time step
  set( S, 'Visible', 'off' );
  %delete(S)
  
  % if at least three cells exists
  if numCells > 3
    % delaunay triangulation visualzation
    tri = delaunayTriangulation( matPos(:,1), matPos(:,2), matPos(:,3) );
    
    % number of total links in current time step
    numTotalLinks = size( edges(tri), 1 );
    
    if drawDelaunay == 1
      tetramesh(tri, 'FaceColor', cc( 1, : ), 'FaceAlpha', 0.9 );
    end
    
    if drawEllipsoids == 1
      % draw an ellipsoid for each cell
      for c=1:numCells
        % get position of current cell
        p1 = [ tri.Points( c, 1 ) tri.Points( c, 2 ) tri.Points( c, 3 ) ];
        
        % get lineage id
        lineageId = getLineageIdFromPos( p1(1), p1(2), p1(3), curT, data );
        
        % and corresponding color
        %color = cm( cellFileMap(lineageId), : );
        if lineageId == 1212 || lineageId == 2233
          color = cm( 1, : );
        else
          color = cm( 2, : );
          if renderMasterFile == 1
            continue;
          end
        end
        
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
        % of neighbors
        M = M./numLinks;
        
        % compute the eigenvectors and eigenvalues of matrix M
        % The columns of Q are the eigenvectors and the diagonal
        % elements of D are the eigenvalues
        [Q,D] = eig(M);
        
        % radii of the ellipsoid
        radii = sqrt( diag(D) );
        
        % determine the angle values between each eigen vector and the
        % coordinate axes
        xAxis = [ 1 0 0 ];
        yAxis = [ 0 1 0 ];
        zAxis = [ 0 0 1 ];
        
        xEigVec = Q(:, 1);
        yEigVec = Q(:, 2);
        zEigVec = Q(:, 3);
        
        % draw the single cell as ellipsoid
        if drawSpheres == 0
          [ x, y, z ] = ellipsoid( 0, 0, 0, radii(1)/2., radii(2)/2., radii(3)/2., 20 );
          X = p1(1) + x*xEigVec(1) + y*yEigVec(1) + z*zEigVec(1);
          Y = p1(2) + x*xEigVec(2) + y*yEigVec(2) + z*zEigVec(2);
          Z = p1(3) + x*xEigVec(3) + y*yEigVec(3) + z*zEigVec(3);
          S(c) = surface( X, Y, Z, 'FaceColor', color, 'EdgeColor', 'none', 'FaceLighting','gouraud' );
          hold on;
        else
          [ x, y, z ] = ellipsoid( p1(1), p1(2), p1(3), 10., 10., 10., 20 );
          S(c) = surface( x, y, z,'FaceColor', color, 'EdgeColor', 'none', 'FaceLighting','gouraud' );
          hold on;
        end
      end
    end
    
    hold off;
    set( f,'nextplot','replacechildren' );
    if drawSpheres == 1
      axis vis3d
    else
      axis( [ minX-offset maxX+offset minY-offset maxY+offset minZ-offset maxZ+offset ] );
      daspect('manual')
    end
    grid off;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title( strcat( 'Time Step ', num2str(curT) ) );
    view(dir);
    % TODO: include optical rotation information into current view
    %rot = [ 2 -5 41 ];
    if cView == 0
      % top view
      camorbit( 0, 80, [ 0 0 1 ] );
      camorbit( 20, 0, [ 0 0 1 ] );
      camorbit( 0, 90, [ 0 0 1 ] );
    elseif cView == 1
      % side view
      camorbit( 0, 80, [ 0 0 1 ] );
      camorbit( 10, 0, [ 0 0 1 ] );
    elseif cView == 2
      % radial view
      camorbit( 0, 80, [ 0 0 1 ] );
      camorbit( 105, 0, [ 0 0 1 ] );
    else
      % 3D view
      view(3);
    end
    
    % first delete last lightsource
    delete(findall(gcf, 'Type', 'light'))
    camlight headlight;
    
    if renderMovie == 1
      writeVideo( writerObj, getframe(f) );
    end
    
    if createImages == 1
      F = getframe(S);
      I = image( F.cdata );
      
      if curT < 10
        digit = '121204_00';
      elseif curT < 100
        digit = '121204_0';
      else
        digit = '121204_';
      end
      
      imwrite(I, strcat( imageDir, digit, num2str(curT), '.png' ), 'png');
    end
  end
end

if renderMovie == 1
  close(writerObj);
end

%close(f);

end
