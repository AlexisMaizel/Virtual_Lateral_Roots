%%%% The mesh deformation is based on the paper of Graner et al,
%%%% Discrete rearranging disordered patterns, part I: Robust
%%%% statistical tools in two or three dimensions, Eur. Phys. J.,
%%%% pp 349 -- 369, 2008

function mainInterpolated()

%%%%% setting of properties %%%%%%
% start with the current time step
startT = 300;
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
renderSingleCellFile = 0;
% master cell file
singleCellFile = 3;
% default view in positive z direction
dir = [ 0, 90 ];
% camera view:
% 0 -> top
% 1 -> side
% 2 -> radial
% 3 -> 3D
cView = 1;

% path to movie output
movieDir = '/home/necrolyte/Data/VLR/VLRDataForMatlab/test_3D_030214_all.avi';

% path to image output
imageDir = '/home/necrolyte/Data/VLR/VLRDataForMatlab/121204/TopView/';

% output format of values
format shortG

% reading raw data
fileID = fopen( '/home/necrolyte/Data/VLR/Virtual_Lateral_Roots/FinalVLRForMatlab/130607_raw.csv');
% format of data sets:
% ObjectID X Y Z Timepoint Radius Precursors Color Lineage TrackID TrackColor TrackGroup Layer DivisionType
formatSpec = '%d %f %f %f %d %d %q %q %d %q %q %d %d %q';
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
% z coord
ZCol = data{4};
% Time Step
TCol = data{5};
% Lineage Tree Id
LCol = data{9};
% Cell File Id
CFCol = data{12};
% Cell Layer Id
CLCol = data{13};
% Divison Type
DCol = data{14};

% get maximum of time steps
maxT = max( TCol );

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

% boundary offset
offset = 130;

% window properties
f = figure( 'Name', 'Mesh Deformation', 'Position', [100 100 800 800] );
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

% number of cell files
numCellFiles = double( maxCF-minCF+1 );
ticks = [ minCF:maxCF ];
% set colormap depending on the number of cell files
cm = hsv( numCellFiles );
colormap( hsv( numCellFiles ) );
h = colorbar;
set( h, 'YTickMode', 'manual' );
set( h, 'YTickLabel', ticks );

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
  % last ellipsoids in the previous time step
  set( S, 'Visible', 'off' );
  %delete(S)
  
  % if at least three cells exists
  if numCells > 3
    % delaunay triangulation visualzation
    tri = delaunayTriangulation( matPos(:,1), matPos(:,2), matPos(:,3) );
    
    % number of total links in current time step
    %numTotalLinks = size( edges(tri), 1 );
    
    if drawDelaunay == 1
      tetramesh(tri, 'FaceColor', cc( 1, : ), 'FaceAlpha', 0.9 );
    end
    
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
