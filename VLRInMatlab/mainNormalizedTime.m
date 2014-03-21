%%%% The mesh deformation is based on the paper of Graner et al,
%%%% Discrete rearranging disordered patterns, part I: Robust
%%%% statistical tools in two or three dimensions, Eur. Phys. J.,
%%%% pp 349 -- 369, 2008

addpath( '/home/necrolyte/Data/VLR/Virtual_Lateral_Roots/VLRInMatlab/geom3d' );

%%%%% setting of properties %%%%%%
% color range for different data sets
% { 'm' 'k' 'r' 'g' 'b' 'c' };
colors = [ [ 1 0 1 ]; [ 0 0 0 ]; [ 1 0 0 ]; [ 0 1 0 ]; [ 0 0 1 ] ];
% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
cView = 2;
% decide which term should be included in the time evolution
renderTermType = 3;
termTypeStr = { 'B' 'T' 'All' };
% startIndex
startI = 1;
% endIndex
endI = 20;
% deltaT value based on the paper mentioned above
% have to be a divider of the max index
deltaI = 1;
% draw delaunay tri?
drawDelaunay = 0;
% if set to one then only a single time step
% is rendered given by startT
exportType = 2;
% vector of data strings
exportTypeStr = { 'SingleFigure' 'AsImages' };
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
% render contour (=convexHull) instead of ellipses
visualizationType = { 'Ellipsoids' 'Ellipses' 'Contour' };
visType = 2;
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
triangulationType = 2;
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
endData = 1;
% num of data
numData = 5;
% number of cell files
numCellFiles = 6;%double( maxCF-minCF+1 );

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

% apply preprocessing step of data
[ cellDatas, dimData, maxT, numCellsPerTimeStep, centerPosPerTimeStep, totalMinAxes, totalMaxAxes, cellFileMap ] =...
  prepareData( dataStr, startData, endData, numData, visualizationType( 1, visType ), renderSingleCellFile, cView );

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

% this variable is set to zero after the first traversal of the loop
% since it indicated the begin of the loop; this is required in order to
% set the information of the next time step to the infos of the current
% one in each time step which prevent redundant computations, but this is
% valid for all loop traversals except the first one
begin = 1;

% the variables of two successive time steps are added a 'C' (Current)
% or 'N' (Next) at the end
% number of cells for the current and next time step
numCellsC = zeros( numData, 1 );
numCellsN = zeros( numData, 1 );

% matrix of positions for current and next time step
matPosC = cell( numData );
matPosN = cell( numData );

% vector of object ids in order to access the cell
% file information in the cell file map
cellIdsC = cell( numData );
cellIdsN = cell( numData );

% vector of strings storing the precursors
cellPrecursorsN = cell( numData );

% delaunay triangulation
triC = cell( numData );
triN = cell( numData );

% number of current time step for each data
curTC = zeros( numData, 1 );
curTN = zeros( numData, 1 );

% unique edges in the triangulation
uniqueEdgesC = cell( numData );
uniqueEdgesN = cell( numData );

% number of total links in current and next time step
numTotalLinksC = zeros( numData, 1 );
numTotalLinksN = zeros( numData, 1 );
numTotalAveragedLinks = zeros( numData, 1 );

% loop over all normalized steps
curI=startI;
while curI < endI-deltaI+1
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
    
    % NEW: only continue with this time step that is synchronized in
    % number of cells for each data set
    if begin ~= 1
      numNormCellsN = getNormalizedCellNumber( curI+deltaI, 18, 143, 1, 20 );
      curTC(dataIndex) = curTN(dataIndex);
      curTN(dataIndex) = 1;
      for j=1:maxT(dataIndex)
        if numNormCellsN - epsilon < numCellsPerTimeStep{dataIndex}(j,1)...
            && numNormCellsN + epsilon > numCellsPerTimeStep{dataIndex}(j,1)
          curTN(dataIndex) = j;
          break;
        end
      end
    else
      numNormCellsC = getNormalizedCellNumber( curI, 18, 143, 1, 20 );
      numNormCellsN = getNormalizedCellNumber( curI+deltaI, 18, 143, 1, 20 );
      curTC(dataIndex) = 1;
      for j=1:maxT(dataIndex)
        if numNormCellsC - epsilon < numCellsPerTimeStep{dataIndex}(j,1)...
            && numNormCellsC + epsilon > numCellsPerTimeStep{dataIndex}(j,1)
          curTC(dataIndex) = j;
          break;
        end
      end
      curTN(dataIndex) = 1;
      for j=1:maxT(dataIndex)
        if numNormCellsN - epsilon < numCellsPerTimeStep{dataIndex}(j,1)...
            && numNormCellsN + epsilon > numCellsPerTimeStep{dataIndex}(j,1)
          curTN(dataIndex) = j;
          break;
        end
      end
    end
    
    % update cell information
    if begin ~= 1
      % number of cells for current and next time step and data set
      numCellsC(dataIndex) = numCellsN(dataIndex);
      numCellsN(dataIndex) = numCellsPerTimeStep{dataIndex}(curTN(dataIndex), 1);
      matPosC{dataIndex} = matPosN{dataIndex};
      matPosN{dataIndex} = [];
      cellIdsC{dataIndex} = cellIdsN{dataIndex};
      cellIdsN{dataIndex} = [];
      
      % matrix of positions for current time step
      matPos2 = zeros( numCellsN(dataIndex), 3 );
      cellIds2 = zeros( numCellsN(dataIndex), 3 );
      cellPrecursors2 = cell( numCellsN(dataIndex),1 );
      
      nc2 = 1;
      for j=1:dimData( dataIndex )
        if cellDatas{ dataIndex }{ j, 5 } == curTN(dataIndex)
          pos = [ cellDatas{ dataIndex }{ j, 2 } cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
          matPos2(nc2, :) = pos;
          cellIds2(nc2, :) = cellDatas{ dataIndex }{ j, 1 };
          cellPrecursors2{nc2} = cellDatas{ dataIndex }{ j, 8 };
          nc2 = nc2 + 1;
        end
      end
    else
      % number of cells for current and next time step and data set
      numCellsC(dataIndex) = numCellsPerTimeStep{dataIndex}(curTC(dataIndex), 1);
      numCellsN(dataIndex) = numCellsPerTimeStep{dataIndex}(curTN(dataIndex), 1);
      
      % matrix of positions for current time step
      matPos1 = zeros( numCellsC(dataIndex), 3 );
      matPos2 = zeros( numCellsN(dataIndex), 3 );
      cellIds1 = zeros( numCellsC(dataIndex), 3 );
      cellIds2 = zeros( numCellsN(dataIndex), 3 );
      cellPrecursors2 = cell( numCellsN(dataIndex),1 );
      
      nc1 = 1;
      nc2 = 1;
      for j=1:dimData( dataIndex )
        if cellDatas{ dataIndex }{ j, 5 } == curTC(dataIndex)
          pos = [ cellDatas{ dataIndex }{ j, 2 } cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
          matPos1(nc1, :) = pos;
          cellIds1(nc1, :) = cellDatas{ dataIndex }{ j, 1 };
          nc1 = nc1 + 1;
        elseif cellDatas{ dataIndex }{ j, 5 } == curTN(dataIndex)
          pos = [ cellDatas{ dataIndex }{ j, 2 } cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
          matPos2(nc2, :) = pos;
          cellIds2(nc2, :) = cellDatas{ dataIndex }{ j, 1 };
          cellPrecursors2{nc2} = cellDatas{ dataIndex }{ j, 8 };
          nc2 = nc2 + 1;
        end
      end
      
      matPosC{dataIndex} = matPos1;
      matPosN{dataIndex} = matPos2;
      cellIdsC{dataIndex} = cellIds1;
      cellIdsN{dataIndex} = cellIds2;
      cellPrecursorsN{dataIndex} = cellPrecursors2;
    end
    
    % if at least three cells exists
    if numCellsC(dataIndex) > 3 && numCellsN(dataIndex) > 3
      if begin ~= 1
        uniqueEdgesC{dataIndex} = uniqueEdgesN{dataIndex};
        triC{dataIndex} = triN{dataIndex};
        if triangulationType == 1
          % delaunay triangulation
          triN{dataIndex} = delaunayTriangulation( matPosN{dataIndex}(:,1), matPosN{dataIndex}(:,2), matPosN{dataIndex}(:,3) );
          uniqueEdgesN{dataIndex} = edges(triN{dataIndex});
        else
          % alpha shape triangulation
          [VolN,ShapeN] = alphavol( [ matPosN{dataIndex}(:,1) matPosN{dataIndex}(:,2) matPosN{dataIndex}(:,3) ], sqrt( alphaRadiiVector( curTN(dataIndex), 1 )) );
          triN{dataIndex} = ShapeN.tri;
          uniqueEdgesN{dataIndex} = getUniqueEdges( triN{dataIndex} );
        end
        numTotalLinksC(dataIndex) = numTotalLinksN(dataIndex);
        numTotalLinksN(dataIndex) = size( uniqueEdgesN{dataIndex}, 1 );
        numTotalAveragedLinks(dataIndex) = ( numTotalLinksC(dataIndex) + numTotalLinksN(dataIndex) )/2.;
      else
        if triangulationType == 1
          % delaunay triangulation
          triC{dataIndex} = delaunayTriangulation( matPosC{dataIndex}(:,1), matPosC{dataIndex}(:,2), matPosC{dataIndex}(:,3) );
          uniqueEdgesC{dataIndex} = edges( triC{dataIndex} );
          triN{dataIndex} = delaunayTriangulation( matPosN{dataIndex}(:,1), matPosN{dataIndex}(:,2), matPosN{dataIndex}(:,3) );
          uniqueEdgesN{dataIndex} = edges( triN{dataIndex} );
        else
          % alpha shape triangulation
          [VolC,ShapeC] = alphavol( [ matPosC{dataIndex}(:,1) matPosC{dataIndex}(:,2) matPosC{dataIndex}(:,3) ], sqrt( alphaRadiiVector( curTC(dataIndex), 1 )) );
          triC{dataIndex} = ShapeC.tri;
          uniqueEdgesC{dataIndex} = getUniqueEdges( triC{dataIndex} );
          [VolN,ShapeN] = alphavol( [ matPosN{dataIndex}(:,1) matPosN{dataIndex}(:,2) matPosN{dataIndex}(:,3) ], sqrt( alphaRadiiVector( curTN(dataIndex), 1 )) );
          triN{dataIndex} = ShapeN.tri;
          uniqueEdgesN{dataIndex} = getUniqueEdges( triN{dataIndex} );
        end
        numTotalLinksC(dataIndex) = size( uniqueEdgesC, 1 );
        numTotalLinksN(dataIndex) = size( uniqueEdgesN, 1 );
        numTotalAveragedLinks(dataIndex) = ( numTotalLinksC(dataIndex) + numTotalLinksN(dataIndex) )/2.;
      end
      
      % determine deltaT
      deltaT = curTN(dataIndex) - curTC(dataIndex);
      
      % compute time evolution for current deltaT and time step
      [cellsInFileMat, lineColorIndex, linePos, minMaxEigenValueIndex,...
        positiveEigenvalueVector, minMaxSemiAxisVector, centerEllipse, indexColorSet ]...
        = computeTimeEvolution( uniqueEdgesC{dataIndex}, uniqueEdgesN{dataIndex},...
        cellIdsC{dataIndex}, cellIdsN{dataIndex}, numCellsN(dataIndex),...
        triC{dataIndex}, triN{dataIndex}, matPosC{dataIndex},...
        matPosN{dataIndex}, cellPrecursorsN{dataIndex}, triangulationType,...
        termTypeStr( 1, renderTermType ), dataStr( 1, dataIndex ),...
        planePos, u, v, TF, deltaT,...
        renderSingleCellFile, singleCellFile, cellFileMap{dataIndex} );
      
      
      
%       
%       % compute center of data
%       centerData = [ 0 0 0 ];
%       numConsideredCells = 0;
%       for c=1:numCells
%         % get position of current cell
%         if triangulationType == 1
%           p1 = [ tri.Points( c, 1 ) tri.Points( c, 2 ) tri.Points( c, 3 ) ];
%         else
%           p1 = [ matPos( c, 1 ) matPos( c, 2 ) matPos( c, 3 ) ];
%         end
%         
%         % get color for ellipsoid
%         if renderSingleCellFile == 1
%           if singleCellFile ~= cellFileMap{ dataIndex }( idToCF(c,1) )
%             continue;
%           else
%             numConsideredCells = numConsideredCells + 1;
%           end
%         else
%           numConsideredCells = numConsideredCells + 1;
%         end
%         
%          centerData = centerData + p1;
%       end
%       
%       if numConsideredCells > 0
%         centerData = centerData./numConsideredCells;
%       end
%       
%       cellFileMat = zeros( numConsideredCells, 3 );
%       linePos = zeros( numConsideredCells, 18 );
%       minMaxSemiAxisVector = zeros( numConsideredCells, 6 );
%       centerEllipse = zeros( numConsideredCells, 3 );
%       
%       nc = 1;
%       % draw an ellipsoid for each cell
%       for c=1:numCells
%         % get position of current cell
%         if triangulationType == 1
%           p1 = [ tri.Points( c, 1 ) tri.Points( c, 2 ) tri.Points( c, 3 ) ];
%         else
%           p1 = [ matPos( c, 1 ) matPos( c, 2 ) matPos( c, 3 ) ];
%         end
%         
%         p1 = p1 - centerData;
%         
%         % get color for ellipsoid
%         if renderSingleCellFile == 1
%           if singleCellFile ~= cellFileMap{ dataIndex }( idToCF(c,1) )
%             continue;
%           end
%         end
%         
%         cellFileMat(nc, :) = [ p1(1) p1(2) p1(3) ];
%         
%         % get neighbor for specific vertex ID = c
%         if triangulationType == 1
%           nVec = getNeighbors( c, tri, numCells );
%         else
%           nVec = getAlphaShapeNeighbors( c, tri );
%         end
% 
%         if size( nVec, 1 ) == 0
%           continue;
%         end
%         
%         % compute the static link matrix
%         M = computeStaticLink( p1, nVec, centerData, tri, matPos, triangulationType );
%         
%         % compute the eigenvectors and eigenvalues of matrix M
%         % The columns of Q are the eigenvectors and the diagonal
%         % elements of D are the eigenvalues
%         [Q,D] = eig(M);
%         
%         % check if the eigenvalues are smaller then zero; if so, then do
%         % not draw a line and consider the absolute value of it -> TODO
%         positiveEigenvalue = [ 1 ; 1 ; 1 ];
%         radii = diag(D);
%         for e=1:3
%           if radii(e) < 0.
%             radii(e) = -radii(e);
%             positiveEigenvalue( e, 1 ) = 0;
%           end
%         end
%         
%         % radii of the ellipsoid
%         radii = sqrt( radii );
%         
%         xEigVec = Q(:, 1);
%         yEigVec = Q(:, 2);
%         zEigVec = Q(:, 3);
%         
%         % draw the single cell as ellipsoid
%         [ x, y, z ] = ellipsoid( 0, 0, 0, radii(1)/2., radii(2)/2., radii(3)/2., nEllip );
%         X = p1(1) + x*xEigVec(1) + y*yEigVec(1) + z*zEigVec(1);
%         Y = p1(2) + x*xEigVec(2) + y*yEigVec(2) + z*zEigVec(2);
%         Z = p1(3) + x*xEigVec(3) + y*yEigVec(3) + z*zEigVec(3);
%         
%         if strcmp( visualizationType( 1, visType ), 'Ellipsoids' )
%           S(c) = surface( X, Y, Z, 'FaceColor', colors( dataIndex, : ),...
%             'EdgeColor', 'none', 'EdgeAlpha', 0,...
%             'FaceLighting', 'gouraud' );
%         end
%         
%         semiLines = zeros( 1, 18 );
%         % draw the three major axes in the origin which are then
%         % rotated according to the eigenvectors
%         for l=1:3
%           if l == 1
%             sX = [ -radii(1)/2., radii(1)/2. ];
%             sY = [ 0, 0 ];
%             sZ = [ 0, 0 ];
%           elseif l == 2
%             sX = [ 0, 0 ];
%             sY = [ -radii(2)/2., radii(2)/2. ];
%             sZ = [ 0, 0 ];
%           else
%             sX = [ 0, 0 ];
%             sY = [ 0, 0 ];
%             sZ = [ -radii(3)/2., radii(3)/2. ];
%           end
%           lineX = p1(1) + sX*xEigVec(1) + sY*yEigVec(1) + sZ*zEigVec(1);
%           lineY = p1(2) + sX*xEigVec(2) + sY*yEigVec(2) + sZ*zEigVec(2);
%           lineZ = p1(3) + sX*xEigVec(3) + sY*yEigVec(3) + sZ*zEigVec(3);
%           
%           projLine1 = applyTransformations( [ lineX(1) lineY(1) lineZ(1) ], planePos, u, v, TF, dataStr( 1, dataIndex ) );
%           projLine2 = applyTransformations( [ lineX(2) lineY(2) lineZ(2) ], planePos, u, v, TF, dataStr( 1, dataIndex ) );
%           
%           % and store the start/end points of the lines in linePos
%           index = 6*(l-1) + 1;
%           semiLines( 1, index:index+2 ) = projLine1;
%           semiLines( 1, index+3:index+5 ) = projLine2;
%         end
%         
%         % update semi axes in 3D
%         linePos(nc, :) = semiLines;
%         
%         % project each vertex of the ellipsoid onto the plane
%         dimP = size( X, 1 );
%         for q=1:dimP
%           for p=1:dimP
%             curPos = applyTransformations( [ X(p,q) Y(p,q) Z(p,q) ], planePos, u, v, TF, dataStr( 1, dataIndex ) );
%             X(p,q) = curPos(1);
%             Y(p,q) = curPos(2);
%             Z(p,q) = curPos(3);
%           end
%         end
%         
% %         PS(c) = surface( X, Y, Z, 'FaceColor', colors( dataIndex, : ),...
% %           'EdgeColor', 'none', 'EdgeAlpha', 0,...
% %           'FaceLighting', 'none' );
%         
%         p1 = applyTransformations( p1, planePos, u, v, TF, dataStr( 1, dataIndex ) );
%         centerEllipse(nc, :) = p1;
%         % the direction is now the normal of the x-y plane
%         minMaxS = determineAxes( X, Y, Z, p1, [ 0 0 1 ] );
%         minMaxSemiAxisVector(nc, :) = minMaxS;
%         nc = nc + 1;
%       end
%       
%       % draw principal components
%       if renderPrincipalComponents == 1
%         start = mean( cellFileMat );
%         %coeff = pca( cellFileMat )
%         if strcmp( visualizationType( 1, visType ), 'Ellipses' ) ||...
%             strcmp( visualizationType( 1, visType ), 'Contour' )
%           start = applyTransformations( start, planePos, u, v, TF, dataStr( 1, dataIndex ) );
%         end
%         arrowLength = 150;
%         for a=1:3
%           if strcmp( visualizationType( 1, visType ), 'Ellipses' ) ||...
%               strcmp( visualizationType( 1, visType ), 'Contour' )
%             endP = applyTransformations( coeff(:, a), planePos, u, v, TF, dataStr( 1, dataIndex ) );
%           elseif strcmp( visualizationType( 1, visType ), 'Ellipsoids' )
%             endP = coeff(:, a);
%           end
%           hold on;
%           P(a) = quiver3( start(1), start(2), start(3),...
%             endP(1), endP(2), endP(3),...
%             arrowLength, 'LineWidth', 3,...
%             'Color', cm( a, : ) );
%         end
%       end
%       
%       if overlapping == 1
%         zOffset = 0.5;
%       else
%         zOffset = 0.;
%       end
%       
%       if strcmp( visualizationType( 1, visType ), 'Ellipses' )
%         dimL = size( centerEllipse, 1 );
%         % draw ellipses and lines
%         for l=1:dimL
%           minSemiPoint = [minMaxSemiAxisVector( l, 1 )...
%             minMaxSemiAxisVector( l, 2 )...
%             minMaxSemiAxisVector( l, 3 )];
%           maxSemiPoint = [minMaxSemiAxisVector( l, 4 )...
%             minMaxSemiAxisVector( l, 5 )...
%             minMaxSemiAxisVector( l, 6 )];
%           c = centerEllipse( l, : );
%           minLength = norm( minSemiPoint - c );
%           maxLength = norm( maxSemiPoint - c );
%           
%           % line of major semi axis
%           lineMaxX = [ maxSemiPoint(1), c(1) + c(1)-maxSemiPoint(1) ];
%           lineMaxY = [ maxSemiPoint(2), c(2) + c(2)-maxSemiPoint(2) ];
%           lineMaxZ = [ maxSemiPoint(3), c(3) + c(3)-maxSemiPoint(3) ];
%           
%           % line of minor semi axis
%           % direction have to be here the normal of the x-y plane
%           minorSemiAxes = cross( [ 0 0 1 ], maxSemiPoint-c );
%           minorSemiAxes = normalizeVector3d( minorSemiAxes );
%           lineMinX = [ c(1) - minLength*minorSemiAxes(1), c(1) + minLength*minorSemiAxes(1) ];
%           lineMinY = [ c(2) - minLength*minorSemiAxes(2), c(2) + minLength*minorSemiAxes(2) ];
%           %lineMinZ = [ c(3) - minLength*minorSemiAxes(3), c(3) + minLength*minorSemiAxes(3) ];
%           lineMinZ = [ 0., 0. ];
%           
%           theta = acos( dot(maxSemiPoint - c, [ 1 0 0 ])/norm(maxSemiPoint - c) );
%           theta = theta*180/pi;
%           
%           % check y values because this is required for
%           % a correct rotation of the ellipses around the z axis
%           if lineMaxY(1) <= lineMaxY(2)
%             theta = -theta;
%           end
%           
%           % draw ellipse
%           hold on;
%           [ ELLIP(l) ELLIPPATCH(l) ] = drawEllipse3d( c(1), c(2), c(3)+l*zOffset+0.2, maxLength, minLength, 0, theta );
%           set( ELLIP(l), 'color', [ 0 0 0 ], 'LineWidth', lineWidth );
%           set( ELLIPPATCH(l), 'FaceColor', [ 1 1 1 ], 'FaceLighting', 'none' );
%           
%           if strcmp( lineStr( 1, lineRenderType ), 'renderLargest2DElongation' )
%             hold on;
%             MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', 'k', 'LineWidth', lineWidth );
%           elseif strcmp( lineStr( 1, lineRenderType ), 'renderAll2DElongation' )
%             hold on;
%             MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', 'k', 'LineWidth', lineWidth );
%             hold on;
%             MIN(l) = line( lineMinX, lineMinY, lineMinZ+l*zOffset+0.2, 'Color', 'k', 'LineWidth', lineWidth );
%           elseif strcmp( lineStr( 1, lineRenderType ), 'renderLargest3DElongation' )
%             lineMaxX = [ linePos( l, 13 ), linePos( l, 16 ) ];
%             lineMaxY = [ linePos( l, 14 ), linePos( l, 17 ) ];
%             lineMaxZ = [ linePos( l, 15 ), linePos( l, 18 ) ];
%             hold on;
%             MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', 'r', 'LineWidth', lineWidth );
%           elseif strcmp( lineStr( 1, lineRenderType ), 'renderAll3DElongation' )
%             lineMinX = [ linePos( l, 1 ), linePos( l, 4 ) ];
%             lineMinY = [ linePos( l, 2 ), linePos( l, 5 ) ];
%             lineMinZ = [ linePos( l, 3 ), linePos( l, 6 ) ];
%             lineMidX = [ linePos( l, 7 ), linePos( l, 10 ) ];
%             lineMidY = [ linePos( l, 8 ), linePos( l, 11 ) ];
%             lineMidZ = [ linePos( l, 9 ), linePos( l, 12 ) ];
%             lineMaxX = [ linePos( l, 13 ), linePos( l, 16 ) ];
%             lineMaxY = [ linePos( l, 14 ), linePos( l, 17 ) ];
%             lineMaxZ = [ linePos( l, 15 ), linePos( l, 18 ) ];
%             hold on;
%             MAX(l) = line( lineMaxX, lineMaxY, lineMaxZ+l*zOffset+0.2, 'Color', 'r', 'LineWidth', lineWidth );
%             hold on;
%             MID(l) = line( lineMidX, lineMidY, lineMidZ+l*zOffset+0.2, 'Color', 'k', 'LineWidth', lineWidth );
%             hold on;
%             MIN(l) = line( lineMinX, lineMinY, lineMinZ+l*zOffset+0.2, 'Color', 'b', 'LineWidth', lineWidth );
%           end
%         end
%       elseif strcmp( visualizationType( 1, visType ), 'Contour' )
%         if size( centerEllipse, 1 ) > 2
%           if triangulationType == 1
%             K = convhull( centerEllipse( :, 1 ), centerEllipse( :, 2 ) );
%             CONTOUR(dataIndex) = line( centerEllipse( K, 1 ), centerEllipse( K, 2 ), centerEllipse( K, 3 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
%           else
%             [VolC,ShapeC] = alphavol( [ centerEllipse( :, 1 ), centerEllipse( :, 2 ) ], sqrt( alphaRadiiVector( curT, 1 )) );
%             K = ShapeC.bnd(:,1);
%             dimK = size( K, 1 );
%             if dimK > 1
%               K(dimK+1,:) = K(1,:);
%               CONTOUR(dataIndex) = line( centerEllipse( K, 1 ), centerEllipse( K, 2 ), centerEllipse( K, 3 ), 'Color', colors( dataIndex, : ), 'LineWidth', lineWidth );
%             end
%           end
%         end
%       end
    else
      % if too few cells are available for PCA then just continue
      curI = curI + deltaI;
      continue;
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
  %leg = legend( '120830', '121204', '121211', '130508', '130607', '131203' );
  leg = legend( '120830', '121204', '121211', '130508', '130607' );
  set(leg, 'Location', 'NorthWestOutside');
  linecolors = { [ 1 0 1 ], [ 0 1 1 ], [ 1 0 0 ], [ 0 1 0 ], [ 0 0 1 ] };
  legendlinestyles( leg, {}, {}, linecolors );
  
  grid off;
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  title( strcat( 'Normalized Step ', num2str(curI+deltaI), '-Cells', num2str(numNormCellsN) ) );
  
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
    filePath = strcat( imageDir, digit, num2str(curI+deltaI), '-Cells', num2str(numNormCellsN), '.png' );
    
    saveas( gcf, char(filePath) );
  end
  % at last increment the current index loop parameter
  curI = curI + deltaI;
end
