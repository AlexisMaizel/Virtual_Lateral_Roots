geomPath = strcat( pwd, '/geom3d' );
splinesPath = strcat( pwd, '/hobbysplines' );
bezierPath = strcat( pwd, '/BezierPatchSurface' );
exportLibPath = strcat( pwd, '/export_fig' );
exportSVGPath = strcat( pwd, '/plot2svg' );
addpath( geomPath );
addpath( splinesPath );
addpath( bezierPath );
addpath( exportLibPath );
addpath( exportSVGPath );

setenv('LC_ALL','C')

%%%%% setting of properties %%%%%%
% color range for different data sets
% { 'm' 'y' 'r' 'g' 'b' 'c' 'dr' };
colors = [ [ 1 0 1 ]; [ 1 1 0 ]; [ 1 0 0 ]; [ 0 1 0 ]; [ 0 0 1 ]; [ 0 1 1 ]; [ 0.5 0 0 ] ];
switchColor = 0; % switch color between white (0) or black (1)
% camera view which is later set by chaning the camera orbit:
% 1 -> top
% 2 -> side
% 3 -> radial
cView = 2;
% startIndex
startI = 19;
% endIndex
endI = 20;
% step size
deltaI = 1;
% min and max index
minI = 1;
maxI = 20;
renderSingleContours = 0;
renderNuclei = 1;
renderContour = 0;
renderMasterContour = 0;
contourEps = 20.;
% exclude nuclei outliers
excludeOutliers = 1;
% render only master file?
renderMasterFile = 0;
renderLineType = 2;
lineStr = { 'renderAllDeformationsInGrid'...
  'renderOnlyAveragedDeformationsInGrid' };
% line width of contour
lineWidth = 2;
% enable z overlapping
overlapping = 1;
% radius of ellipses
radEllip = 3;
% use triangulation based on delaunay or alpha shape
% 1 -> delaunay
% 2 -> alpha shape
triangulationType = 1;
% vector of data strings
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' '131203_raw' };
pureDataStr = { '120830' '121204' '121211' '130508' '130607' };
% vector of view strings
viewStr = { 'Top' 'Side' 'Radial' };
% data Index:
% 1 -> 120830 -> m
% 2 -> 121204 -> y
% 3 -> 121211 -> r
% 4 -> 130508 -> g
% 5 -> 130607 -> b
% 6 -> 131203
% start id of data
startData = 1;
% end id of data
endData = 1;
% num of data
numData = 5;
% resolution of grid
resGrid = 30;
% register data sets based on base instead of dome tip
registerBase = 0;
% apply the registration to all existing cell ranges and not only the one
% that all data sets share (cells in [18,143])
considerAllCells = 0;
% represent the deformation only based on the longest deformation or all
% principal components of the deformation
onlyLongestDeformation = 1;

% figure properties
f = figure( 'Name', 'Mesh Deformation', 'Position', [0 0 1600 1200] );
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
[ divisionProperties, cellDatas, dimData, maxT, numCellsPerTimeStep,...
  centerPosPerTimeStep, totalMinAxes, totalMaxAxes, cellFileMap ] =...
  prepareData( dataStr, startData, endData, numData, 'Ellipses',...
  renderMasterFile, cView, excludeOutliers, 0 );

% drawing handles
CONTOUR = [];
L = [];
SD = [];
ELLIP = [];
ELLIPPATCH = [];

% path to image output
imageDir = strcat( 'images/Deformation/' );
mkdir( char(imageDir) );

% output format of values
format longG

% set colormap and colorbar depending on the number of cell files
cm = jet( 3 );
colormap( cm );

% gca is the current axes handle
set( gca,'nextplot','replacechildren' );
% gcf is the current figure handle
%lighting phong
axis on
set( gcf, 'Renderer', 'zbuffer' );
lighting gouraud
shading interp
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.5;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.6;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
set( gcf, 'Renderer', 'OpenGL' );
set( gcf,'nextplot','replacechildren' );

% initialize 2D grid
[ rows, columns ] =...
  generate2DGrid(...
  [totalMinAxes(1) totalMinAxes(2)],...
  [totalMaxAxes(1) totalMaxAxes(2)], resGrid );

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

% number of total cells per step
contourCounter = 1;
lineDeformationCounter = 1;
nucleiCounter = 1;
averageDeltaT = 0;

% loop over all registered time steps
for curI=startI:deltaI:endI-1
  % check handles and if they exist hide them for redrawing
  hideHandle( CONTOUR );
  hideHandle( L );
  hideHandle( SD );
  hideHandle( ELLIP );
  hideHandle( ELLIPPATCH );
  
  % initialize tile grid for storing start and end points of longest
  % deformation of a cell between two subsequent time steps which includes
  % the orientation and the magnitude of deformation
  tileGrid = cell( rows*columns, 1 );
  tileGridM = cell( rows*columns, 1 );
  
  curCellsC = zeros( numData, 1 );
  curCellsN = zeros( numData, 1 );
  allCellsC = zeros( numData, 1 );
  allCellsN = zeros( numData, 1 );
  
  % store all cell positions for all data sets such that
  % the average contour is generated for all these cells
  allCells = [];
  allMasterCells = [];
  timeStepDiff = 0;
  
  % loop over all data sets
  for dataIndex=startData:endData
    % init transformation of data sets to ensure a specific view and
    % registration of data sets
    [ u, v, dir, planePos, TF ] =...
      initTransformations( dataStr( 1, dataIndex ), cView, registerBase );
    
    % set viewing and lighting direction
    [az, el] = view( [ dir(1), dir(2), dir(3) ] );
    lightangle( az, el );
    
    % get the alpha shape radii for all time steps
    if triangulationType == 2
      alphaRadiiVector = getAlphaRadius( dataStr( 1, dataIndex ) );
    end
    
    % get the corresponding time steps for the registered steps
    [ curTC(dataIndex), numNormCellsC ] = getCorrespondingTimeStep( curI, minI, maxI,...
      maxT(dataIndex), numCellsPerTimeStep{dataIndex},...
      numCellsPerTimeStep{dataIndex}(maxT(dataIndex),1), considerAllCells );
    
    [ curTN(dataIndex), numNormCellsN ] = getCorrespondingTimeStep( curI+deltaI, minI, maxI,...
      maxT(dataIndex), numCellsPerTimeStep{dataIndex},...
      numCellsPerTimeStep{dataIndex}(maxT(dataIndex),1), considerAllCells );
    
    timeStepDiff = timeStepDiff + curTN(dataIndex) - curTC(dataIndex);
    
    % update cell information
    % number of cells for current and next time step and data set
    numCellsC(dataIndex) = numCellsPerTimeStep{dataIndex}(curTC(dataIndex), 1);
    numCellsN(dataIndex) = numCellsPerTimeStep{dataIndex}(curTN(dataIndex), 1);
    
    % matrix of positions for current time step
    matPos1 = zeros( numCellsC(dataIndex), 3 );
    matPos2 = zeros( numCellsN(dataIndex), 3 );
    cellIds1 = zeros( numCellsC(dataIndex), 1 );
    cellIds2 = zeros( numCellsN(dataIndex), 1 );
    cellPrecursors2 = cell( numCellsN(dataIndex),1 );
    
    nc1 = 1;
    nc2 = 1;
    for j=1:dimData( dataIndex )
      % this is a special case for the data set 130508 for which we
      % ignore the two cells that arise in the master cell file with
      % lineage ID 5
      if excludeOutliers == 1
        if strcmp( dataStr( 1, dataIndex ), '130508_raw' ) &&...
            cellDatas{dataIndex}{j, 6} == 5
          continue;
        end
      end
      if cellDatas{ dataIndex }{ j, 5 } == curTC(dataIndex)
        pos = [ cellDatas{ dataIndex }{ j, 2 }...
          cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
        matPos1(nc1, :) = pos;
        cellIds1(nc1, :) = cellDatas{ dataIndex }{ j, 1 };
        nc1 = nc1 + 1;
      end
      if cellDatas{ dataIndex }{ j, 5 } == curTN(dataIndex)
        pos = [ cellDatas{ dataIndex }{ j, 2 }...
          cellDatas{ dataIndex }{ j, 3 } cellDatas{ dataIndex }{ j, 4 } ];
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
    
    if renderNuclei == 1
      for m=1:size(matPosN{dataIndex}, 1)
        p = matPosN{dataIndex}(m,:) - centerPosPerTimeStep{dataIndex}(curTN(dataIndex),:);
        [ x, y, z ] = ellipsoid( p(1), p(2), p(3), 5, 5, 5, 20 );
        surf(x, y, z,'FaceColor',[ 85/255 189/255 189/255 ], 'FaceAlpha', 1, 'LineStyle', 'none', 'FaceLighting', 'gouraud');
      end
    end
    
    % draw principal components
    coeff = getNormalizedPrincipalComponents( dataStr( 1, dataIndex ), 1 );
    arrowLength = 80;
    colors = [ 1 0 0 ; 0 1 0 ; 0 0 1 ];
    for a=1:3
      if a == 1
        coeff(:,a) = -coeff(:,a);
      elseif a == 2
       coeff(:,a) = coeff(:,a);
      end
      P(a) = quiver3( 0, 0, 0,...
        coeff(1,a), coeff(2,a), coeff(3,a),...
        arrowLength, 'LineWidth', 5,...
        'Color', colors( a, : ) );
      hold on;
    end
  end
  
  % consider the average
  timeStepDiff = timeStepDiff/(endData-startData+1);
  averageDeltaT = averageDeltaT + timeStepDiff;
  
  % render only first and last frame
  if curI == startI || curI == endI
    render = 1;
  else
    render = 0;
  end
  
  hold off;
  set( f,'nextplot','replacechildren' );
  daspect( [ 1 1 1 ] );
  grid off;
  axis off;
end