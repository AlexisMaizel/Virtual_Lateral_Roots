function [NS, cellCenters, maxIntensities, thresholdStack] = determineLocalMaxima(...
  imageStack, cellRadius, randomIntensities, membraneFileName )

height = size( imageStack, 1 );
width = size( imageStack, 2 );
slices = size( imageStack, 3 );
connectivity = 26;

% read membrane channel
memImageStack = readTIFstack( char(membraneFileName) );

%thresholdType = 0 is really slow (dilation part).....
thresholdType = 1;

% perform simple global thresholding as a preprocessing step
% for each pixel object set an 0 else background is 1
% apply gauss filtering
imageStack = imgaussfilt3(imageStack, 1.);
maxIntens = max( imageStack(:) );

if thresholdType == 0
  threshold = 1000;%double(maxIntens)/2.;
  disp( strcat( 'Using threshold = ', num2str(threshold) ) );
  thresStack = zeros( height, width, slices );
  [ indFind ] = find( imageStack(:) < threshold );
  thresStack( ind2sub( size(imageStack), indFind ) ) = 1;
  %restoredefaultpath
  %imshow(thresStack)
  %setWorkingPathProperties()
  
  % compute the Euclidean distance between each zero and its nearest
  % non-zero element
  D = bwdist( thresStack );
  
  % find to local maxima in the distance matrix
  length = 1+cellRadius*2;
  msk = true( length, length, length );
  mid = cellRadius+1;
  msk(mid, mid, mid) = false;
  
  disp( 'Applying dilation' );
  
  D_dil = imdilate(D, msk);
  DCC = zeros( size(D,1), size(D,2), size(D,3) );
  [ indDFind ] = find( D(:) > 2 & D(:) >= D_dil(:) );
  DCC( ind2sub( size(D), indDFind ) ) = 1;
else
  disp( 'Determining regional maxima' );
  threshold = 750;%double(maxIntens)/2.;
  disp( strcat( 'Using threshold = ', num2str(threshold) ) );
  [ indFind ] = find( imageStack(:) < threshold );
  imageStack( ind2sub( size(imageStack), indFind ) ) = 0;
  thresholdStack = imageStack;
  DCC = imregionalmax( thresholdStack );
end

disp( 'Determining connected components' );
NCC = bwconncomp( DCC, connectivity );
NS = regionprops( NCC, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );
numCCs = size(NS, 1);
disp( strcat( 'Found cc =', num2str(numCCs) ) );

radiusSearch = 5;
minCellWallVoxels = 7 + 26*(radiusSearch-1);
minCellWallVoxels = 1;%minCellWallVoxels/2
samplingStep = 0.1;

centers = zeros( numCCs, 3 );
for i=1:numCCs
  centers( i, : ) = NS(i, :).Centroid;
end
%diary 'I:\SegmentationResults\Matlab\Segmentation\20160427\Membrane\tmp.txt'
disp( 'Merging closest connected components' );
distMat = triu( squareform( pdist( centers ) ), 0 );
distMat(distMat==0) = Inf;
[ minVal, minIndex ] = min( distMat(:) );
while minVal < cellRadius
  [ minI, minJ ] = ind2sub( size(distMat), minIndex );
  pos1 = centers( minI, : );
  pos2 = centers( minJ, : );
  centers( minI, : ) = (pos1 + pos2)./2;
  centers( minJ, : ) = [];
  distMat = triu( squareform( pdist( centers ) ), 0 );
  % set zero values to inf so that only non-zero entries are considered
  distMat(distMat==0) = Inf;
  [ minVal, minIndex ] = min( distMat(:) );
end
numCCs = size(centers, 1);
disp( strcat( 'Remaining cc =', num2str(numCCs) ) );

disp( 'Applying nearest neighbor search' );
cellCenters = [];
[ ccIndex, dist ] = rangesearch( centers, centers, cellRadius*2.5 );
for i=1:numCCs
  numCellsInRange = size( ccIndex{ i, 1 }, 2 );
  % each cell will have at least one neighbor in
  % range which is the cell itself
  pos1 = centers( ccIndex{ i, 1 }(1), : );
  if numCellsInRange > 1
    cellWall = 0;
    for c=2:numCellsInRange
      pos2 = centers( ccIndex{ i, 1 }(c), : );
      % if at least one cell wall was detected between the center cell
      % pos1 and another cell in range pos2 then the current center cell
      % is added to cellCenters
      %cenPos1 = [ pos1(1) -pos1(2)+height pos1(3) ];
      %cenPos2 = [ pos2(1) -pos2(2)+height pos2(3) ];
      cellWall = bitor( cellWall, determineCellWallIntersection( pos1, pos2,...
        memImageStack, radiusSearch, samplingStep, minCellWallVoxels ) );
    end
    if cellWall == 1
      cellCenters = [ cellCenters; pos1 ];
    end
  else
    % else only the cell itself is its neighbor and we add it as a cell
    cellCenters = [ cellCenters; pos1 ];
  end
end
%diary off

% at last generate intensities for each cell center
maxIntensities = zeros( size( cellCenters, 1 ), 1 );
for i=1:size( cellCenters, 1 )
  if randomIntensities == 1
    a = 800;
    b = 1300;
    r = (b-a).*rand(1,1) + a;
    maxIntensities( i, 1 ) = int64(r);
  else
    maxIntensities( i, 1 ) = max( imageStack( NS(i, :).PixelIdxList ) );
  end
end
