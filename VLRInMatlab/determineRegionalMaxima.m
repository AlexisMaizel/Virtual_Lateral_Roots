function [S, cellCenters, maxIntensities] = determineRegionalMaxima(...
  fileName, cellRadius, randomIntensities, membraneFileName )

imageStack = readTIFstack( char(fileName) );
% create 3D array of binary image data
height = size( imageStack, 1 );
width = size( imageStack, 2 );
slices = size( imageStack, 3 );

% apply gauss filtering
imageStack = imgaussfilt3(imageStack, 1.);

connectivity = 26;
%regMax = imregionalmax( imageStack( :, :, : ), connectivity );
[ Maxima, MaxPos, Minima, MinPos ] = MinimaMaxima3D( double(imageStack( :, :, : )), 1, 0 );
imageStack = zeros( height, width, slices, 'uint16' );
size(Maxima, 1)
for m=1:size(Maxima, 1)
  imageStack( MaxPos( m, 1 ), MaxPos( m, 2 ), MaxPos( m, 3 ) ) = 1;
end

CC = bwconncomp( imageStack( :, :, : ), connectivity );
S = regionprops( CC, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );

% read membrane channel
memImageStack = readTIFstack( char(membraneFileName) );

radiusSearch = 5;
minCellWallVoxels = 7 + 26*(radiusSearch-1);
minCellWallVoxels = minCellWallVoxels/2;
samplingStep = 0.05;

numCells = size(S, 1)
cellCenters = [];
centers = zeros( numCells, 3 );
for i=1:numCells
  centers( i, : ) = S(i, :).Centroid;
end

[ idx, dist ] = rangesearch( centers, centers, cellRadius );

for i=1:numCells
  numCellsInRange = size( idx( i, 1 ), 2 );
  % each cell will has at least one neighbor in
  % range which is the cell itself
  if numCellsInRange > 1
    needMerge = 1;
    numUpdatedCells = numCellsInRange;
    cen = zeros( numUpdatedCells, 3 );
    for c=1:numUpdatedCells
      cen( c, : ) = centers( idx( i, c ), : );
    end
    while needMerge == 1
      cellWall = zeros( numUpdatedCells, numUpdatedCells );
      separatedCells = 1;
      for j=1:numUpdatedCells
        pos1 = cen( j, : );
        for k=1:numUpdatedCells
          if k>j
            pos2 = cen( k, : );
            cellWall( j, k ) = determineCellWallIntersection( pos1, pos2,...
              memImageStack, radiusSearch, samplingStep, minCellWallVoxels );
            separatedCells = bitand( separatedCells, cellWall( j, k ) );
          end
        end
      end
      if separatedCells == 1
        for c=1:numUpdatedCells
          cellCenters = [ cellCenters; cen( c, : ) ];
        end
        needMerge = 0; % or break;
      else
        if numUpdatedCells > 1
          found = 0;
          for j=1:numUpdatedCells
            for k=1:numUpdatedCells
              if k>j
                % find first event of having no cell wall and merge cell
                % centers
                if cellWall( j, k ) == 0
                  cen( j, : ) = (cen( j, : ) + cen( k, : ))./2;
                  cen( k, : ) = [];
                  found = 1;
                  break;
                end
              end
            end
            if found == 1
              break;
            end
          end
          numUpdatedCells = numUpdatedCells - 1;
        else
          cellCenters = [ cellCenters; cen( 1, : ) ];
          needMerge = 0; % or break;
        end
      end
    end
  else
    % regional maxima has no other maxima within cellRadius
    cellCenters = [ cellCenters;  S(i, :).Centroid ];
  end
end

% at last generate intensities for each cell center
maxIntensities = zeros( size( cellCenters, 1 ), 1 );
for i=1:size( cellCenters, 1 )
  if randomIntensities == 1
    a = 800;
    b = 1300;
    r = (b-a).*rand(1,1) + a;
    maxIntensities( i, 1 ) = int64(r);
  else
    maxIntensities( i, 1 ) = max( imageStack( S(i, :).PixelIdxList ) );
  end
end
