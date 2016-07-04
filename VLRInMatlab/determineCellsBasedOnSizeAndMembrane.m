function [S, cellCenters, maxIntensities] =...
  determineCellsBasedOnSizeAndMembrane( imageStack, voxelCountPerCell,...
  randomIntensities, minVoxelCount, membraneFileName )
% determine connected components in image
% can be 6, 18 or 26
connectivity = 26;
CC = bwconncomp( imageStack( :, :, : ), connectivity );
% determine properties of ccs
S = regionprops( CC, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );

% read membrane channel
memImageStack = readTIFstack( char(membraneFileName) );

%thicknessThreshold = 10;
cellCenters = [];
maxIntensities = [];
for i=1:size(S, 1)
  k = 2;
  areaVal = S(i, :).Area;
  
  % ignore ccs with pixel areas smaller than some threshold
  if areaVal < minVoxelCount
    continue;
  end
  
  % determine max thickness values of single cc
  %   [ xmax, ymax, zmax ] = determineCCThickness( S(i,:).PixelList, S(i,:).BoundingBox );
  %   % ignore ccs with small thickness values (e.g. cell wall)
  %   if xmax < thicknessThreshold || ymax < thicknessThreshold || zmax < thicknessThreshold
  %     continue;
  %   end
  
  if randomIntensities == 1
    a = 800;
    b = 1300;
    r = (b-a).*rand(1,1) + a;
    maxInt = int64(r);
  else
    maxInt = max( imageStack( S(i, :).PixelIdxList ) );
  end
  
  centroid = S(i, :).Centroid;
  % if the area size of the cc is bigger than an average cell size
  % apply k means clustering to distinguish between the cells and to
  % find the cell's center
  if areaVal > k*voxelCountPerCell
    while areaVal > (k+1)*voxelCountPerCell
      k = k + 1;
    end
    maxK = k;
    % now apply for each k > 1 (beginning with the largest) a clustering
    % and identification of cluster centroids; if between the connection
    % of these centroids a cell wall was found in the membrane channel
    % then there ARE at least two cells
    startM = [];
    for mm=1:maxK
      startM = [ startM; centroid ];
    end
    [ ~, C ] = kmeans( S(i,:).PixelList, maxK, 'Start', startM );
    radiusSearch = 5;
    minCellWallVoxels = 7 + 26*(radiusSearch-1);
    minCellWallVoxels = minCellWallVoxels/2;
    samplingStep = 0.05;
    reduceCluster = 1;
    while reduceCluster == 1
      cellWall = 1;
      for ii=1:maxK
        for jj=1:maxK
          if jj>ii
            pos1 = [ C(ii, 1) C(ii, 2) C(ii, 3) ];
            pos2 = [ C(jj, 1) C(jj, 2) C(jj, 3) ];
            cellWall = bitand( cellWall, determineCellWallIntersection( pos1, pos2,...
              memImageStack, radiusSearch, samplingStep, minCellWallVoxels ) );
          end
        end
      end
      if cellWall == 1
        reduceCluster = 0;
      else
        if maxK > 2
          maxK = maxK - 1;
          startM = [];
          for mm=1:maxK
            startM = [ startM; centroid ];
          end
          [ ~, C ] = kmeans( S(i,:).PixelList, maxK, 'Start', startM );
        elseif maxK > 1
          maxK = maxK - 1;
          break;
        else
          break;
        end
      end
    end
    
    if maxK == 1
      cellCenters = [ cellCenters; S(i, :).Centroid ];
      maxIntensities = [ maxIntensities ; maxInt ];
    else
      for c=1:maxK
        cellCenters = [ cellCenters; C(c, 1) C(c, 2) C(c, 3) ];
        maxIntensities = [ maxIntensities ; maxInt ];
      end
    end
  else
    cellCenters = [ cellCenters; S(i, :).Centroid ];
    maxIntensities = [ maxIntensities ; maxInt ];
  end
end