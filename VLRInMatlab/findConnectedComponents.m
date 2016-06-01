function [S, cellCenters, maxIntensities, width, height, slices] =...
  findConnectedComponents( fileName, minVoxelCount, voxelCountPerCell, anisotropyZ )
imageStack = readTIFstack( char(fileName) );

% create 3D array of binary image data
height = size( imageStack, 1 );
width = size( imageStack, 2 );
slices = size( imageStack, 3 );

% generate convex hull/alpha shape of cc and determine new cc size to apply the
% clustering based on that size but to the original cc structure
%useSurfaceGeneration = 0;

thicknessThreshold = 10;

considerCellSize = 1;

randomIntensities = 1;

r = 30;

% determine connected components in image
% can be 6, 18 or 26
connectivity = 26;
CC = bwconncomp( imageStack( :, :, : ), connectivity );

% get the centroids of each cc
S = regionprops( CC, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );

% generate adjacency matrix for each cc
%ccAdjacencyMap = generateAdjacencyMatrixOfCC( CC, S, connectivity );
if considerCellSize == 0
  ccLocalMaximaMap = determineNumberOfLocalMaximaOfCC( CC, S, minVoxelCount, connectivity, anisotropyZ );
end

%diary('matlabOutput.txt')
cellCenters = [];
maxIntensities = [];
%subplot( 2, 2, 2 );
%hold on;
for i=1:size(S, 1)
  k = 2;
  areaVal = S(i, :).Area;
  
  % ignore cc with pixeal area smaller than some threshold
  if areaVal < minVoxelCount
    %disp('current center');
    %cen = S(i,:).Centroid
    continue;
  end
  
  % determine max thickness values of single cc
%   [ xmax, ymax, zmax ] = determineCCThickness( S(i,:).PixelList, S(i,:).BoundingBox );
%   % ignore ccs with small thickness values (e.g. cell wall)
%   if xmax < thicknessThreshold || ymax < thicknessThreshold || zmax < thicknessThreshold
%     continue;
%   end
  
  %   if useSurfaceGeneration == 1 && areaVal > 10
  %     coords = S(i,:).PixelList;
  %     %cen = S(i,:).Centroid
  %     alphaRadius = 1000;
  %     shp = alphaShape( coords(:,1), -coords(:,2), coords(:,3), alphaRadius );
  %     %[ K, areaVal ] = convhulln( S(i,:).PixelList );
  %     areaVal = volume(shp);
  %     plot(shp)
  %   end
  
  %[center, radii, evecs, pars ] = ellipsoid_fit( S(i,:).PixelList );
  %center
  %radii
  
  if randomIntensities == 1
    a = 800;
    b = 1300;
    r = (b-a).*rand(1,1) + a;
    maxInt = int64(r);
  else
    maxInt = max( imageStack( S(i, :).PixelIdxList ) );
  end
  
  if considerCellSize == 1
    %areaVal
    centroid = S(i, :).Centroid;
    % if the area size of the cc is bigger than an average cell size
    % apply k means clustering to distinguish between the cells and to
    % find the cell's center
    if areaVal > k*voxelCountPerCell
      % first check based on the current area size of the cc how many
      % clusters shoud be generated
      %disp('previous center');
      %cen = S(i,:).Centroid
      startM = [ centroid ; centroid ];
      while areaVal > (k+1)*voxelCountPerCell
        k = k + 1;
        startM = [ startM; centroid ];
      end
      %k
      [ ~, C ] = kmeans( S(i,:).PixelList, k, 'Start', startM );
      
      for c=1:k
        cellCenters = [ cellCenters; C(c, 1) C(c, 2) C(c, 3) ];
        maxIntensities = [ maxIntensities ; maxInt ];
        %disp('new centers');
        %C( c, : )
      end
      % should be a single cell
      %k
    else
      %disp('current center');
      %cen = S(i,:).Centroid
      cellCenters = [ cellCenters; S(i, :).Centroid ];
      maxIntensities = [ maxIntensities ; maxInt ];
      %k=1
    end
  else
    centroid = S(i, :).Centroid;
    %center = ccLocalMaximaMap(i);
%     for c=1:size( center, 1 )
%       cellCenters = [ cellCenters; center(c, 1) center(c, 2) center(c, 3) ];
%       maxIntensities = [ maxIntensities ; maxInt ];
%       %disp('new centers');
%       %C( c, : )
%     end
    %k = size( center, 1 );
    k = ccLocalMaximaMap(i);
    startM = [];
    if k == 0
      continue;
    else
      for c=1:k
        startM = [ startM; S(i, :).Centroid ];
      end
    end
    
    [ ~, C ] = kmeans( S(i,:).PixelList, k, 'Start', startM );
    for c=1:k
      cellCenters = [ cellCenters; C(c, 1) C(c, 2) C(c, 3) ];
      maxIntensities = [ maxIntensities ; maxInt ];
      %disp('new centers');
      %C( c, : )
    end
    
%     if areaVal > k*voxelCountPerCell
%       % CalinskiHarabasz, DaviesBouldin, gap, silhouette
%       eva = evalclusters(S(i,:).PixelList,'kmeans','CalinskiHarabasz','KList',[1:kmax]);
%       k = eva.OptimalK;
%       [ ~, C ] = kmeans( S(i,:).PixelList, k );
%       for c=1:k
%         cellCenters = [ cellCenters; C(c, 1) C(c, 2) C(c, 3) ];
%         maxIntensities = [ maxIntensities ; maxInt ];
%         %disp('new centers');
%         %C( c, : )
%       end
%     else
%       cellCenters = [ cellCenters; S(i, :).Centroid ];
%       maxIntensities = [ maxIntensities ; maxInt ];
%     end
  end
end
%diary off
