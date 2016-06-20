function [S, cellCenters, maxIntensities, width, height, slices] =...
  findConnectedComponents( fileName, membraneFileName, minVoxelCount,...
  voxelCountPerCell, anisotropyZ )
imageStack = readTIFstack( char(fileName) );

% create 3D array of binary image data
height = size( imageStack, 1 );
width = size( imageStack, 2 );
slices = size( imageStack, 3 );

% generate convex hull/alpha shape of cc and determine new cc size to apply the
% clustering based on that size but to the original cc structure
%useSurfaceGeneration = 0;

thicknessThreshold = 10;

% type of method how to map a connected component (cc) onto a cell
% type == 0 -> nearest cc are merged and these represent a cell
% type == 1 -> consider a user-chosen cell size to check if a cc consists
% of one or more cells and also ignore cc smaller than some threshold
% type == 2 -> consider a distance mask and its number of local maxima to
% check if a cc consists of one or more cells and also ignore cc smaller
% than some threshold
% type == 3 -> same as 1 but also consider membrane channel to check if for
% a specific cc a cell wall exists in between
% type == 4 -> combination of type 0 and 3
type = 3;

randomIntensities = 1;

r = 50;

% determine connected components in image
% can be 6, 18 or 26
connectivity = 26;
CC = bwconncomp( imageStack( :, :, : ), connectivity );

% get the centroids of each cc
S = regionprops( CC, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );

% generate adjacency matrix for each cc
%ccAdjacencyMap = generateAdjacencyMatrixOfCC( CC, S, connectivity );
if type == 0
  % merge nearest cc
  if size(S, 1) > 1
    X = zeros( size(S, 1), 3 );
    for i=1:size(S, 1)
      X( i, : ) = S(i,:).Centroid;
    end
    %Y = pdist(X);
    Z = linkage( X, 'centroid' );
    T = cluster( Z, 'cutoff', r, 'criterion', 'distance' );
    numC = max( T );
    cellCenters = zeros( numC, 3 );
    maxIntensities = zeros( numC, 1 );
    for i=1:size(T, 1)
      if randomIntensities == 1
        a = 800;
        b = 1300;
        r = (b-a).*rand(1,1) + a;
        maxInt = int64(r);
      else
        maxInt = max( imageStack( S(i, :).PixelIdxList ) );
      end
      cellCenters( T(i,1), : ) = cellCenters( T(i,1), : ) + S(i, :).Centroid;
      maxIntensities( T(i,1), 1 ) = maxIntensities( T(i,1), 1 ) + maxInt;
    end
    for c=1:numC
      cellCenters( c, : ) = cellCenters( c, : ) / nnz(T==c);
      maxIntensities( c, 1 ) = maxIntensities( c, 1 ) / nnz(T==c);
    end
  end
  return;
elseif type == 1
  ccLocalMaximaMap = determineNumberOfLocalMaximaOfCC( CC, S, minVoxelCount, connectivity, anisotropyZ );
elseif type == 3
  % read membrane channel
  memImageStack = readTIFstack( char(membraneFileName) );
elseif type == 4
  % at first merge nearest ccs
  if size(S, 1) > 1
    X = zeros( size(S, 1), 3 );
    for i=1:size(S, 1)
      X( i, : ) = S(i,:).Centroid;
    end
    Z = linkage( X, 'centroid' );
    T = cluster( Z, 'cutoff', r, 'criterion', 'distance' );
    numC = max( T );
    centroidPerCC = zeros( numC, 3 );
    areaPerCC = zeros( numC, 1 );
    pixelListPerCC = cell( numC, 1 );
    for i=1:size(T, 1)
      centroidPerCC( T(i,1), : ) = centroidPerCC( T(i,1), : ) + S(i, :).Centroid;
      pixelListPerCC{ T(i,1), 1 } = [ pixelListPerCC{ T(i,1), 1 } ; S(i,:).PixelList ];
      areaPerCC( T(i,1), 1 ) = areaPerCC( T(i,1), 1 ) + S(i, :).Area;
    end
    for c=1:numC
      centroidPerCC( c, : ) = centroidPerCC( c, : ) / nnz(T==c);
    end
  end
  % read membrane channel
  memImageStack = readTIFstack( char(membraneFileName) );
end

%diary('matlabOutput.txt')
cellCenters = [];
maxIntensities = [];
if type ~= 4
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
    
    if type == 1
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
    elseif type == 2
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
    elseif type == 3
      %diary('matlabOutput.txt')
      areaVal
      centroid = S(i, :).Centroid;
      % get the maximum number of assumed cells based on the considered cell
      % size
      if areaVal > k*voxelCountPerCell
        %disp('previous center');
        %cen = S(i,:).Centroid
        while areaVal > (k+1)*voxelCountPerCell
          k = k + 1;
        end
        maxK = k;
        % now apply for each k > 1 (beginning with the largest) a clustering
        % and identification of cluster centroids; if between the connection
        % of these centroids a cell wall was found in the membrane channel
        % then there ARE at least two cells
        %for m=maxK:-1:2
        startM = [];
        for mm=1:maxK
          startM = [ startM; centroid ];
        end
        [ ~, C ] = kmeans( S(i,:).PixelList, maxK, 'Start', startM );
        %end
        %maxK
        radiusSearch = 5;
        minCellWallVoxels = 7 + 26*(radiusSearch-1);
        minCellWallVoxels = minCellWallVoxels/2;
        samplingStep = 0.05;
        reduceCluster = 1;
        while reduceCluster == 1
          %maxK
          cellWall = 1;
          for ii=1:maxK
            for jj=1:maxK
              if jj>ii
                pos1 = [ C(ii, 1) C(ii, 2) C(ii, 3) ]
                pos2 = [ C(jj, 1) C(jj, 2) C(jj, 3) ]
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
            %disp('new centers');
            %C( c, : )
          end
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
      %diary off
    end
  end
else
  % type == 4
  for i=1:numC
    k = 2;
    
    % max intensities
    a = 800;
    b = 1300;
    r = (b-a).*rand(1,1) + a;
    maxInt = int64(r);
    
    centroid = centroidPerCC( i, : );
    areaVal = areaPerCC( i, : );
    % get the maximum number of assumed cells based on the considered cell
    % size
    if areaVal > k*voxelCountPerCell
      %disp('previous center');
      %cen = S(i,:).Centroid
      while areaVal > (k+1)*voxelCountPerCell
        k = k + 1;
      end
      maxK = k;
      % now apply for each k > 1 (beginning with the largest) a clustering
      % and identification of cluster centroids; if between the connection
      % of these centroids a cell wall was found in the membrane channel
      % then there ARE at least two cells
      %for m=maxK:-1:2
      startM = [];
      for mm=1:maxK
        startM = [ startM; centroid ];
      end
      [ ~, C ] = kmeans( pixelListPerCC{ i, 1 }, maxK, 'Start', startM );
      %end
      %maxK
      reduceCluster = 1;
      while reduceCluster == 1
        cellWall = 1;
        for ii=1:maxK
          for jj=1:maxK
            if jj>ii
              pos1 = [ C(ii, 1) C(ii, 2) C(ii, 3) ];
              pos2 = [ C(jj, 1) C(jj, 2) C(jj, 3) ];
              cellWall = bitand( cellWall, determineCellWallIntersection( pos1, pos2, memImageStack, 0.05, 30 ) );
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
            [ ~, C ] = kmeans( pixelListPerCC{ i, 1 }, maxK, 'Start', startM );
          elseif maxK > 1
            maxK = maxK - 1;
            break;
          else
            break;
          end
        end
      end
      
      if maxK == 1
        cellCenters = [ cellCenters; centroid ];
        maxIntensities = [ maxIntensities ; maxInt ];
      else
        for c=1:maxK
          cellCenters = [ cellCenters; C(c, 1) C(c, 2) C(c, 3) ];
          maxIntensities = [ maxIntensities ; maxInt ];
          %disp('new centers');
          %C( c, : )
        end
      end
      % should be a single cell
      %k
    else
      %disp('current center');
      %cen = S(i,:).Centroid
      cellCenters = [ cellCenters; centroid ];
      maxIntensities = [ maxIntensities ; maxInt ];
      %k=1
    end
  end
end
%diary off
