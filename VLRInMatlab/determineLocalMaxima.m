function [S, cellCenters, maxIntensities] = determineLocalMaxima(...
  imageStack, randomIntensities, minVoxelCount, membraneFileName )

ultEro = bwulterode( imageStack( :, :, : ) );
nnz(ultEro)

% apply gauss filtering
imageStack = imgaussfilt3(imageStack, 1.);

% determine connected components in image
% can be 6, 18 or 26
connectivity = 26;
%CC = bwconncomp( imageStack( :, :, : ), connectivity );
% determine properties of ccs
%S = regionprops( CC, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );
%ccLocalMaximaMap = determineNumberOfLocalMaximaOfCC( CC, S, minVoxelCount, connectivity );


% TODO: maybe easier with ultimate erosion
ultEro = bwulterode( imageStack( :, :, : ) );
%restoredefaultpath
%imshow(ultEro)
%setWorkingPathProperties()
nnz(ultEro)
CC = bwconncomp( ultEro, connectivity );
% get the centroids of each cc
S = regionprops( CC, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );
% read membrane channel
memImageStack = readTIFstack( char(membraneFileName) );

cellCenters = [];
maxIntensities = [];
return
for i=1:size(S, 1)
  areaVal = S(i, :).Area;
  
  % ignore cc with pixeal area smaller than some threshold
  if areaVal < minVoxelCount
    %disp('current center');
    %cen = S(i,:).Centroid
    continue;
  end
  
  if randomIntensities == 1
    a = 800;
    b = 1300;
    r = (b-a).*rand(1,1) + a;
    maxInt = int64(r);
  else
    maxInt = max( imageStack( S(i, :).PixelIdxList ) );
  end
  
  centroid = S(i, :).Centroid;
  k = ccLocalMaximaMap(i);
  startM = [];
  if k == 0
    continue;
  else
    for c=1:k
      startM = [ startM; centroid ];
    end
  end
  
  [ ~, C ] = kmeans( S(i,:).PixelList, k, 'Start', startM );
  for c=1:k
    cellCenters = [ cellCenters; C(c, 1) C(c, 2) C(c, 3) ];
    maxIntensities = [ maxIntensities ; maxInt ];
  end
end