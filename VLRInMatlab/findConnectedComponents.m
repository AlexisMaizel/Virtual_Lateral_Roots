function [S, cellCenters, maxIntensities, imageStack] = findConnectedComponents( fileName, minVoxelCount, voxelCountPerCell )
imageStack = readTIFstack( char(fileName) );

% create 3D array of binary image data
height = size( imageStack, 1);
width = size( imageStack, 2);
slices = size( imageStack, 3);

% determine connected components in image
CC = bwconncomp( imageStack( :, :, : ) );

% get the centroids of each cc
S = regionprops( CC, 'Centroid', 'Area', 'PixelList', 'PixelIdxList' );

% testing
%diary('matlabOutput.txt')
cellCenters = [];
maxIntensities = [];
for i=1:size(S, 1)
  k = 2;
  areaVal = S(i, :).Area;
  maxInt = max( imageStack( S(i, :).PixelIdxList ) );
  % ignore cc with pixeal area smaller than some threshold
  if areaVal < minVoxelCount
    %disp('current center');
    %cen = S(i,:).Centroid
    continue;
  % if the area size of the cc is bigger than an average cell size
  % apply k means clustering to distinguish between the cells and to 
  % find the cell's center
  elseif areaVal > k*voxelCountPerCell
    % first check based on the current area size of the cc how many
    % clusters shoud be generated
    %disp('previous center');
    %cen = S(i,:).Centroid
    startM = [ S(i, :).Centroid ; S(i, :).Centroid ];
    while areaVal > (k+1)*voxelCountPerCell
      k = k + 1;
      startM = [ startM; S(i, :).Centroid ];
    end
    %k
    [ ~, C ] = kmeans( S(i,:).PixelList, k, 'Start', startM );
    for c=1:k
      cellCenters = [ cellCenters; C(c, 1) C(c, 2) C(c, 3) ];
      maxIntensities = [ maxIntensities ; maxInt ];
      %disp('new centers');
      %C( c, : )
    end
  else
    %disp('current center');
    %cen = S(i,:).Centroid
    cellCenters = [ cellCenters; S(i, :).Centroid ];
    maxInt = max( imageStack( S(i, :).PixelIdxList ) );
    maxIntensities = [ maxIntensities ; maxInt ];
  end
end
%diary off
