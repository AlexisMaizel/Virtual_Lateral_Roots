function [ mergeImage, borderCenter ] = mergeWeakBorderRegions( gradientImage, borderImage,...
  WSImage, borderThreshold )
borderData = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
borderCoord = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
indices = find( borderImage ~= 0 );
for i=1:size( indices, 1 )
  [ ii, jj, kk ] = ind2sub( size( borderImage ), indices( i, 1 ) );
  % determine all adjacent labeled regions
  labels = determineAdjacentLabels( WSImage, [ ii, jj, kk ] );
  % get intensity at border position
  intens = gradientImage( ii, jj, kk );
  % add intensity value to each border label if none adjacent label belong
  % to the background (=1); this means that we only consider "inner"
  % borders of regions withtin possible cell structures
  validBorder = all( labels ~= 1 );
  if validBorder == 1
    for l=1:size( labels, 2 )
      lab = labels( 1, l );
      % ignore additional border labels
      if lab ~= 0
        if isKey( borderData, lab ) == 1
          borderData( lab ) = [ borderData( lab ) intens ];
          borderCoord( lab ) = [ borderCoord( lab ) ; [ ii, jj, kk ] ];
        else
          borderData( lab ) = intens;
          borderCoord( lab ) = [ ii, jj, kk ];
        end
      end
    end
  end
end

% compute the mean value for each border
keySet = keys( borderData );
valueSet = values( borderData );
numTotalLabels = size( keySet, 2 );
meanBorders = zeros( numTotalLabels, 1 );
borderCenter = zeros( numTotalLabels, 3 );
for k=1:size( keySet, 2 )
  sum = 0;
  numLabels = size( valueSet{ 1, k }, 2 );
  for s=1:numLabels
    sum = sum + valueSet{ 1, k }( 1, s );
  end
  sum = sum / numLabels;
  meanBorders( k, 1 ) = sum;
end

% remove "weak" borders
mergeImage = borderImage;
valueCoordSet = values( borderCoord );
for m=1:numTotalLabels
  centroid = [ 0, 0, 0 ];
  numCoords = size( valueCoordSet{ 1, m }, 1 );
  if meanBorders( m, 1 ) < borderThreshold
    for c=1:numCoords
      coords = valueCoordSet{ 1, m }( c, : );
      centroid = centroid + coords;
      mergeImage( coords(1,1), coords(1,2), coords(1,3) ) = 0;
    end
    centroid = centroid ./ numCoords;
  else
    for c=1:numCoords
      coords = valueCoordSet{ 1, m }( c, : );
      centroid = centroid + coords;
    end
    centroid = centroid ./ numCoords;
  end
  borderCenter( m, : ) = int64(centroid);
end
