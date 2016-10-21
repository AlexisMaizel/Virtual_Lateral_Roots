function [ mergeImage, mergedWS, innerBorderImage ] = mergeWeakBorderRegions( gradientImage, borderImage,...
  WSImage, borderThreshold )
borderData = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
borderCoord = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
adjacentRegions = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
indices = find( borderImage ~= 0 );
for i=1:size( indices, 1 )
  [ ii, jj, kk ] = ind2sub( size( borderImage ), indices( i, 1 ) );
  % determine all adjacent labeled regions
  labels = determineAdjacentLabels( WSImage, [ ii, jj, kk ] );
  % get intensity at border position
  intens = gradientImage( ii, jj, kk );
  % add intensity value to each border label if none adjacent label belong
  % to the background (=1); this means that we only consider "inner"
  % borders of regions within possible cell structures
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
        
        % add adjacent regions to each occurring label
        if isKey( adjacentRegions, lab ) == 1
          adjacentRegions( lab ) = [ adjacentRegions( lab ) labels ];
        else
          adjacentRegions( lab ) = labels;
        end
      end
    end
  end
end

ARKeys = keys( adjacentRegions );
ARValues = values( adjacentRegions );
for v=1:size( ARKeys, 2 )
  %label = ARKeys{ 1, v };
  vals = unique( ARValues{ 1, v } );
  adjacentRegions( ARKeys{ 1, v } ) = vals;
end

% compute the mean value for each border
keySet = keys( borderData );
valueSet = values( borderData );
numTotalLabels = size( keySet, 2 );
meanBorders = zeros( numTotalLabels, 1 );
for k=1:size( keySet, 2 )
  sum = 0;
  numIntens = size( valueSet{ 1, k }, 2 );
  for s=1:numIntens
    sum = sum + valueSet{ 1, k }( 1, s );
  end
  sum = sum / numIntens;
  meanBorders( k, 1 ) = sum;
end
%meanBorders

% remove "weak" borders
mergedWS = WSImage;
mergeImage = borderImage;
innerBorderImage = borderImage;
% just set all values to zero and store later only the detected inner
% borders
innerBorderImage( innerBorderImage > 0 ) = 0;
valueCoordSet = values( borderCoord );
for m=1:numTotalLabels
  % current label value
  lab = keySet{ 1, m };
  % get liner indices of border coordinates
  linInd = sub2ind( size( mergeImage ), valueCoordSet{ 1, m }( :, 1 ),...
      valueCoordSet{ 1, m }( :, 2 ), valueCoordSet{ 1, m }( :, 3 ) );
  % case of having a weak border
  if meanBorders( m, 1 ) < borderThreshold
    % merge regions if the border has disappeared
    adjRegionLabels = adjacentRegions( lab );
    for n=1:size( adjRegionLabels, 2 )
      nLabel = adjRegionLabels( 1, n );
      if nLabel ~= 0 && nLabel ~= lab
        % replace labels by merged region label
        ind = find( mergedWS == nLabel );
        mergedWS(ind) = lab;
        % replace previous regions by new merged region in the map
        curRegionValues = values( adjacentRegions );
        for r=1:size( curRegionValues, 2 )
          curRegionValues{ 1, r }( curRegionValues{ 1, r } == nLabel ) = lab;
          adjacentRegions( ARKeys{ 1, r } ) = curRegionValues{ 1, r };
        end
      end
    end
    
    % remove coords from borderImage if the border was removed
    mergeImage( linInd ) = 0;
    % replace borders in merged watershed image by merged labeled region
    mergedWS( linInd ) = lab;
  end
  % store all coords of inner borders
  innerBorderImage( linInd ) = 1;
end
