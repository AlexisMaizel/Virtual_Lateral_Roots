% different thresholding types applied to an imageStack
function [ thresholdStack ] = applyThresholding( imageStack )
% thresholdType == 0: manual global thresholding set by user-defined parameter
% (1 is object, 0 is background)
% thresholdType == 1: manual global thresholding set by user-defined parameter
% (0 is background plus remaining object pixels)
% thresholdType == 2: automatic Richard-Calvard thresholding (iterative isodata method)
% thresholdType == 3: automatic balanced histogram thresholding
% thresholdType == 4: automatic Otsu's method
% thresholdType == 5: H-maxima transform
% thresholdType == 6: extended H-maxima transform
thresholdType = 0;

% first apply gauss filtering for blurring and removing noise
imageStack = imgaussfilt3(imageStack, 1.);
height = size( imageStack, 1 );
width = size( imageStack, 2 );
slices = size( imageStack, 3 );
mi = min( imageStack(:) );
ma = max( imageStack(:) );
threshold = 750;
h = mi;%mean( imageStack(:) );

if thresholdType == 0
  thresholdStack = ones( height, width, slices );
  [ indFind ] = find( imageStack(:) < threshold );
  thresholdStack( ind2sub( size(imageStack), indFind ) ) = 0;
elseif thresholdType == 1
  [ indFind ] = find( imageStack(:) < threshold );
  imageStack( ind2sub( size(imageStack), indFind ) ) = 0;
  thresholdStack = imageStack;
elseif thresholdType == 2
  [threshold, thresholdStack] = isodata( imageStack );
  %[threshold, thresholdStack] = isodata( imageStack, 'log' );
elseif thresholdType == 3
  % TODO
  [counts, binLocations] = imhist3( imageStack );
  thresholdStack = imageStack;
elseif thresholdType == 4
  thresholdStack = ones( height, width, slices );
  threshold = 65535*graythresh( imageStack );
  [ indFind ] = find( imageStack(:) < threshold );
  thresholdStack( ind2sub( size(imageStack), indFind ) ) = 0;
elseif thresholdType == 5  
  thresholdStack = imhmax( imageStack, h );
  disp( strcat( 'Using level = ', num2str(h) ) );
elseif thresholdType == 6
  thresholdStack = imextendedmax( imageStack, h );
  disp( strcat( 'Using level = ', num2str(h) ) );
end
disp( strcat( 'Using threshold = ', num2str(threshold) ) );
