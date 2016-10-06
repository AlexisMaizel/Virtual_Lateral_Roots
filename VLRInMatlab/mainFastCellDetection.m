% Applying nuclei segmentation method of Buggenthin et al (2013)
% to our data sets
setWorkingPathProperties()

%figure;
%plotIndex = 1;

for t=3:3
  if t < 10
    digit = '00';
  elseif t < 100
    digit = '0';
  else
    digit = '';
  end
  
  nucleiFileName = strcat( 'I:\SegmentationResults\MIPsRawData\20160426\MIP_small_cropped_t', digit, num2str(t), '.jpg' );
  %nucleiFileName = strcat( 'I:\NewDatasets\Zeiss\20160426\red\cropped_spim_TL', digit, num2str(t), '_Angle1.tif' );
  %membraneFileName = strcat( 'I:\NewDatasets\Zeiss\20160426\green\cropped_spim_TL006_Angle1.tif' );
  
  % create 3D array of binary image data
%   imageStack = readTIFstack( char(nucleiFileName) );
%   height = size( imageStack, 1 );
%   width = size( imageStack, 2 );
%   slices = size( imageStack, 3 );
  
  tiledim = 50;
  lambda = 5;
  minSizeMSER = 500;
  maxSizeMSER = 4000;
  maxVariation = 100;
  maxEcc = .7;
  minSizeSplit = 500;
  maxSizeSplit = 2500;

  imageStack = imread( char( nucleiFileName ) );
  bw = segmentImage( imcomplement( im2uint8(imageStack) ), 'tiledim', tiledim,...
    'lambda', lambda, 'minSizeMSER', minSizeMSER, ...
    'maxSizeMSER', maxSizeMSER, 'maxVariation', maxVariation,...
    'maxEcc', maxEcc, 'minSizeSplit', minSizeSplit, ...
    'maxSizeSplit', maxSizeSplit, 'visualize', true );
  setWorkingPathProperties
  
%   subplot( 2, 2, plotIndex )
%   h = showMIP( imageStack );
%   plotIndex = plotIndex + 1;
%   subplot( 2, 2, plotIndex )
%   h2 = showMIP( nImageStack );
%   plotIndex = plotIndex + 1; 
end