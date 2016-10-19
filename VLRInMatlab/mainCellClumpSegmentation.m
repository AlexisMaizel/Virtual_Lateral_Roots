% Method to identify and separate cells that are part of a cell clump
setWorkingPathProperties()

f = figure( 'Name', 'NucleiSegmentation', 'Position', [ 50 50 1200 900 ] );
numPlotsX = 2;
numPlotsY = 3;
numPlots = 6;
connectivity = 26;
smallVoxelCount = 100;
plotIndex = 1;

% morphological parameters
% maskType: 0: cube, 1: sphere, 2: cross
maskType = 1;
length = 5;
% perform morphological operation at the end
% for now only cube or sphere mask
if maskType == 0
  msk = ones( length, length, length );
elseif maskType == 1
  r = (length-1)/2;
  [ xs, ys, zs ] = ndgrid( -r:r, -r:r, -r:r );
  msk = (xs).^2 + (ys).^2 + (zs).^2 <= r.^2;
elseif maskType == 2
  msk = zeros( length, length, length );
  bW = 1;
  mid = int16(length/2);
  msk( :, mid-bW:mid+bW, mid-bW:mid+bW ) = 1;
  msk( mid-bW:mid+bW, :, mid-bW:mid+bW ) = 1;
  msk( mid-bW:mid+bW, mid-bW:mid+bW, : ) = 1;
end
se = strel( 'arbitrary', msk );

for t=3:3
  if t < 10
    digit = '00';
  elseif t < 100
    digit = '0';
  else
    digit = '';
  end
  
  imgData = cell( numPlots, 1 );
  titles = cell( numPlots, 1 );
  
  nucleiFileName = strcat( 'I:\SegmentationResults\MIPsRawData\20160426\MIP_small_cropped_t', digit, num2str(t), '.jpg' );
  nucleiFileName3D = strcat( 'I:\NewDatasets\Zeiss\20160426\red\test_cropped_spim_TL', digit, num2str(t), '_Angle1.tif' );
  %membraneFileName = strcat( 'I:\NewDatasets\Zeiss\20160426\green\cropped_spim_TL006_Angle1.tif' );
  
  % create 3D array of binary image data
  I_orig = readTIFstack( char(nucleiFileName3D) );
  height = size( I_orig, 1 );
  width = size( I_orig, 2 );
  slices = size( I_orig, 3 );
  
  %I_orig = imread( char( nucleiFileName ) );
  imgData{ plotIndex, 1 } = I_orig;
  titles{ plotIndex, 1 } = 'Original';
  plotIndex = plotIndex + 1;
  
  I_orig = imcomplement( im2uint8(I_orig) );
  I_orig = im2double( I_orig );
  IBC = zeros( height, width, slices );
  ICE = zeros( height, width, slices );
  ITH = zeros( height, width, slices );
  % compute the background for each slice of the image
  parfor s=1:slices
    ISlice = I_orig( :, :, s );
    bg = bgest( ISlice, 50 );
    
    % correct the image (2D only)
    I = ISlice./bg;
    I(I>1) = 1;
    I(I<0) = 0;
    
    I = im2uint8(I);
    if all( I == I(1) )
      I(:) = 0;
      ISlice = I;
    else
      ISlice = imcomplement( I );
    end
    
    IBC( :, :, s ) = ISlice;
    
    % Perform contrast-limited adaptive histogram equalization (CLAHE) (2D only)
    I_contrast = adapthisteq( ISlice );
    ICE( :, :, s ) = I_contrast;
    
    % Convert image to binary image, based on threshold (2D only)
    I_thres = im2bw( I_contrast, graythresh(I_contrast) );
    ITH( :, :, s ) = I_thres;
    %I_thres = imbinarize( I_contrast, graythresh(I_contrast) );
  end
  
  imgData{ plotIndex, 1 } = IBC;
  titles{ plotIndex, 1 } = 'Background corrected';
  plotIndex = plotIndex + 1;
  imgData{ plotIndex, 1 } = ICE;
  titles{ plotIndex, 1 } = 'Contrast enhanced';
  plotIndex = plotIndex + 1;
  imgData{ plotIndex, 1 } = ITH;
  titles{ plotIndex, 1 } = 'Thresholded';
  plotIndex = plotIndex + 1;
  
  imageDir = strcat( 'I:\SegmentationResults\CellClump\Thres_T', digit, num2str(t), '.tif' );
  outThresStack = uint16(ITH).*65535;
  options.message = true;
  writeTIFstack( outThresStack, char(imageDir), 2^31, options );
  
%   bw2 = imfill( I_thres, 'holes' );
%   bw3 = imopen( bw2, se );
%   % Remove small objects from binary image smaller than number of pixels
%   bw4 = bwareaopen( bw3, smallVoxelCount );
%   % Find perimeter of objects in binary image
%   bw4_perim = bwperim( bw4 );
%   % Fills the input image with a solid color where the input binary
%   % mask is true.
%   overlay1 = imoverlay(I_contrast, bw4_perim, [.3 1 .3]);
%   imgData{ plotIndex, 1 } = overlay1;
%   titles{ plotIndex, 1 } = 'Perimeter of Contrast';
%   plotIndex = plotIndex + 1;
  
  % Regional maxima of the H-maxima transform of the contrast image
%   mask_em = imextendedmax(I_contrast, 80);
%   mask_em = imclose( mask_em, se );
%   mask_em = imfill( mask_em, 'holes' );
%   % Remove small objects from binary image smaller than number of pixels
%   mask_em = bwareaopen( mask_em, smallVoxelCount );
%   imgData{ plotIndex, 1 } = mask_em;
%   titles{ plotIndex, 1 } = 'Regional Maxima';
%   plotIndex = plotIndex + 1;

  % Computes the Euclidean distance transform of the binary image
  ITH = imdilate( ITH, se );
  D = bwdist( ~ITH );
  imgData{ plotIndex, 1 } = D;
  titles{ plotIndex, 1 } = 'Dilation';
  plotIndex = plotIndex + 1;
  
  %overlay2 = imoverlay(I_contrast, bw4_perim | mask_em, [.3 1 .3]);
%   imgData{ plotIndex, 1 } = overlay1;
%   titles{ plotIndex, 1 } = 'Perimeter of Contrast/Max';
%   plotIndex = plotIndex + 1;
  
%   I_contrast_c = imcomplement(I_contrast);
%   % Modifies the image using morphological reconstruction so it only has
%   % regional minima wherever ~bw4 | mask_em is nonzero
%   I_mod = imimposemin( I_contrast_c, ~bw4 | mask_em );
%   L = watershed( I_mod );
%   imgData{ 8, 1 } = label2rgb(L);
%   titles{ 8, 1 } = 'Final';
  
  % Complement the distance transform, and force pixels that don't belong
  % to the objects to be at -Inf.
  D = -D;
  D( ~ITH ) = -Inf;
  L = watershed( D );
  
  NS = regionprops( L, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );
  numCCs = 0;
  % set background (=1) and boundary (=0) to zero
  zer = find( L <= 1 );
  L( zer ) = 0;
  for cc=1:size( NS, 1 )
    if NS( cc, 1 ).Area < smallVoxelCount
      NS( cc, 1 ).Area;
      coord = NS( cc, 1 ).PixelList;
      % set smaller cc to background label (=0)
      L( coord( :, 2 ), coord( :, 1 ), coord( :, 3 ) ) = 0;
    else
      numCCs = numCCs + 1;
    end
  end
  outp = sprintf( 'After split: %d', numCCs );
  disp( outp )
  
  %imgData{ plotIndex, 1 } = label2rgb(L);
  imgData{ plotIndex, 1 } = L;
  titles{ plotIndex, 1 } = 'Final';
  plotIndex = plotIndex + 1;
  
  imageDir = strcat( 'I:\SegmentationResults\CellClump\Final_T', digit, num2str(t), '.tif' );
  outThresStack = uint16(L).*65535;
  options.message = true;
  writeTIFstack( outThresStack, char(imageDir), 2^31, options );
  
%   for p=1:numPlots
%     subplot( numPlotsX, numPlotsY, p )
%     showMIP( imgData{ p, 1 } );
%     %restoredefaultpath
%     %imshow( imgData{ p, 1 } )
%     %title( titles{ p, 1 } )
%     setWorkingPathProperties()
%   end
end