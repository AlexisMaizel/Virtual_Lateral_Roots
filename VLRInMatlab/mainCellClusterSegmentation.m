% Method to identify and separate cells that are part of a cell clump
setWorkingPathProperties()

f = figure( 'Name', 'NucleiSegmentation', 'Position', [ 50 50 1200 900 ] );
numPlotsX = 2;
numPlotsY = 6;
numPlots = 12;
connectivity = 26;
smallVoxelCount = 20000;
borderThreshold = 200;
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
  %I_orig = readTIFstack( char(nucleiFileName3D) );
  height = size( I_orig, 1 );
  width = size( I_orig, 2 );
  slices = size( I_orig, 3 );
  I_orig = imread( char( nucleiFileName ) );
  
  imgData{ plotIndex, 1 } = I_orig;
  titles{ plotIndex, 1 } = 'Original';
  plotIndex = plotIndex + 1;
  
  % gauss filtering
  I_orig = imgaussfilt3(I_orig, 1.);
  imgData{ plotIndex, 1 } = I_orig;
  titles{ plotIndex, 1 } = 'Gauss Filtered';
  plotIndex = plotIndex + 1;
  
  % compute the gradient magnitude image
  [ Gmag, Gazi, Gelevation ] = imgradient3( I_orig );
  imgData{ plotIndex, 1 } = Gmag;
  titles{ plotIndex, 1 } = 'Gradient Magnitude';
  plotIndex = plotIndex + 1;
  
  % Regional maxima of opening-closing by reconstruction
  Io = imopen( I_orig, se );
  Ie = imerode( I_orig, se );
  Iobr = imreconstruct( Ie, I_orig );
  Ioc = imclose( Io, se );
  Iobrd = imdilate( Iobr, se );
  Iobrcbr = imreconstruct( imcomplement(Iobrd), imcomplement(Iobr) );
  Iobrcbr = imcomplement( Iobrcbr );
%   mask_em = imregionalmax( Iobrcbr );
  
  % foreground markers by h-maxima transformation
  % the lower h, the more seeds (therefore more per object) are generated;
  % it is important that each object has at least one seed
  h = 2;
  fore_mark = imextendedmax(I_orig, h);
  %fore_mark = imclose( mask_em, se );
  %fore_mark = imfill( mask_em, 'holes' );
  % Remove small objects from binary image smaller than number of pixels
  %fore_mark = bwareaopen( mask_em, smallVoxelCount );
  imgData{ plotIndex, 1 } = fore_mark;
  titles{ plotIndex, 1 } = 'Foreground Seeds';
  plotIndex = plotIndex + 1;
  
  CC = bwconncomp( fore_mark( :, :, : ), connectivity );
  outp = sprintf( 'Foreground markers found: %d', CC.NumObjects );
  disp( outp )
  
  % determine properties of ccs
  %S = regionprops( CC, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );
  
  % mark the background
%   back_mark = imbinarize( Iobrcbr );
%   back_mark = imcomplement( back_mark );
%   imgData{ plotIndex, 1 } = back_mark;
%   titles{ plotIndex, 1 } = 'Background Seeds';
%   plotIndex = plotIndex + 1;
  
  % compute the background seeds by regional minima of grad magnitude image
  back_mark = imextendedmin( Gmag, 80 );
  % Remove small objects within cells that also feature a regional minima
  back_mark = bwareaopen( back_mark, smallVoxelCount );
  back_mark = imcomplement( back_mark );
  back_mark = imclose( back_mark, se );
  back_mark = imfill( back_mark, 'holes' );
  back_mark = bwareaopen( back_mark, 100 );
  back_mark = imcomplement( back_mark );
  imgData{ plotIndex, 1 } = back_mark;
  titles{ plotIndex, 1 } = 'Background Seeds';
  plotIndex = plotIndex + 1;
  
%   D = bwdist(back_mark);
%   DL = watershed( D );
%   backW = DL == 0;
%   imgData{ plotIndex, 1 } = backW;
%   titles{ plotIndex, 1 } = 'Watershed ridge lines';
%   plotIndex = plotIndex + 1;
  
  gradmag2 = imimposemin( Gmag, back_mark | fore_mark );
  imgData{ plotIndex, 1 } = gradmag2;
  titles{ plotIndex, 1 } = 'Imposed Image';
  plotIndex = plotIndex + 1;
  
  L = watershed( gradmag2 );
  imgData{ plotIndex, 1 } = label2rgb(L);
  titles{ plotIndex, 1 } = 'Watershed';
  plotIndex = plotIndex + 1;
  
  CCL = bwconncomp( L( :, :, : ), connectivity );
  outp = sprintf( 'Watershed regions found: %d', CCL.NumObjects );
  disp( outp )
  
  WS = L == 0;
  %WS = imdilate( WS, ones(3, 3, 3) );
  imgData{ plotIndex, 1 } = WS;
  titles{ plotIndex, 1 } = 'PreFinal';
  plotIndex = plotIndex + 1;
  
  % merge watershed regions that share a "weak" border
  [ MergedBorders, borderCenter ] = mergeWeakBorderRegions( Gmag, WS, L, borderThreshold );
  borderCenter
  
  imgData{ plotIndex, 1 } = MergedBorders;
  titles{ plotIndex, 1 } = 'Merged Borders';
  plotIndex = plotIndex + 1;
  
  borderPos = MergedBorders;
  borderPos( MergedBorders ) = 0;
  for p=1:size( borderCenter, 1 )
    borderPos( borderCenter( p, 1 ), borderCenter( p, 2 ), borderCenter( p, 3 ) ) = 1;
  end
  imgData{ plotIndex, 1 } = borderPos;
  titles{ plotIndex, 1 } = 'Border Positions';
  plotIndex = plotIndex + 1;
  
  BinIm = L;
  BinIm(L < 2) = 0;
  BinIm(L > 1) = 1;
  D = bwdist( ~BinIm );
  imgData{ plotIndex, 1 } = D;
  titles{ plotIndex, 1 } = 'Distance Transform';
  plotIndex = plotIndex + 1;
  
  % Complement the distance transform, and force pixels that don't belong
  % to the objects to be at -Inf.
  D = -D;
  D( ~BinIm ) = -Inf;
  DL = watershed( D );
  imgData{ plotIndex, 1 } = label2rgb(DL);
  titles{ plotIndex, 1 } = 'Watershed Distance';
  plotIndex = plotIndex + 1;
  
  for p=1:numPlots
    subplot( numPlotsX, numPlotsY, p )
    %showMIP( imgData{ p, 1 } );
    restoredefaultpath
    imshow( imgData{ p, 1 }, [] )
    title( titles{ p, 1 } )
    setWorkingPathProperties()
  end
end