% Method to identify and separate cells that are part of a cell
% clump/cluster
% The idea is inspired by the paper of Wählby et al. (2004)
setWorkingPathProperties()

% simple figure settings and output parameters
f = figure( 'Name', 'NucleiSegmentation', 'Position', [ 50 50 1200 900 ] );
numPlotsX = 2;
numPlotsY = 5;
numPlots = 10;
plotIndex = 1;

% export 2D image plot?
storePNG = 1;
% plot the number of labeled regions in the 2D plot
showRegionText = 0;

% important parameters for algortihm:
% number of neighbors in 3D
connectivity = 26;

% radius of gauss filtering of original image
gauR = 1.;

% h value of h-maxima transformation of original image to determine the
% foreground markers
% the lower h, the more seeds (therefore more per object) are generated;
% it is important that each object has at least one seed such that
% watershed will detect the object
hMax = 2; % hMax = 2 for 2D, hMax = 30 for 3D

% h value of h-minima transformation of gradient magnitude image to
% determine the background markers
hMin = 80; % hMin = 80 for 2D, hMin = 480 for 3D

% objects that consists of less voxel are removed (only used for detecting
% the background markers)
smallVoxelCount = 20000;

% intensity value: inner borders with less mean intensity than this value
% are removed and therefore watershed regions are merged
borderThreshold = 200; % borderThreshold = 200 for 2D, borderThreshold = 4000 for 3D

% NOT USED YET
% intensity value: inner borders with less mean intensity than this value
% are removed and therefore watershed regions are merged for the splitting step
splittingThreshold = 200;

% morphological parameters
% maskType: 0 -> cube, 1 -> sphere, 2 -> cross
maskType = 1;
% size of mask
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

tic
% loop over all time steps
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
  
  % create 3D array of binary image data
  %I_orig = readTIFstack( char(nucleiFileName3D) );
  I_orig = imread( char( nucleiFileName ) );
  slices = size( I_orig, 3 );
  if slices == 1
    Data2D = 1;
  else
    Data2D = 0;
  end
  
  imgData{ plotIndex, 1 } = I_orig;
  titles{ plotIndex, 1 } = 'Original';
  plotIndex = plotIndex + 1;
  
  % gauss filtering
  I_orig = imgaussfilt3(I_orig, gauR);
  imgData{ plotIndex, 1 } = I_orig;
  titles{ plotIndex, 1 } = 'Gauss Filtered';
  plotIndex = plotIndex + 1;
  disp( 'Gauss Filtering applied' );
  
  % compute the gradient magnitude image
  [ Gmag, Gazi, Gelevation ] = imgradient3( I_orig );
  imgData{ plotIndex, 1 } = Gmag;
  titles{ plotIndex, 1 } = 'Gradient Magnitude';
  plotIndex = plotIndex + 1;
  disp( 'Gradient image computed' );
  
  % This is also one way to extract the foreground and background markers
  % but I decided to use another approach due to output quality
  % Regional maxima of opening-closing by reconstruction
%   Io = imopen( I_orig, se );
%   Ie = imerode( I_orig, se );
%   Iobr = imreconstruct( Ie, I_orig );
%   Ioc = imclose( Io, se );
%   Iobrd = imdilate( Iobr, se );
%   Iobrcbr = imreconstruct( imcomplement(Iobrd), imcomplement(Iobr) );
%   Iobrcbr = imcomplement( Iobrcbr );
%   mask_em = imregionalmax( Iobrcbr );
% mark the background
%   back_mark = imbinarize( Iobrcbr );
%   back_mark = imcomplement( back_mark );
%   imgData{ plotIndex, 1 } = back_mark;
%   titles{ plotIndex, 1 } = 'Background Seeds';
%   plotIndex = plotIndex + 1;

  
  % foreground markers by h-maxima transformation
  fore_mark = imextendedmax(I_orig, hMax);
  imgData{ plotIndex, 1 } = fore_mark;
  titles{ plotIndex, 1 } = 'Foreground Seeds';
  plotIndex = plotIndex + 1;
  disp( 'Foreground Seeds determined' );
  
  CC = bwconncomp( fore_mark( :, :, : ), connectivity );
  outp = sprintf( 'Foreground markers found: %d', CC.NumObjects );
  disp( outp )
  
  % compute the background seeds by regional minima of grad magnitude image
  back_mark = imextendedmin( Gmag, hMin );
  % Remove small objects within cells that also feature a regional minima
  back_mark = bwareaopen( back_mark, smallVoxelCount );
  % perform operations to fill gaps and holes on the complement;
  % close operation is a dilation followed by an erosion
  back_mark = imcomplement( back_mark );
  back_mark = imclose( back_mark, se );
  back_mark = imfill( back_mark, 'holes' );
  % remove small artifacts
  back_mark = bwareaopen( back_mark, 100 );
  back_mark = imcomplement( back_mark );
  imgData{ plotIndex, 1 } = back_mark;
  titles{ plotIndex, 1 } = 'Background Seeds';
  plotIndex = plotIndex + 1;
  disp( 'Background Seeds determined' );
  
  % modifies the gradient magnitude image using morphological reconstruction
  % so it only has regional minima wherever back_mark or fore_mark is nonzero
  gradmag = imimposemin( Gmag, back_mark | fore_mark );
  % then apply watershed
  WS = watershed( gradmag );
  CCL = bwconncomp( WS( :, :, : ), connectivity );
  S = regionprops( CCL, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );
  outp = sprintf( 'Watershed regions found: %d', CCL.NumObjects );
  disp( outp )
  
  % if 2D data add the region labels to the plot
  if Data2D == 1
    if showRegionText == 1
      pos = [];
      vals = [];
      for c=1:CCL.NumObjects
        pos = [ pos ; S( c, 1 ).Centroid ];
        vals = [ vals ; WS( S( c, 1 ).PixelIdxList( 1, 1 ) ) ];
      end
      textImg = insertText( label2rgb(WS), pos, vals, 'AnchorPoint','LeftBottom' );
      imgData{ plotIndex, 1 } = textImg;
    else
      imgData{ plotIndex, 1 } = label2rgb(WS);
    end
    titles{ plotIndex, 1 } = 'Watershed';
    plotIndex = plotIndex + 1;
  end
  
  % only extract the borders which are labeled by zero from watershed
  borderImage = WS == 0;
  borderImage = imdilate( borderImage, ones(3, 3, 3) );
  imgData{ plotIndex, 1 } = borderImage;
  titles{ plotIndex, 1 } = 'All Borders';
  plotIndex = plotIndex + 1;
  
  % merge watershed regions that share a "weak" border
  [ MergedBorders, mergedWS, innerBorderImage ] =...
    mergeWeakBorderRegions( Gmag, borderImage, WS, borderThreshold );
  disp( 'Weak borders determined' );
  
  % highlight only the inner borders
  if Data2D == 1
    overlay = imoverlay(borderImage, innerBorderImage, [.2 1 .2]);
    imgData{ plotIndex, 1 } = overlay;
    titles{ plotIndex, 1 } = 'Inner Borders';
    plotIndex = plotIndex + 1;
  end
  
  % show the result of the removed borders therefore merged regions
  imgData{ plotIndex, 1 } = MergedBorders;
  titles{ plotIndex, 1 } = 'Merged Borders';
  plotIndex = plotIndex + 1;
  
  CCL = bwconncomp( mergedWS( :, :, : ), connectivity );
  S = regionprops( CCL, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );
  outp = sprintf( 'Merged watershed regions: %d', CCL.NumObjects );
  disp( outp )
  
  % if 2D data add the region labels to the plot
  if Data2D == 1
    if showRegionText == 1
      pos = [];
      vals = [];
      for c=1:CCL.NumObjects
        pos = [ pos ; S( c, 1 ).Centroid ];
        vals = [ vals ; mergedWS( S( c, 1 ).PixelIdxList( 1, 1 ) ) ];
      end
      textImg = insertText( label2rgb(mergedWS), pos, vals, 'AnchorPoint','LeftBottom' );
      textImg( textImg == 1 ) = 0;
      imgData{ plotIndex, 1 } = textImg;
    else
      imgData{ plotIndex, 1 } = label2rgb(mergedWS);
    end
    titles{ plotIndex, 1 } = 'Merged Watershed';
    plotIndex = plotIndex + 1;
  end
  
  % export labeled watershed result as TIFF
  binWS = mergedWS;
  imageDir = strcat( 'I:\SegmentationResults\CellClump\Final_WS_T', digit, num2str(t), '.tif' );
  outThresStack = uint8(binWS);
  options.message = true;
  writeTIFstack( outThresStack, char(imageDir), 2^31, options );
  
  % convert to binary and also export as TIFF
  binWS( binWS > 1 ) = 0;
  binWS = ~binWS;
  imageDir = strcat( 'I:\SegmentationResults\CellClump\Final_Bin_T', digit, num2str(t), '.tif' );
  outThresStack = uint16(binWS).*65535;
  options.message = true;
  writeTIFstack( outThresStack, char(imageDir), 2^31, options );
  
  % Here are some approaches to split cell clusters based on the distance
  % transform with respect to the foreground and background markers; the
  % idea is quite the same as above: the resulting watershed will have too
  % many regions (so over-segmentation) and based on the intensities of the
  % gradient magnitude image weak borders are removed; however it is kind
  % of similar to the previous output and maybe I am missing something
  % here; therefore it is not used yet
%   DWS = bwdist( binWS );
%   imgData{ plotIndex, 1 } = DWS;
%   titles{ plotIndex, 1 } = 'Distance Transform';
%   plotIndex = plotIndex + 1;

%   DWS = -DWS;
%   DWS( logical(binWS) ) = -Inf;
%   splitImg = imimposemin( uint8(DWS), fore_mark | back_mark );
%   SL = watershed( splitImg );
%   imgData{ plotIndex, 1 } = label2rgb(SL);
%   titles{ plotIndex, 1 } = 'Splitting Watershed';
%   plotIndex = plotIndex + 1;
%   
%   % only extract the borders which are labeled by zero from watershed
%   borderImage = SL == 0;
%   borderImage = imdilate( borderImage, ones(3, 3, 3) );
%   imgData{ plotIndex, 1 } = borderImage;
%   titles{ plotIndex, 1 } = 'All Borders';
%   plotIndex = plotIndex + 1;
%   
%   % merge watershed regions that share a "weak" border
%   [ MergedBorders, mergedWS, innerBorderImage ] =...
%     mergeWeakBorderRegions( Gmag, borderImage, SL, splittingThreshold );
%   
%   CCL = bwconncomp( mergedWS( :, :, : ), connectivity );
%   S = regionprops( CCL, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );
%   outp = sprintf( 'Merged watershed regions: %d', CCL.NumObjects );
%   disp( outp )
%   
%   pos = [];
%   vals = [];
%   for c=1:CCL.NumObjects
%     pos = [ pos ; S( c, 1 ).Centroid ];
%     vals = [ vals ; mergedWS( S( c, 1 ).PixelIdxList( 1, 1 ) ) ];
%   end
%   textImg = insertText( label2rgb(mergedWS), pos, vals, 'AnchorPoint','LeftBottom' );
%   textImg( textImg == 1 ) = 0;
%   imgData{ plotIndex, 1 } = textImg;
%   titles{ plotIndex, 1 } = 'Merged Watershed';
%   plotIndex = plotIndex + 1;
  
  % Complement the distance transform, and force pixels that don't belong
  % to the objects to be at -Inf.
%   D = -D;
%   D( ~splitImg ) = -Inf;
%   DL = watershed( D );
%   imgData{ plotIndex, 1 } = label2rgb(DL);
%   titles{ plotIndex, 1 } = 'Watershed Distance';
%   plotIndex = plotIndex + 1;
  
  if Data2D == 1
    for p=1:numPlots
      subplot( numPlotsX, numPlotsY, p )
      %showMIP( imgData{ p, 1 } );
      restoredefaultpath
      imshow( imgData{ p, 1 }, [] )
      title( titles{ p, 1 } )
      setWorkingPathProperties()
    end
    if storePNG == 1
      filePath = strcat( 'I:\SegmentationResults\CellClump\2DData_T', digit, num2str(t), '.png' );
      export_fig( gcf, char(filePath), '-m4', '-png' );
    end
  end
end
toc