function [ gradImage ] = generate3DGradientMagnitudeImage( imageStack )
imageStack = imgaussfilt3(imageStack, 1.);
% filterType == 0 -> simple 3Dgradient filter
% fitlerType == 1 -> 3D Sobel filtering
filterType = 1;
% create gradient magnitude image
if filterType == 0
  xfilt = zeros( 3, 1, 1 );
  yfilt = zeros( 1, 3, 1 );
  zfilt = zeros( 1, 1, 3 );
  xfilt( :, 1, 1 ) = [ -1 0 1 ];
  yfilt( 1, :, 1 ) = [ -1 0 1 ];
  zfilt( 1, 1, : ) = [ -1 0 1 ];
elseif filterType == 1
  xfilt = zeros( 3, 3, 3 );
  yfilt = zeros( 3, 3, 3 );
  zfilt = zeros( 3, 3, 3 );
  xfilt( :, :, 1 ) = [ -1 0 1 ; -2 0 2 ; -1 0 1 ];
  xfilt( :, :, 2 ) = [ -2 0 2 ; -4 0 4 ; -2 0 2 ];
  xfilt( :, :, 3 ) = [ -1 0 1 ; -2 0 2 ; -1 0 1 ];
  yfilt( :, :, 1 ) = [ -1 -2 -1 ; 0 0 0 ; 1 2 1 ];
  yfilt( :, :, 2 ) = [ -2 -4 -2 ; 0 0 0 ; 2 4 2 ];
  yfilt( :, :, 3 ) = [ -1 -2 -1 ; 0 0 0 ; 1 2 1 ];
  zfilt( :, :, 1 ) = [ 1 2 1 ; 2 4 2 ; 1 2 1 ];
  zfilt( :, :, 2 ) = [ 0 0 0 ; 0 0 0 ; 0 0 0 ];
  zfilt( :, :, 3 ) = [ -1 -2 -1 ; -2 -4 -2 ; -1 -2 -1 ];
end
xI = imfilter( imageStack, xfilt, 'conv' );
yI = imfilter( imageStack, yfilt, 'conv' );
zI = imfilter( imageStack, zfilt, 'conv' );
% simple approximation instead of using sqrt of squared values
gradImage = uint16(sqrt( double(xI.*xI + yI.*yI + zI.*zI) ));