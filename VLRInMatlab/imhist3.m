function [counts, binLocations, mi, ma] = imhist3( imageStack )
m = size( imageStack );
imageStack = reshape(imageStack, m(1)*m(2)*m(3), 1);
mi = min( imageStack(:) );
ma = max( imageStack(:) );
imageStack = im2uint16( double( double(imageStack-mi)./double(ma-mi) ) );
[counts, binLocations] = imhist( imageStack );