function [ x, y ] = invPairingFunction( z )
w = floor( (sqrt( 8*z+1 )-1)/2 );
y = z - (w*w + w)/2;
x = w - y;