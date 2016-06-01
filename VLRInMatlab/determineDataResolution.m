function [width, height, slices] = determineDataResolution( t, inputPath )
if t < 10
  digit = '00';
elseif t < 100
  digit = '0';
else
  digit = '';
end
fileName = strcat( inputPath, digit, num2str(t), '.tif' );
imageStack = readTIFstack( char(fileName) );
height = size( imageStack, 1 );
width = size( imageStack, 2 );
slices = size( imageStack, 3 );