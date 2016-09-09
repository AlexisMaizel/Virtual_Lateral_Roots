function [ sampImage ] = resampleZDimInImage( imageStack, sz )
osize = size( imageStack, 3 );
nsize = floor( sz*osize );
dimX = size( imageStack, 1 );
dimY = size( imageStack, 2 );
sampImage = zeros( dimX, dimY, nsize );
uu = 1:(osize-1)/(nsize-1):osize;
lower = floor( uu );
upper = ceil( uu );
lambda = upper - uu;
for i=1:dimX
  for j=1:dimY
    oldZVals = reshape( double(imageStack( i, j, : )), [], osize );
    newZVals = lambda .* oldZVals( 1, upper ) + (1-lambda) .* oldZVals( 1, lower );
    sampImage( i, j, : ) = floor(newZVals);
  end
end