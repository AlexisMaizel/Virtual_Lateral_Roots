% Modified version of DSLT method published by
% Kawase et al (2015), A direction-selective local-thresholding method,
% DSLT, in combination with a dye-based method for automated
% three-dimensional segmentation of cells and airspaces in developing leaves
function [ binImageStack ] = applyDSLT( imageStack, psiStart, psiEnd, psiStep,...
  thetaStart, thetaEnd, thetaStep, lineHalfLength, lineSampling, lineSteps,...
  alpha, cconst)
height = size( imageStack, 1 );
width = size( imageStack, 2 );
slices = size( imageStack, 3 );
numRots = (((psiEnd-psiStart)/psiStep)+1) * (((thetaEnd-thetaStart)/thetaStep)+1);
rots = zeros( height, width, slices, numRots, 'uint16' );
angles = zeros( 2, numRots );
counter = 1;
for theta=thetaStart:thetaStep:thetaEnd
  for psi=psiStart:psiStep:psiEnd
    angles( 1, counter ) = psi;
    angles( 2, counter ) = theta;
    counter = counter + 1;
  end
end

parfor r=1:numRots
  psi = angles( 1, r );
  theta = angles( 2, r );
  %X = sprintf( 'psi = %f, theta = %f', psi, theta );
  %disp(X)
  % initialize kernel with zeros
  rotKernel = zeros( lineSteps, lineSteps, lineSteps );
  % generate rotated line segment kernel
  for l=-lineHalfLength:lineSampling:lineHalfLength
    x = int16(1 + l * cos(theta) * cos(psi));
    y = int16(1 + l * cos(theta) * sin(psi));
    z = int16(1 + l * sin(theta));
    rotKernel( lineHalfLength+x, lineHalfLength+y, lineHalfLength+z ) = 1;
  end
  % define center of kernel
  rotKernel( lineHalfLength+1, lineHalfLength+1, lineHalfLength+1 ) = 2;
  rotKernel = rotKernel./nnz(rotKernel);
  
  % apply kernel on original image
  rots( :, :, :, r ) = imfilter( imageStack, rotKernel, 'symmetric' );
end

% determine minimum of weighted sums among all line rotations
[ minValues, minIndices ] = min( rots(:, :, :, :), [], 4 );
% determine the theta pitch angle for which the weighted sum is minimal
[ xx, yy, zz, ww ] = ind2sub( size( rots(:, :, :, :) ), minIndices );
thetamins = zeros( height, width, slices );
for i=1:height
  for j=1:width
    for k=1:slices
      thetamins( i, j, k ) = angles( 2, ww( i, j, k ) );
    end
  end
end
% correction term
%cor = cconst * ( 1 - alpha * abs(2*thetamin/pi) );
thetamins = abs(thetamins .* 2 ./pi);
thetamins = thetamins .* alpha;
cor = cconst .* ( -thetamins + 1 );

binImageStack = imageStack > (double(minValues) + cor);