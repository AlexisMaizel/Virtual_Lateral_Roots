setWorkingPathProperties()

chosenData = 6;
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };
rawDataStr = { '120830' '121204' '121211' '130508' '130607' '20160427' '20160428' '20160426' };
startT = 10;
endT = 10;

showDSLT = 0;
storePNG = 0;

thetaStart = -pi/2;
thetaStep = pi/4;
thetaEnd = pi/2-thetaStep;

psiStart = -pi;
psiStep = pi/4;
psiEnd = -psiStep;


lineHalfLength = 15;
lineSampling = 1;
lineSteps = lineHalfLength*2 + 1;

%profile on
% start measuring elapsed time
tic
for t=startT:endT
  if t < 10
    digit = '00';
  elseif t < 100
    digit = '0';
  else
    digit = '';
  end
  
  % input path
  if chosenData == 6
    inputPath = strcat( 'I:\NewDatasets\Zeiss\20160427\green\cropped_spim_TL', digit, num2str(t), '_Angle1.tif' );
  elseif chosenData == 7
    inputPath = strcat( 'I:\NewDatasets\2016-04-28_17.35.59_JENS\Tiffs\membrane\left\cropped_Ch0_CamL_T00', digit, num2str(t), '.tif' );
  elseif chosenData == 8
    inputPath = strcat( 'I:\NewDatasets\Zeiss\20160426\green\cropped_spim_TL', digit, num2str(t), '_Angle1.tif' );
  end
  % path to image output
  imageDir = strcat( 'I:\SegmentationResults\Matlab\Segmentation\', rawDataStr( 1, chosenData ), '\Membrane\' );
  mkdir( char(imageDir) );
  
  imageStack = readTIFstack( char(inputPath) );
  height = size( imageStack, 1 );
  width = size( imageStack, 2 );
  slices = size( imageStack, 3 );
  
  % apply gauss filtering
  imageStack = imgaussfilt3(imageStack, 1.0);
  
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
    rotKernel( lineHalfLength+1, lineHalfLength+1, lineHalfLength+1 ) = 3;
    rotKernel = rotKernel./nnz(rotKernel);
    
    % apply kernel on original image
    rots( :, :, :, r ) = imfilter( imageStack, rotKernel, 'symmetric' );
  end
  
  % determine minimum of weighted sums among all line rotations
  % initialize binary image stack
  minValues = min( rots(:, :, :, :), [], 4 );
  
  binImageStack = uint16(imageStack > minValues);
  binImageStack = binImageStack.*65535;
  
  outputPath = strcat( imageDir, rawDataStr( 1, chosenData ), '_Membrane_T', digit, num2str(t), '.tif' );
  writeTIFstack( binImageStack, char(outputPath), 2^31 );
  
  %restoredefaultpath
  %h = imshow( char( outputPath ) );
  %setWorkingPathProperties()
  
  if showDSLT == 1
    radEllip = 0.05;
    f = figure( 'Name', 'MembraneSegmentation', 'Position', [ 50 50 height height ] );
    totalMinSum = 65535;
    totalMaxSum = 0;
    for i=1:height
      for j=1:width
        pos = lineArray{ i, j };
        if size(pos, 2) == 5
          minT = pos(1,5);
          if minT < totalMinSum
            totalMinSum = minT;
          end
          if minT > totalMaxSum
            totalMaxSum = minT;
          end
        end
      end
    end
    
    hold on
    for i=1:height
      for j=1:width
        pos = lineArray{ i, j };
        if size(pos, 2) == 5
          length = 0.5+double(pos(1,5) - totalMinSum)/double(totalMaxSum - totalMinSum);
          
          p1 = [ double(pos(1,1)) double(pos(1,2)) ];
          p2 = [ double(pos(1,3)) double(pos(1,4)) ];
          
          cen = (p1+p2)/2.;
          dir = [ double(p2(1)-p1(1)) double(p2(2)-p1(2)) ];
          dir = normalize( dir );
          
          dir = dir.*length/2.;
          lineX = [ cen(1)-dir(1), cen(1)+dir(1) ];
          lineY = [ cen(2)-dir(2), cen(2)+dir(2) ];
          
          line( lineY, lineX, 'LineWidth', 0.1 );
          
          %line( [ pos(1,2) pos(1,4) ], [ pos(1,1) pos(1,3) ], 'LineWidth', 0.1 );
        end
      end
    end
  end
  
  if storePNG == 1
    tit = strcat( 'DSLT\_', 'T', {' '}, num2str(t) );
    title( char(tit) );
    filePath = strcat( imageDir, rawDataStr( 1, chosenData ), '_DSLT_blurred_T', digit, num2str(t), '.png' );
    export_fig( gcf, char(filePath), '-m20', '-png' );
  end
end
% print elapsed time
toc

%profile viewer