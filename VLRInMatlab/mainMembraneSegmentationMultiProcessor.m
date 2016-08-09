setWorkingPathProperties()

chosenData = 6;
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };
rawDataStr = { '120830' '121204' '121211' '130508' '130607' '20160427' '20160428' '20160426' '20160706' };
startT = 10;
endT = 10;

%showDSLT = 0;
storePNG = 0;
extractCellShapes = 1;
showMIPs = 0;

thetaStart = -pi/2;
thetaStep = pi/4;
thetaEnd = pi/2-thetaStep;

psiStart = -pi;
psiStep = pi/4;
psiEnd = -psiStep;

% DSLT parameters
lineHalfLength = 5;
lineSampling = 1;
lineSteps = lineHalfLength*2 + 1;
% alpha is in [0, 1]
alpha = 1;
% const was set in paper to 0.004
cconst = 0.004;

% morphological parameters
cubeMask = 0;
length = 5;

numPlots = 5;
plotIndex = 1;

if showMIPs ~= 1
  numPlots = numPlots - 2;
end

if extractCellShapes == 1
  numPlots = 1;
end

%profile on
% start measuring elapsed time
tic
for lineHalfLength = 11:2:11
  for length = 9:2:9
    lineSteps = lineHalfLength*2 + 1;
    plotIndex = 1;
    L = lineHalfLength
    R = length
    for t=startT:endT
      X = sprintf( 'Time step %d', t );
      disp(X)
      if t < 10
        digit = '00';
      elseif t < 100
        digit = '0';
      else
        digit = '';
      end
      
      % input path
      if chosenData == 6
        %inputPath = strcat( 'I:\NewDatasets\Zeiss\20160427\green\cropped_spim_TL', digit, num2str(t), '_Angle1.tif' );
        %inputPath = strcat( 'I:\NewDatasets\ilastikWorkshopData\20160427\membrane\slice_small_cropped_membrane_T', digit, num2str(t), '_Angle1.tif' );
        inputPath = strcat( 'I:\NewDatasets\ilastikWorkshopData\20160427\membrane\small_cropped_membrane_T', digit, num2str(t), '_Angle1.tif' );
      elseif chosenData == 7
        inputPath = strcat( 'I:\NewDatasets\2016-04-28_17.35.59_JENS\Tiffs\membrane\left\cropped_Ch0_CamL_T00', digit, num2str(t), '.tif' );
      elseif chosenData == 8
        inputPath = strcat( 'I:\NewDatasets\Zeiss\20160426\green\cropped_spim_TL', digit, num2str(t), '_Angle1.tif' );
      elseif chosenData == 9
        inputPath = strcat( 'I:\NewDatasets\20160706\cropped_beamExp_Ch2_CamL_T00', digit, num2str(t), '.tif' );
      end
      % path to image output
      imageDir = strcat( 'I:\SegmentationResults\Matlab\Segmentation\', rawDataStr( 1, chosenData ), '\Membrane\' );
      mkdir( char(imageDir) );
      
      imageStack = readTIFstack( char(inputPath) );
      height = size( imageStack, 1 );
      width = size( imageStack, 2 );
      slices = size( imageStack, 3 );
      
      if slices ~= 1
        numPlots = 1;
      end
      
      if lineHalfLength == 5 && length == 5
        winWidth = 1800;
        f = figure( 'Name', 'MembraneSegmentation', 'Position', [ 100 50 winWidth winWidth/numPlots ] );
      end
      
      xMinMax = [ 0 width ];
      yMinMax = [ 0 height ];
      zMinMax = [ 0 slices ];
      for p=1:numPlots
        subplot( 1, numPlots, p )
        axis( [ xMinMax(1) xMinMax(2) yMinMax(1) yMinMax(2) zMinMax(1) zMinMax(2) 0 1 ] );
        axis on
        daspect( [ 1 1 1 ] );
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        if extractCellShapes == 1
          camproj( 'perspective' );
        else
          camproj( 'orthographic' );
        end
      end
      
      if slices == 1 && showMIPs == 1
        subplot( 1, numPlots, plotIndex );
        restoredefaultpath
        h1 = showMIP( imageStack );
        setWorkingPathProperties()
        tit = strcat( 'Original\_', 'T', {' '}, num2str(t) );
        title( char(tit) );
        plotIndex = plotIndex + 1;
      end
      
      % apply gauss filtering
      imageStack = imgaussfilt3(imageStack, 1.);
      
      if slices == 1 && showMIPs == 1
        subplot( 1, numPlots, plotIndex );
        restoredefaultpath
        h2 = showMIP( imageStack );
        setWorkingPathProperties()
        tit = strcat( 'Gaussian\_', 'T', {' '}, num2str(t) );
        title( char(tit) );
        plotIndex = plotIndex + 1;
      end
      
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
      
      if slices == 1 && extractCellShapes ~= 1
        subplot( 1, numPlots, plotIndex );
        restoredefaultpath
        h3 = showMIP( binImageStack );
        setWorkingPathProperties()
        tit = strcat( 'Binary\_', 'T', {' '}, num2str(t) );
        title( char(tit) );
        plotIndex = plotIndex + 1;
      end
      
      % perform morphological operation at the end
      % for now only cube or sphere mask
      if cubeMask == 1
        msk = ones( length, length, length );
      else
        r = (length-1)/2;
        [ xs, ys, zs ] = ndgrid( -r:r, -r:r, -r:r );
        msk = (xs).^2 + (ys).^2 + (zs).^2 <= r.^2;
      end
      se = strel( 'arbitrary', msk );
      invImageStack = imclose( binImageStack, se );
      %invImageStack = imdilate( binImageStack, se );
      invImageStack = imcomplement( invImageStack );
      binImageStack = uint16(binImageStack).*65535;
      
      if slices == 1 && extractCellShapes ~= 1
        subplot( 1, numPlots, plotIndex );
        restoredefaultpath
        h4 = showMIP( invImageStack );
        setWorkingPathProperties()
        tit = strcat( 'MorphInverse\_', 'T', {' '}, num2str(t) );
        title( char(tit) );
        plotIndex = plotIndex + 1;
      end
      
      if extractCellShapes == 1
        if slices == 1
          subplot( 1, numPlots, plotIndex );
        else
          subplot( 1, numPlots, 1 );
        end
        plotIndex = plotIndex + 1;
        connectivity = 26;
        %invOutputPath = strcat( imageDir, rawDataStr( 1, chosenData ), '_InvMembrane_T', digit, num2str(t), '.tif' );
        %invMImageStack = uint16(invImageStack).*65535;
        %writeTIFstack( invMImageStack, char(invOutputPath), 2^31 );
        CC = bwconncomp( invImageStack, connectivity );
        S = regionprops( CC, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );
        numCCs = size(S, 1);
        cm = colormap( jet(numCCs) );
        remainingCCs = 0;
        hold on
        shading interp
        light
        lighting phong
        for i=1:numCCs
          disp( strcat( num2str(i), '/', num2str(numCCs) ) );
          voxels = S(i, :).PixelList;
          numVoxels = size( voxels, 1 );
          if numVoxels > 50 && numVoxels < 500000
            k = boundary( voxels );
            %k = convhull( voxels );
            if slices == 1
              h5 = plot( voxels(k,1), -voxels(k,2)+height );
            else
              color = cm( i, : );
              h5 = trisurf( k, voxels(:,1), -voxels(:,2)+height, voxels(:,3),...
                'Facecolor', color, 'FaceLighting', 'gouraud', 'LineStyle', 'none', 'FaceAlpha', 1 );
            end
            remainingCCs = remainingCCs + 1;
          end
        end
        remainingCCs
        tit = strcat( 'Shape\_', 'T', {' '}, num2str(t) );
        title( char(tit) );
        hold off
      end
      
      % perform morphological operation at the end
      %se = strel( 'cube', 1 );
      %binImageStack = imclose( binImageStack, se );
      % erosion is not really usable
      %binImageStack = imerode( binImageStack, se );
      %binImageStack = imopen( binImageStack, se );
      % fill holes
      %binImageStack = imfill( binImageStack, 'holes' );
      % remove noise
      %binImageStack = medfilt2( binImageStack );
      
      
      outputPath = strcat( imageDir, rawDataStr( 1, chosenData ), '_Membrane_T', digit, num2str(t), '.tif' );
      writeTIFstack( binImageStack, char(outputPath), 2^31 );
      
      %   if showDSLT == 1
      %     radEllip = 0.05;
      %     totalMinSum = 65535;
      %     totalMaxSum = 0;
      %     for i=1:height
      %       for j=1:width
      %         pos = lineArray{ i, j };
      %         if size(pos, 2) == 5
      %           minT = pos(1,5);
      %           if minT < totalMinSum
      %             totalMinSum = minT;
      %           end
      %           if minT > totalMaxSum
      %             totalMaxSum = minT;
      %           end
      %         end
      %       end
      %     end
      %
      %     hold on
      %     for i=1:height
      %       for j=1:width
      %         pos = lineArray{ i, j };
      %         if size(pos, 2) == 5
      %           length = 0.5+double(pos(1,5) - totalMinSum)/double(totalMaxSum - totalMinSum);
      %
      %           p1 = [ double(pos(1,1)) double(pos(1,2)) ];
      %           p2 = [ double(pos(1,3)) double(pos(1,4)) ];
      %
      %           cen = (p1+p2)/2.;
      %           dir = [ double(p2(1)-p1(1)) double(p2(2)-p1(2)) ];
      %           dir = normalize( dir );
      %
      %           dir = dir.*length/2.;
      %           lineX = [ cen(1)-dir(1), cen(1)+dir(1) ];
      %           lineY = [ cen(2)-dir(2), cen(2)+dir(2) ];
      %
      %           line( lineY, lineX, 'LineWidth', 0.1 );
      %
      %           %line( [ pos(1,2) pos(1,4) ], [ pos(1,1) pos(1,3) ], 'LineWidth', 0.1 );
      %         end
      %       end
      %     end
      %   end
      
      if storePNG == 1
        filePath = strcat( imageDir, rawDataStr( 1, chosenData ), '_DSLT_T', digit,...
          num2str(t), '_L', num2str(lineHalfLength), '_R', num2str(length), '.png' );
        export_fig( gcf, char(filePath), '-m2', '-png' );
      end
    end
  end
end
% print elapsed time
toc

%profile viewer