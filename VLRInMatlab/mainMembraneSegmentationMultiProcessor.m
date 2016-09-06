setWorkingPathProperties()

chosenData = 6;
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };
rawDataStr = { '120830' '121204' '121211' '130508' '130607' '20160427' '20160428' '20160426' '20160706' };
startT = 50;
endT = 50;

%showDSLT = 0;
storePNG = 0;
show3DCellShapes = 1;
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
% maskType: 0: cube, 1: sphere, 2: cross
maskType = 2;
length = 5;

if chosenData < 6
  anisotropyZ = 2.;
elseif chosenData == 6
  anisotropyZ = 1.;
elseif chosenData == 7
  anisotropyZ = 4.;
elseif chosenData == 8
  anisotropyZ = 2;%7.5;
end

numPlots = 4;
plotIndex = 1;

if show3DCellShapes == 1
  numPlots = 1;
end

%profile on
% start measuring elapsed time
tic
for lineHalfLength = 11:2:11
  for length = 9:2:9
    lineSteps = lineHalfLength*2 + 1;
    plotIndex = 1;
    %L = lineHalfLength
    %R = length
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
        %inputPath = strcat( 'I:\NewDatasets\ilastikWorkshopData\20160427\membrane\small_cropped_membrane_T', digit, num2str(t), '_Angle1.tif' );
        inputPath = strcat( 'I:\NewDatasets\ilastikWorkshopData\20160427\membrane\division_cropped_membrane_T', digit, num2str(t), '_Angle1.tif' );
      elseif chosenData == 7
        inputPath = strcat( 'I:\NewDatasets\2016-04-28_17.35.59_JENS\Tiffs\membrane\left\cropped_Ch0_CamL_T00', digit, num2str(t), '.tif' );
      elseif chosenData == 8
        inputPath = strcat( 'I:\NewDatasets\Zeiss\20160426\green\test_cropped_spim_TL', digit, num2str(t), '_Angle1.tif' );
      elseif chosenData == 9
        inputPath = strcat( 'I:\NewDatasets\20160706\cropped_beamExp_Ch2_CamL_T00', digit, num2str(t), '.tif' );
      end
      % path to image output
      imageDir = strcat( 'I:\SegmentationResults\Matlab\Segmentation\', rawDataStr( 1, chosenData ), '\Membrane\' );
      [status,message,messageid] = mkdir( char(imageDir) );
      
      imageStack = readTIFstack( char(inputPath) );
      %imageStack = resampleZDimInImage( imageStack, anisotropyZ );
      height = size( imageStack, 1 );
      width = size( imageStack, 2 );
      slices = size( imageStack, 3 );
      
      if showMIPs == 1 || show3DCellShapes == 1
        winWidth = 500;
        winHeight = 500;
        if showMIPs == 1
          winHeight = winHeight/4;
        end
        f = figure( 'Name', 'MembraneSegmentation', 'Position', [ 100 50 winWidth winHeight ] );
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
        if show3DCellShapes == 1
          camproj( 'perspective' );
        else
          camproj( 'orthographic' );
        end
      end
      
      if showMIPs == 1
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
      
      if showMIPs == 1
        subplot( 1, numPlots, plotIndex );
        restoredefaultpath
        h2 = showMIP( imageStack );
        setWorkingPathProperties()
        tit = strcat( 'Gaussian\_', 'T', {' '}, num2str(t) );
        title( char(tit) );
        plotIndex = plotIndex + 1;
      end
      
      % apply DSLT method to extract cell walls
      [ binImageStack ] = applyDSLT( imageStack, psiStart, psiEnd,...
        psiStep, thetaStart, thetaEnd, thetaStep, lineHalfLength,...
        lineSampling, lineSteps, alpha, cconst);
      
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
      
      % generate inner 3D skeleton of cell walls to get thinner cell walls
      % prior to image closing operation
      binImageStack = imclose( binImageStack, se );
      %binImageStack = Skeleton3D( binImageStack );
      
      if showMIPs == 1
        subplot( 1, numPlots, plotIndex );
        restoredefaultpath
        h3 = showMIP( binImageStack );
        setWorkingPathProperties()
        tit = strcat( 'Binary\_', 'T', {' '}, num2str(t) );
        title( char(tit) );
        plotIndex = plotIndex + 1;
      end
      
      if showMIPs == 1
        subplot( 1, numPlots, plotIndex );
        restoredefaultpath
        h4 = showMIP( invImageStack );
        setWorkingPathProperties()
        tit = strcat( 'MorphInverse\_', 'T', {' '}, num2str(t) );
        title( char(tit) );
        plotIndex = plotIndex + 1;
      end
      
      % generate and show extracted 3D cell shapes
      if show3DCellShapes == 1
        subplot( 1, numPlots, plotIndex );
        plotIndex = plotIndex + 1;
        
        % if it is the first time step then initially generate a cell array
        % containing the cells ids for each cell track
        first = 0;
        if t == startT
          first = 1;
          cellTracks = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
        end
        
        %invImageStack = imopen( binImageStack, se );
        %invImageStack = imdilate( binImageStack, se );
        invImageStack = imcomplement( binImageStack );
        %msk = ones( 1, 1, 1 );
        %se2 = strel( 'arbitrary', msk );
        %invImageStack = imdilate( invImageStack, se2 );
        txtPath = strcat( imageDir, 'CellShapesLOD_', rawDataStr( 1, chosenData ), '_T', digit, num2str(t), '.txt' );
        [ h5, cellTracks ] = generate3DCellShapes( invImageStack, txtPath, slices, first, t, cellTracks );
      end
      
      outputPath = strcat( imageDir, rawDataStr( 1, chosenData ), '_Membrane_T', digit, num2str(t), '.tif' );
      binImageStack = uint16(binImageStack).*65535;
      options.message = true;
      writeTIFstack( binImageStack, char(outputPath), 2^31, options );
      %outputinvPath = strcat( imageDir, rawDataStr( 1, chosenData ), '_inv_Membrane_T', digit, num2str(t), '.tif' );
      %invImageStack = uint16(invImageStack).*65535;
      %writeTIFstack( invImageStack, char(outputinvPath), 2^31 );
      
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