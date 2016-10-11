% This algo applies a membrane segmentation inspired by the idea of Kawase et al (2015):
% A direction-selective local-thresholding method, DSLT, in combination with a dye-based method
% for automated three-dimensional segmentation of cells and airspaces in developing leaves.
% Additionally, a cell shape extractor is included that enables the detection and generation of
% cell shapes based on the segmented membrane information.
setWorkingPathProperties()

chosenData = 6;
dataStr = { '120830_raw' '121204_raw_2014' '121211_raw' '130508_raw' '130607_raw' };
rawDataStr = { '120830' '121204' '121211' '130508' '130607' '20160427' '20160428' '20160426' '20160706' };
startT = 10;
endT = 10;

% store the result as a png file
storePNG = 0;

% after applying DSLT show the extraced 3D cell shapes: WARNING: Needs a
% lot of time due to lots of triangles for the meshes.
% if this is set to one the other 2D plots are not shown.
show3DCellShapes = 0;

% show the MIPs of the DSLT result
showMIPs = 1;

% angle ranges for DSLT step
thetaStart = -pi/2;
thetaStep = pi/4;
thetaEnd = pi/2-thetaStep;
psiStart = -pi;
psiStep = pi/4;
psiEnd = -psiStep;

% DSLT parameters
lineHalfLength = 11;
lineSampling = 1;
lineSteps = lineHalfLength*2 + 1;
% Two parameters given in the paper which I actually set to default values
% and therefore not really affecting the outcome
% alpha is in [0, 1]
alpha = 1;
% const was set in paper to 0.004
cconst = 0.004;

% morphological parameters
% maskType: 0: cube, 1: sphere, 2: cross
maskType = 2;
length = 5;

% the parameters lineHalfLength and length fundamentally affect the outcome
% lineHalfLength: radius/length of line segment kernel used for DSLT in
% membrane segmentation
% length: radius of morphological filter (i used a sphere filter) to
% improve segmentation result

% number of plots used: should not be changed
numPlots = 4;
plotIndex = 1;

if show3DCellShapes == 1
  numPlots = 1;
end

%profile on
% start measuring elapsed time
tic
for t=startT:endT
  plotIndex = 1;
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
    inputPath = strcat( 'I:\NewDatasets\ilastikWorkshopData\20160427\membrane\slice_small_cropped_membrane_T', digit, num2str(t), '_Angle1.tif' );
    %inputPath = strcat( 'I:\NewDatasets\ilastikWorkshopData\20160427\membrane\small_cropped_membrane_T', digit, num2str(t), '_Angle1.tif' );
    %inputPath = strcat( 'I:\NewDatasets\ilastikWorkshopData\20160427\membrane\division_cropped_membrane_T', digit, num2str(t), '_Angle1.tif' );
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
    winWidth = 1000;
    winHeight = 1000;
    if showMIPs == 1
      winHeight = winHeight/2;
    end
    title = strcat( 'MembSeg_T', num2str(t) );
    f = figure( 'Name', char(title), 'Position', [ 100 50 winWidth winHeight ] );
    % activate orbit rotation by default
    cameratoolbar( 'SetMode', 'orbit' );
    % activate none coord system by default for not resetting the camera up
    % vector when rotating
    cameratoolbar( 'SetCoordSys', 'none' );
    % show camera toolbar by default
    cameratoolbar( 'Show' );
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
    clear title
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
    clear title
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
    tit = strcat( 'MorphBinary\_', 'T', {' '}, num2str(t) );
    clear title
    title( char(tit) );
    plotIndex = plotIndex + 1;
  end
  
  invImageStack = imcomplement( binImageStack );
  
  if showMIPs == 1
    subplot( 1, numPlots, plotIndex );
    restoredefaultpath
    h4 = showMIP( invImageStack );
    setWorkingPathProperties()
    tit = strcat( 'MorphInverse\_', 'T', {' '}, num2str(t) );
    clear title
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
    %invImageStack = imcomplement( binImageStack );
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
  
  if storePNG == 1
    filePath = strcat( imageDir, rawDataStr( 1, chosenData ), '_DSLT_T', digit,...
      num2str(t), '_L', num2str(lineHalfLength), '_R', num2str(length), '.png' );
    export_fig( gcf, char(filePath), '-m2', '-png' );
  end
end
% print elapsed time
toc