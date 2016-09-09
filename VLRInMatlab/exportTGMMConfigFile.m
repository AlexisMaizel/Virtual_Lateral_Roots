function [rawDataPath] = exportTGMMConfigFile( dataStr, backgroundThreshold,...
  persistanceSegmentationTau, minNucleiSize, maxNucleiSize, dataType )
% dataType decides which kind of data is used for TGMM: either the raw data
% (=0) or the preprocessed data (=1)
mkdir( 'C:\Jens\TGMM_Supplementary_Software_1_0\build\AutoTGMMConfig\' );
fileName = strcat( 'C:\Jens\TGMM_Supplementary_Software_1_0\build\AutoTGMMConfig\TGMM_configFile', dataStr, '.txt' );
fileID = fopen( char(fileName), 'w' );
dataPathStr = 'imgFilePattern=I:\';
inputPreprocessedPath = strcat( 'SegmentationResults\Preprocessing\', dataStr, '\preprocessed_t???' );
if strcmp(dataStr, '120830')
  if dataType == 0
    dataPathStr = strcat(dataPathStr, 'FrankfurtLSFMDatasets\20120830-pGATA_H2B_Wave\AllTimePoints\Not_Registered\cropped_pGata23_120830_TL0???_CHN00');
  else
    dataPathStr = strcat(dataPathStr, inputPreprocessedPath);
  end
  rawDataPath = 'I:\SegmentationResults\Preprocessing\120830\changed_t001_M.tif';
  anisotropyZ = 2.;
elseif strcmp(dataStr, '121204')
  if dataType == 0
    dataPathStr = strcat(dataPathStr, 'FrankfurtLSFMDatasets\20121204_pGATA_H2B_Wave\driftcorrected_stacks_cropped\ch01\crop_t???_c1');
  else
    dataPathStr = strcat(dataPathStr, inputPreprocessedPath);
  end
  rawDataPath = 'I:\SegmentationResults\Preprocessing\121204\changed_t001_M.tif';
  anisotropyZ = 2.;
elseif strcmp(dataStr, '121211')
  if dataType == 0
    dataPathStr = strcat(dataPathStr, 'FrankfurtLSFMDatasets\20121211_pGATA_H2B_Wave\3Ddrift_Stacks_cropped\cropped_t???_c1');
  else
    dataPathStr = strcat(dataPathStr, inputPreprocessedPath);
  end
  rawDataPath = 'I:\SegmentationResults\Preprocessing\121211\changed_t001_M.tif';
  anisotropyZ = 2.;
elseif strcmp(dataStr, '130508')
  if dataType == 0
    dataPathStr = strcat(dataPathStr, 'FrankfurtLSFMDatasets\130508\crop_t???_c1');
  else
    dataPathStr = strcat(dataPathStr, inputPreprocessedPath);
  end
  rawDataPath = 'I:\SegmentationResults\Preprocessing\130508\changed_t001_M.tif';
  anisotropyZ = 2.;
elseif strcmp(dataStr, '130607')
  if dataType == 0
    dataPathStr = strcat(dataPathStr, 'FrankfurtLSFMDatasets\130607\c1\crop_t???_c1');
  else
    dataPathStr = strcat(dataPathStr, inputPreprocessedPath);
  end
  rawDataPath = 'I:\SegmentationResults\Preprocessing\130607\changed_t001_M.tif';
  anisotropyZ = 2.;
elseif strcmp(dataStr, '20160426')
  if dataType == 0
    dataPathStr = strcat(dataPathStr, 'NewDatasets\Zeiss\20160426\red\cropped_spim_TL???_Angle1');
  else
    dataPathStr = strcat(dataPathStr, inputPreprocessedPath);
  end
  rawDataPath = 'I:\NewDatasets\Zeiss\20160426\red\cropped_spim_TL003_Angle1.tif';
  anisotropyZ = 7.5;
elseif strcmp(dataStr, '20160427')
  if dataType == 0
    dataPathStr = strcat(dataPathStr, 'NewDatasets\Zeiss\20160427\red\cropped_spim_TL???_Angle1');
  else
    dataPathStr = strcat(dataPathStr, inputPreprocessedPath);
  end
  rawDataPath = 'I:\NewDatasets\Zeiss\20160427\red\cropped_spim_TL010_Angle1.tif';
  anisotropyZ = 1.;
elseif strcmp(dataStr, '20160428')
  if dataType == 0
    dataPathStr = strcat(dataPathStr, 'NewDatasets\2016-04-28_17.35.59_JENS\Tiffs\nuclei\left\cropped_Ch1_CamL_T?????');
  else
    dataPathStr = strcat(dataPathStr, inputPreprocessedPath);
  end
  rawDataPath = 'I:\NewDatasets\2016-04-28_17.35.59_JENS\Tiffs\nuclei\left\cropped_Ch1_CamL_T00010.tif';
  anisotropyZ = 4.;
end
fprintf( fileID, '%s\n', char(dataPathStr) );
output = strcat( 'debugPathPrefix=I:\SegmentationResults\TGMM\', dataStr );
fprintf( fileID, '%s\n', char(output) );
anisoStr = strcat( 'anisotropyZ=', num2str( anisotropyZ ) );
fprintf( fileID, '%s\n', char(anisoStr) );
backStr = strcat( 'backgroundThreshold=', num2str( backgroundThreshold ) );
fprintf( fileID, '%s\n', char(backStr) );
tauStr = strcat( 'persistanceSegmentationTau=', num2str( persistanceSegmentationTau ) );
fprintf( fileID, '%s\n', char(tauStr) );
fprintf( fileID, '%s\n', 'betaPercentageOfN_k=0.01' );
fprintf( fileID, '%s\n', 'nuPercentageOfN_k=1.0' );
fprintf( fileID, '%s\n', 'alphaPercentage=100' );
fprintf( fileID, '%s\n', 'maxIterEM=100' );
fprintf( fileID, '%s\n', 'tolLikelihood=1e-6' );
fprintf( fileID, '%s\n', 'regularizePrecisionMatrixConstants_lambdaMin=0.02' );
fprintf( fileID, '%s\n', 'regularizePrecisionMatrixConstants_lambdaMax=0.222' );
fprintf( fileID, '%s\n', 'regularizePrecisionMatrixConstants_maxExcentricity=9.0' );
fprintf( fileID, '%s\n', 'temporalWindowForLogicalRules=5' );
fprintf( fileID, '%s\n', 'thrBackgroundDetectorHigh=0.7' );
fprintf( fileID, '%s\n', 'thrBackgroundDetectorLow=0.2' );
fprintf( fileID, '%s\n', 'SLD_lengthTMthr=5' );
% Positive  integer  number.  This  parameter  provides  the  radius
% (in  pixels)  of  the  median  filter applied before the watershed
% hierarchical segmentation is performed. The noisier the images, the 
% larger this value should be. 
fprintf( fileID, '%s\n', 'radiusMedianFilter=4' );
% Non-negative floating point value. This parameter provides the minimum
% value of ? used for the hierarchical  segmentation  using
% persistence-based  clustering  of  watershed  regions.  The  higher 
% minTau, the larger the minimum super-voxel size that can be generated at
% the lower level of the hierarchical  segmentation.  This  value  should
% be  kept  low  so  as  not  to  compromise  the framework’s capability t
% recover from under-segmentation. 
fprintf( fileID, '%s\n', 'minTau=4' );
fprintf( fileID, '%s\n', 'conn3D=74' );
fprintf( fileID, '%s\n', 'estimateOpticalFlow=0' );
fprintf( fileID, '%s\n', 'maxDistPartitionNeigh=80.0' );
fprintf( fileID, '%s\n', 'deathThrOpticalFlow=-1' );
% Positive  integer  number.  If  the  number  of  voxels  belonging  to
% a  super-voxel  is  less  than minNucleiSize  the  super-voxel  is
% deleted.  This  parameter  is  useful  to  delete  spurious  super-
% voxels representing background intensity. 
minNucleiSizeStr = strcat( 'minNucleiSize=', num2str( minNucleiSize ) );
fprintf( fileID, '%s\n', char(minNucleiSizeStr) );
% Positive  integer  number.  This  parameter  defines  the  maximum
% allowed  size  (in  voxels)  of  a super-voxel  after  applying  Otsu's
% threshold.  If  Otsu's  threshold  generates  an  object  larger  than 
% maxNucleiSize the threshold is increased until the objet size falls below maxNucleiSize. 
maxNucleiSizeStr = strcat( 'maxNucleiSize=', num2str( maxNucleiSize ) );
fprintf( fileID, '%s\n', char(maxNucleiSizeStr) );
fprintf( fileID, '%s\n', 'maxPercentileTrimSV=0.8' );
fprintf( fileID, '%s\n', 'conn3DsvTrim=6' );
% Positive integer number. This parameter defines the maximum number of
% nearest neighbors to consider for each super-voxel when building the
% spatio-temporal graph for tracking. The shorter the nuclear displacement
% between time points, the lower the parameter value can be.
fprintf( fileID, '%s\n', 'maxNumKNNsupervoxel=10' );
% Floating  point  positive  number.  This  parameter  defines  the 
% maximum  distance  (in  pixels)  to consider for each super-voxel when
% building the spatio-temporal graph for tracking. The shorter the
% nuclear displacement between time points, the lower the parameter value can be. 
fprintf( fileID, '%s\n', 'maxDistKNNsupervoxel=80.0' );
fprintf( fileID, '%s\n', 'thrSplitScore=-1.0' );
% The  feature  calculates  the  distance  (in  pixels)  between  mother
% cell  and  the midplane defined by the two daughter cells. If the value
% is above thrCellDivisionPlaneDistance, the  cell  division  is
% considered  a  false  positive  and  the  linkage  between  mother
% and  furthest daughter is removed.
fprintf( fileID, '%s\n', 'thrCellDivisionPlaneDistance=10.0' );
%fprintf( fileID, '%s\n', 'thrCellDivisionWithTemporalWindow=0.3' );
fclose( fileID );
