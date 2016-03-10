function exportTGMMConfigFile( dataStr,...
  anisotropyZ, backgroundThreshold, persistanceSegmentationTau )
mkdir( 'AutoTGMMConfig\' );
fileName = strcat( 'AutoTGMMConfig\TGMM_configFile', dataStr, '.txt' );
fileID = fopen( char(fileName), 'w' );
dataPathStr = 'imgFilePattern=I:\FrankfurtLSFMDatasets\';
if strcmp(dataStr, '120830')
  dataPathStr = strcat(dataPathStr, '20120830-pGATA_H2B_Wave\AllTimePoints\Registered\CH01\cropped_120830TM30162244SPC001TL????ANG000FRQ000PH00CM0CHN01.Resampled_reg');
elseif strcmp(dataStr, '121204')
  dataPathStr = strcat(dataPathStr, '20121204_pGATA_H2B_Wave\driftcorrected_stacks_cropped\ch01\crop_t???_c1');
elseif strcmp(dataStr, '121211')
  dataPathStr = 'imgFilePattern=I:\SegmentationResults\Preprocessing\121211\changed_t???_c1';
  %dataPathStr = strcat(dataPathStr, '20121211_pGATA_H2B_Wave\3Ddrift_Stacks_cropped\cropped_t???_c1');
elseif strcmp(dataStr, '130508')
  dataPathStr = strcat(dataPathStr, '130508\crop_t???_c1');
elseif strcmp(dataStr, '130607')
  dataPathStr = strcat(dataPathStr, '130607\c1\crop_t???_c1');
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
fprintf( fileID, '%s\n', 'betaPercentageOfN_k=1.0' );
fprintf( fileID, '%s\n', 'nuPercentageOfN_k=1.0' );
fprintf( fileID, '%s\n', 'alphaPercentage=5.0' );
fprintf( fileID, '%s\n', 'maxIterEM=100' );
fprintf( fileID, '%s\n', 'tolLikelihood=1e-6' );
fprintf( fileID, '%s\n', 'regularizePrecisionMatrixConstants_lambdaMin=0.002' );
fprintf( fileID, '%s\n', 'regularizePrecisionMatrixConstants_lambdaMax=5.0' );
fprintf( fileID, '%s\n', 'regularizePrecisionMatrixConstants_maxExcentricity=10.0' );
fprintf( fileID, '%s\n', 'temporalWindowForLogicalRules=5' );
fprintf( fileID, '%s\n', 'thrBackgroundDetectorHigh=0.8' );
fprintf( fileID, '%s\n', 'thrBackgroundDetectorLow=0.2' );
fprintf( fileID, '%s\n', 'SLD_lengthTMthr=5' );
fprintf( fileID, '%s\n', 'radiusMedianFilter=2' );
fprintf( fileID, '%s\n', 'minTau=2' );
fprintf( fileID, '%s\n', 'conn3D=74' );
fprintf( fileID, '%s\n', 'estimateOpticalFlow=0' );
fprintf( fileID, '%s\n', 'maxDistPartitionNeigh=80.0' );
fprintf( fileID, '%s\n', 'deathThrOpticalFlow=-1' );
fprintf( fileID, '%s\n', 'minNucleiSize=5' );
fprintf( fileID, '%s\n', 'maxNucleiSize=16000' );
fprintf( fileID, '%s\n', 'maxPercentileTrimSV=0.4' );
fprintf( fileID, '%s\n', 'conn3DsvTrim=6' );
fprintf( fileID, '%s\n', 'maxNumKNNsupervoxel=10' );
fprintf( fileID, '%s\n', 'maxDistKNNsupervoxel=41.0' );
fprintf( fileID, '%s\n', 'thrSplitScore=-1.0' );
fprintf( fileID, '%s\n', 'thrCellDivisionPlaneDistance=3.2' );
fprintf( fileID, '%s\n', 'thrCellDivisionWithTemporalWindow=0.3' );
fclose( fileID );
