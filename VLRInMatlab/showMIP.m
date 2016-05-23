function h = showMIP( dataStr, MIPtype, timestep )
dataPathStr = 'I:\SegmentationResults\Preprocessing\';
% raw MIP
if MIPtype == 0
  mipT = 'raw';
 % threshold MIP
else
  mipT = 'threshold';
end
dataPathStr = strcat(dataPathStr, dataStr, '\MIP', mipT, '_T', num2str(timestep), '.jpg' );

restoredefaultpath
h = imshow( char( dataPathStr ) );
setWorkingPathProperties()