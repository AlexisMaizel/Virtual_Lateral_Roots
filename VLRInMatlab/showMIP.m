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

% if MIPtype == 0
%   dataPathStr = strcat(dataPathStr, dataStr, '\MIP', mipT, '_T', num2str(timestep), '.jpg' );
% else
%   if timestep < 10
%     digit = '00';
%   elseif timestep < 100
%     digit = '0';
%   else
%     digit = '';
%   end
%   dataPathStr = strcat(dataPathStr, dataStr, '\MIP', '_trainedThresholding\MIP_trainedThres_t', digit, num2str(timestep), '.jpg' );
% end

restoredefaultpath
h = imshow( char( dataPathStr ) );
setWorkingPathProperties()