function [curS, curCellCenters] = applyCellSegmentation( t, minVoxelCount,...
  voxelCountPerCell, inputPath, storeTIFF )
cellRadius = 15;
%   if endT ~= startT
% 		curTRatio = double(t-1)/double(abs( endT - startT ));
%   else
% 		curTRatio = 0;
%   end

%voxelCountPerCell = voxelCountPerCellStart + curTRatio*double(abs( voxelCountPerCellEnd - voxelCountPerCellStart ));
%msg = strcat( msg, {' '}, 'VoxelSize', {' '}, num2str(voxelCountPerCell) );
%disp( msg );

% apply cell segmentation based only on current time step
if t < 10
  digit = '00';
elseif t < 100
  digit = '0';
else
  digit = '';
end
fileName = strcat( inputPath, '\changed_t', digit, num2str(t), '.tif' );
%fileName = strcat( inputPath, '\manThres_t', digit, num2str(t), '.tif' );
%fileName = strcat( 'I:\NewDatasets\Zeiss\20160427\red\cropped_spim_TL', digit, num2str(t), '_Angle1.tif' );
membraneFileName = strcat( inputPath, '\changed_membrane_t', digit, num2str(t), '.tif' );
newFileName = strcat( inputPath, '\preprocessed_t', digit, num2str(t), '.tif' );
[curS, curCellCenters, curMaxIntensities, width, height, slices] =...
  identifyCellObjects( fileName, membraneFileName, minVoxelCount, voxelCountPerCell );
if storeTIFF == 1
  generateCellShape( width, height, slices, curCellCenters, curMaxIntensities, cellRadius, newFileName );
end


