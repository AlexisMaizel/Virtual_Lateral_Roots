# open TIFF files and store them as MIPs
import time
import ij.IJ as IJ
from ij.io import FileSaver
import ij.plugin.ZProjector as ZProjector

chosenData = 8 #[6,8]
if chosenData == 6:
	sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\red\\spim_TL'
	outputFolder = 'I:\\SegmentationResults\\MIPsRawData\\20160427\\MIP_t'
	startT = 10
	endT = 168
	appendix = '_Angle1.tif'
elif chosenData == 7:
	sourceFolder = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\nuclei\\left\\_Ch1_CamL_T00'
	outputFolder = 'I:\\SegmentationResults\\MIPsRawData\\20160428\\MIP_t'
	startT = 10
	endT = 34
	appendix = '.tif'
elif chosenData == 8:
	#sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160426\\red\\partiallyNormalized\\cropped_spim_TL'
	sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160426\\red\\test_cropped_spim_TL'
	outputFolder = 'I:\\SegmentationResults\\MIPsRawData\\20160426\\MIP_small_cropped_t'
	startT = 3
	endT = 23
	appendix = '_Angle1.tif'

for i in range(startT, endT+1):		
	if i < 10:
		digit = '00' + str(i)
	elif i < 100:
		digit = '0' + str(i)
	elif i < 1000:
		digit = str(i)

	imagePath = sourceFolder + digit + appendix
	outputPath = outputFolder + digit + '.jpg'

	# open image file
	print 'Opening image', imagePath
	imp = IJ.openImage( imagePath )

	zp = ZProjector(imp)
	zp.setMethod( ZProjector.MAX_METHOD )
	zp.doProjection()
	MIPimp = zp.getProjection()
	fs = FileSaver( MIPimp )
	fs.saveAsJpeg( outputPath )

print 'Done'