import time
import ij.IJ as IJ
from ij.io import FileSaver
import ij.plugin.ZProjector as ZProjector

sourceFolder = 'I:\\SegmentationResults\\Preprocessing\\20160426\\changed_t'
outputFolder = 'I:\\SegmentationResults\\Preprocessing\\20160426\\MIP_trainedThres_t'
appendix = '.tif'

startT = 3
endT = 23

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