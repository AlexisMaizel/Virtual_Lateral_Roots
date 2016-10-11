# Testing of using ultimate erosion in Fiji which has different results compared to the results
# in MATLAB
import time
import ij.IJ as IJ
import ij.process.StackStatistics as STAT
import ij.ImageStack as IMGS
from ij.io import FileSaver
import ij.plugin.filter.BackgroundSubtracter as TOPHAT
import ij.plugin.filter.EDM as EDM
import ij.ImagePlus as IP
import ij.plugin.ZProjector as ZProjector
import sys
sys.path.append( 'C:\\Jens\\VLRRepository\\FijiScripts' )
from cropVoxelValuesBelowThreshold import *

startT = 3
endT = 3

# data can be in [1,8]
chosenData = 8
# 120830
if chosenData == 1:
	sourceFolder = 'I:\\FrankfurtLSFMDatasets\\20120830-pGATA_H2B_Wave\\AllTimePoints\\Not_Registered\\cropped_pGata23_120830_TL0'
	outputFolder = 'I:\\FrankfurtLSFMDatasets\\20120830-pGATA_H2B_Wave\\AllTimePoints\\Not_Registered\\normalized\\cropped_norm_pGata23_120830_TL0'
	appendix = '_CHN00.tif'
# 121204
elif chosenData == 2:
	sourceFolder = 'I:\\FrankfurtLSFMDatasets\\20121204_pGATA_H2B_Wave\\driftcorrected_stacks_cropped\\ch01\\crop_t'
	outputFolder = 'I:\\FrankfurtLSFMDatasets\\20121204_pGATA_H2B_Wave\\driftcorrected_stacks_cropped\\ch01\\normalized\\crop_norm_t'
	appendix = '_c1.tif'
# 121211
elif chosenData == 3:
	sourceFolder = 'I:\\FrankfurtLSFMDatasets\\20121211_pGATA_H2B_Wave\\3Ddrift_Stacks_cropped\\cropped_t'
	outputFolder = 'I:\\FrankfurtLSFMDatasets\\20121211_pGATA_H2B_Wave\\3Ddrift_Stacks_cropped\\normalized\\cropped_norm_t'
	appendix = '_c1.tif'
# 130508
elif chosenData == 4:
	sourceFolder = 'I:\\FrankfurtLSFMDatasets\\130508\\crop_t'
	outputFolder = 'I:\\FrankfurtLSFMDatasets\\130508\\normalized\\crop_norm_t'
	appendix = '_c1.tif'
# 130607
elif chosenData == 5:
	sourceFolder = 'I:\\FrankfurtLSFMDatasets\\130607\\c1\\crop_t'
	outputFolder = 'I:\\FrankfurtLSFMDatasets\\130607\\c1\\normalized\\crop_norm_t'
	appendix = '_c1.tif'
# 20160427
elif chosenData == 6:
	sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\red\\cropped_spim_TL'
	outputFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\red\\normalized\\cropped_norm_spim_TL'
	appendix = '_Angle1.tif'
	fixStartT = 10
	fixEndT = 168
# 20160428
elif chosenData == 7:
	sourceFolder = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\nuclei\\left\\cropped_Ch1_CamL_T00'
	outputFolder = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\nuclei\\left\\normalized\\cropped_norm_Ch1_CamL_T00'
	appendix = '.tif'
	I2 = 800
	fixStartT = 10
	fixEndT = 34
# 20160426
else:
	sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160426\\red\\cropped_spim_TL'
	outputFolder = 'I:\\NewDatasets\\Zeiss\\20160426\\red\\normalized\\cropped_norm_spim_TL'
	appendix = '_Angle1.tif'
	fixStartT = 3
	fixEndT = 23

t = time.clock()

# normalize each data set based on minI and maxI values
for i in range(startT, endT+1):		
	if i < 10:
		digit = '00' + str(i)
	elif i < 100:
		digit = '0' + str(i)
	elif i < 1000:
		digit = str(i)

	imagePath = sourceFolder + digit + appendix
	outputPath = outputFolder + digit + '.tif'

	# open image file
	print 'Opening image', imagePath
	imp = IJ.openImage( imagePath )

	# removing background noise
	slices = imp.getNSlices()
	stack = imp.getStack()
	newStack = IMGS( imp.getWidth(), imp.getHeight() )
	meanI = 0
	dim = 0
	print 'Applying top hat filter...'
	for s in range( 0, slices ):
		ip = stack.getProcessor(s+1)
		# ip, double radius, boolean createBackground, boolean lightBackground, boolean useParaboloid,
		# boolean doPresmooth, boolean correctCorners
		TOPHAT().rollingBallBackground( ip, 20, False, False, False, False, True )

		for x in range( 0, imp.getWidth() ):
			for y in range( 0, imp.getHeight() ):
				val = ip.get( x, y )
				if val > 500:
					meanI = meanI + val
					dim = dim + 1
		IJ.showProgress(s+1, slices)
	
	#dim = imp.getWidth()*imp.getHeight()*slices
	intens = float(meanI)/float(dim)
	print 'Applying thresholding with intensity =', intens
	newImp = convertToBinaryMask( imp, int(intens) )
	#newImp.show()

	# TODO
	#zp = ZProjector(newImp)
	#zp.setMethod( ZProjector.MAX_METHOD )
	#zp.doProjection()
	#newMIPImp = zp.getProjection()
	#newMIPImp.show()
	#ipMIP = newMIPImp.getProcessor()
	#EDM().toEDM( ipMIP.convertToByte(True) )
	#imgMIP = IP( "MIP", ipMIP )
	#imgMIP.show()
	#break

	print 'Converting to 8 bit...'
	stack = newImp.getStack()
	for s in range( 0, slices ):
		ip = stack.getProcessor(s+1)
		label = stack.getSliceLabel(s+1)
		ip.resetMinAndMax()
		newStack.addSlice( label, ip.convertToByte(True) )
		IJ.showProgress(s+1, slices)

	print 'Applying ultimate erosion...'
	EDMStack = IMGS( imp.getWidth(), imp.getHeight() )
	for s in range( 0, slices ):
		ip = newStack.getProcessor(s+1)
		EDM().toEDM( ip )
		#EDM().run( ip )
		EDMStack.addSlice( label, ip )
		IJ.showProgress(s+1, slices)

	imp.setStack( EDMStack )
	imp.show()
	
	#fs = FileSaver( imp )
	#print "Saving to", outputPath
	#fs.saveAsTiff( outputPath )
	