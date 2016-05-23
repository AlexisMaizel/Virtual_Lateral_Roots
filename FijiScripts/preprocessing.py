import sys
import time
import ij
import ij.IJ as IJ
from ij.io import FileSaver
#import features
import fiji
#import fiji.threshold as threshold
import ij.plugin.ZProjector as ZProjector
import ij.process as process
import ij.plugin as plugin
import ij.ImagePlus
sys.path.append( 'C:\\Jens\\VLRRepository\\FijiScripts' )
from cropVoxelValuesBelowThreshold import *

startT = 0
endT = 34

applyThresholding = 1
manualThresholding = 1
saveThresholding = 1
saveRawMultiImageMIP = 1
saveThresholdMultiImageMIP = 1
saveRawImageMIP = 1
saveThresholdImageMIP = 1
cropImage = 1

showRawMIP = 0
showThresholdMIP = 0

MIPoutput = 'I:\\SegmentationResults\\Preprocessing\\'

# data can be in [1,7]
chosenData = 7
# 120830
if chosenData == 1:
	sourceFolder = 'I:\\FrankfurtLSFMDatasets\\20120830-pGATA_H2B_Wave\\AllTimePoints\\Not_Registered\\cropped_pGata23_120830_TL0'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\120830\\changed_t'
	appendix = '_CHN00.tif'
	MIPoutput = MIPoutput + '120830\\MIP'
	# NOTE: data quality too low
	I1 = 900
	I2 = 900
	I3 = 900
	I4 = 900
# 121204
elif chosenData == 2:
	sourceFolder = 'I:\\FrankfurtLSFMDatasets\\20121204_pGATA_H2B_Wave\\driftcorrected_stacks_cropped\\ch01\\crop_t'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\121204\\changed_t'
	appendix = '_c1.tif'
	MIPoutput = MIPoutput + '121204\\MIP'
	I1 = 550
	I2 = 550
	I3 = 550
	I4 = 400
# 121211
elif chosenData == 3:
	sourceFolder = 'I:\\FrankfurtLSFMDatasets\\20121211_pGATA_H2B_Wave\\3Ddrift_Stacks_cropped\\cropped_t'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\121211\\changed_t'
	appendix = '_c1.tif'
	MIPoutput = MIPoutput + '121211\\MIP'
	I1 = 900
	I2 = 900
	I3 = 800
	I4 = 690
# 130508
elif chosenData == 4:
	sourceFolder = 'I:\\FrankfurtLSFMDatasets\\130508\\crop_t'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\130508\\changed_t'
	appendix = '_c1.tif'
	MIPoutput = MIPoutput + '130508\\MIP'
	I1 = 850
	I2 = 850
	I3 = 750
	I4 = 700
# 130607
elif chosenData == 5:
	sourceFolder = 'I:\\FrankfurtLSFMDatasets\\130607\\c1\\crop_t'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\130607\\changed_t'
	appendix = '_c1.tif'
	MIPoutput = MIPoutput + '130607\\MIP'
	# NOTE: data quality too low
	I1 = 900
	I2 = 900
	I3 = 900
	I4 = 900
elif chosenData == 6:
	#sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\green\\spim_TL'
	sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\red\\spim_TL'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\20160427\\changed_t'
	appendix = '_Angle1.tif'
	MIPoutput = MIPoutput + '20160427\\MIP'
	I1 = 950
	I2 = 900
	I3 = 900
	I4 = 900
else:
	sourceFolder = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\nuclei\\left\\_Ch1_CamL_T00'
	#sourceFolder = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\nuclei\\right\\_Ch1_CamR_T00'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\20160428\\changed_t'
	appendix = '.tif'
	MIPoutput = MIPoutput + '20160428\\MIP'
	I1 = 800
	I2 = 800
	I3 = 800
	I4 = 800

t = time.clock()
for i in range(startT, endT+1):		
	if i < 10:
		imagePath = sourceFolder + '00' + str(i) + appendix
		output = outputFolder + '00' + str(i) + '.tif'
	elif i < 100:
		imagePath = sourceFolder + '0' + str(i) + appendix
		output = outputFolder + '0' + str(i) + '.tif'
	elif i < 1000:
		imagePath = sourceFolder + str(i) + appendix
		output = outputFolder + str(i) + '.tif'

	# set intensity ranges with respect to the current time step i
	if i < 50:
		startI = I1
		endI = I2
	elif i < 100:
		startI = I2
		endI = I3
	elif i < 150:
		startI = I3
		endI = I4
	else:
		startI = I4
		endI = I4

	# open image file
	print 'Opening image', imagePath
	imp = IJ.openImage( imagePath )

	if cropImage == 1:
		if chosenData == 7:
			imp.setRoi( 650, 0, 700, 2048 )
			slices = imp.getNSlices()
			stack = imp.getStack()
			imp.setStack( stack.crop( 650, 0, 0, 700, 2048, slices ) )

	if i == startT and (saveRawMultiImageMIP == 1 or saveThresholdMultiImageMIP == 1):
		width = imp.getWidth()
		height = imp.getHeight()
		rawImageMipStack = ij.ImageStack( width, height, endT-startT+1 )
		thresholdImageMipStack = ij.ImageStack( width, height, endT-startT+1 )
		stackCounter = 1
	
	if showRawMIP == 1 or saveRawMultiImageMIP == 1 or saveRawImageMIP == 1:
		zp = ZProjector(imp)
		zp.setMethod( ZProjector.MAX_METHOD )
		zp.doProjection()
		MIPimp = zp.getProjection()
		if showRawMIP == 1:
			MIPimp.show()
			
		if saveRawMultiImageMIP == 1:
			rawImageMipStack.setProcessor( MIPimp.getProcessor(), stackCounter )

		if saveRawImageMIP == 1:
			fs = FileSaver( MIPimp )
			rawSingleMIPOutput = MIPoutput + 'raw_T' + str(i) + '.jpg'
			fs.saveAsJpeg( rawSingleMIPOutput )

	if applyThresholding == 1:
		if manualThresholding == 1:
			iFactor = i%50
			if iFactor == 0:
				curTRatio = 0
			else:
				curTRatio = float(iFactor-1)/float(50)
			intens = curTRatio * endI + (1-curTRatio) * startI
			
			#if endT != startT:
			#	curTRatio = float(i-1)/float(abs( endT - startT ))
			#else:
			#	curTRatio = 0
			#intens = startI - curTRatio*float(abs( endI - startI ))
			
			print 'Threshold:', intens
			newImp = convertToBinaryMask( imp, int(intens) )
			#newImp = cropVoxelValuesBelowThreshold( imp, intens )
		else:
			# arguments are: ImagePlus imp, String myMethod, boolean noWhite, boolean noBlack, boolean doIwhite,
			# boolean doIset, boolean doIlog , boolean doIstackHistogram
			methods = [ "Default", "Huang", "Intermodes", "IsoData",  "Li", "MaxEntropy","Mean", "MinError(I)", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag" , "Triangle", "Yen" ]
			at = threshold.Auto_Threshold()
			met = methods[15]
			print 'Used thresholding method:', met
			ret = at.exec( imp, met, True, True, True, False, False, True )
			#ret = at.exec( imp, methods[len(methods)-1], True, True, True, False, False, True )
			newImp = ret[1]
			thres = ret[0]
			# print threshold
			print 'Threshold:', thres
			#newImp.setTitle( methods[len(methods)-1] + ' threshold: ' + str(thres) )

		# z project image stack
		if showThresholdMIP == 1 or saveThresholdMultiImageMIP == 1 or saveThresholdImageMIP == 1:
			zp = ZProjector(newImp)
			zp.setMethod( ZProjector.MAX_METHOD )
			zp.doProjection()
			newMIPImp = zp.getProjection()
			if showThresholdMIP == 1:
				newMIPImp.show()
				
			if saveThresholdMultiImageMIP == 1:
				thresholdImageMipStack.setProcessor( newMIPImp.getProcessor(), stackCounter )

			if saveThresholdImageMIP == 1:
				fs = FileSaver( newMIPImp )
				thresholdSingleMIPOutput = MIPoutput + 'threshold_T' + str(i) + '.jpg'
				fs.saveAsJpeg( thresholdSingleMIPOutput )

		# save image file
		#print 'Saving result to', output
		if saveThresholding == 1:
			fs = FileSaver( newImp )
			fs.saveAsTiff( output )

	if saveThresholdMultiImageMIP == 1 or saveRawMultiImageMIP == 1:
		stackCounter = stackCounter + 1

	# process with tubeness
	#print 'Applying tubeness'
	#sigma = 6.
	#useCalibration = True
	#tp = features.TubenessProcessor( sigma, useCalibration )
	#newImp = tp.generateImage( newImp )
	
	# convert to 16 bit image
	#print 'Converting to 16 bit'
	#con = process.StackConverter(newImp)
	#con.convertToGray16()

if saveRawMultiImageMIP == 1:
	rawMIPimp = ij.ImagePlus()
	rawMIPimp.setStack( rawImageMipStack )
	fsRM = FileSaver( rawMIPimp )
	rawMIPOutput = MIPoutput + 'raw_startT' + str(startT) + '_endT' + str(endT) + '.tif'
	fsRM.saveAsTiff( rawMIPOutput )

if saveThresholdMultiImageMIP == 1:
	thresholdMIPimp = ij.ImagePlus()
	thresholdMIPimp.setStack( thresholdImageMipStack )
	fsTM = FileSaver( thresholdMIPimp )
	thresholdMIPOutput = MIPoutput + 'threshold_startT' + str(startT) + '_endT' + str(endT) + '.tif'
	fsTM.saveAsTiff( thresholdMIPOutput )
	
print (time.clock() - t)/60., 'minutes'
print 'Done'
	