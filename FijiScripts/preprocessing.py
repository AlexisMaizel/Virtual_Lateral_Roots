import sys
import time
import ij
import ij.IJ as IJ
from ij.io import FileSaver
#import features
import fiji
import fiji.threshold.Auto_Threshold as AT
import fiji.threshold.Auto_Local_Threshold as ALT
import ij.plugin.ZProjector as ZProjector
import ij.process as process
import ij.process.StackConverter as SCON
import ij.plugin.Converter as CON
import imagescience.image.Image as Image
import imagescience.transform.Scale as Scale
import ij.plugin as plugin
import ij.ImagePlus
import ij.ImageStack as IMGS
import trainableSegmentation.WekaSegmentation as WS
import ij.plugin.Thresholder as THRES
sys.path.append( 'C:\\Jens\\VLRRepository\\FijiScripts' )
from cropVoxelValuesBelowThreshold import *

startT = 10
endT = 168

applyThresholding = 1
# auto, manual or training thresholding: 0, 1, or 2
thresholdingType = 2

saveThresholding = 1
saveRawImageMIP = 1
saveThresholdImageMIP = 1

# save time-dependent stack of MIPs
saveRawMultiImageMIP = 1
saveThresholdMultiImageMIP = 1

autoLocalThresholding = 0
thresholdingRadius = 50

resampleData = 0

cropImage = 0
ROI_xStart = 400
ROI_yStart = 600
ROI_xSize = 1100
ROI_ySize = 800

showRawMIP = 0
showThresholdMIP = 0

MIPoutput = 'I:\\SegmentationResults\\Preprocessing\\'

# data can be in [1,8]
chosenData = 6
# 120830
if chosenData == 1:
	sourceFolder = 'I:\\FrankfurtLSFMDatasets\\20120830-pGATA_H2B_Wave\\AllTimePoints\\Not_Registered\\cropped_pGata23_120830_TL0'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\120830\\changed_t'
	appendix = '_CHN00.tif'
	MIPoutput = MIPoutput + '120830\\MIP'
	# NOTE: data quality too low
	zRatio = 2
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
	zRatio = 2
	I1 = 550
	I2 = 550
	I3 = 550
	I4 = 400
# 121211
elif chosenData == 3:
	sourceFolder = 'I:\\FrankfurtLSFMDatasets\\20121211_pGATA_H2B_Wave\\3Ddrift_Stacks_cropped\\cropped_t'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\121211\\changed_t'
	appendix = '_c1.tif'
	classifierPath = 'I:\\SegmentationResults\\Preprocessing\\121211\\nucleiClassifierT1-50.model'
	MIPoutput = MIPoutput + '121211\\MIP'
	zRatio = 2
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
	zRatio = 2
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
	zRatio = 2
	# NOTE: data quality too low
	I1 = 900
	I2 = 900
	I3 = 900
	I4 = 900
elif chosenData == 6:
	#sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\green\\spim_TL'
	sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\red\\cropped_spim_TL'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\20160427\\changed_t'
	appendix = '_Angle1.tif'
	classifierPath = 'I:\\SegmentationResults\\Preprocessing\\20160427\\nucleiClassifierT10-168.model'
	MIPoutput = MIPoutput + '20160427\\MIP'
	zRatio = 1
	I1 = 950
	I2 = 900
	fixStartT = 10
	fixEndT = 168
elif chosenData == 7:
	sourceFolder = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\nuclei\\left\\cropped_Ch1_CamL_T00'
	#sourceFolder = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\nuclei\\right\\_Ch1_CamR_T00'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\20160428\\changed_t'
	appendix = '.tif'
	MIPoutput = MIPoutput + '20160428\\MIP'
	zRatio = 4
	I1 = 800
	I2 = 800
	fixStartT = 10
	fixEndT = 34
else:
	sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160426\\red\\cropped_spim_TL'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\20160426\\changed_t'
	appendix = '_Angle1.tif'
	classifierPath1 = 'I:\\SegmentationResults\\Preprocessing\\20160426\\nucleiClassifierT3-23.model'
	classifierPath2 = 'I:\\SegmentationResults\\Preprocessing\\20160426\\nucleiClassifierT3-23_2.model'
	MIPoutput = MIPoutput + '20160426\\MIP'
	zRatio = 7.5
	I1 = 900
	I2 = 1100
	fixStartT = 3
	fixEndT = 23

t = time.clock()
for i in range(startT, endT+1):		
	if i < 10:
		digit = '00' + str(i)
	elif i < 100:
		digit = '0' + str(i)
	elif i < 1000:
		digit = str(i)

	imagePath = sourceFolder + digit + appendix
	output = outputFolder + digit + '.tif'

	# set intensity ranges with respect to the current time step i
	if chosenData < 6:
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
	else:
		startI = I1
		endI = I2

	# open image file
	print 'Opening image', imagePath
	imp = IJ.openImage( imagePath )

	if cropImage == 1:
		if chosenData == 7 or chosenData == 8:
			imp.setRoi( ROI_xStart, ROI_yStart, ROI_xSize, ROI_ySize )
			slices = imp.getNSlices()
			stack = imp.getStack()
			imp.setStack( stack.crop( ROI_xStart, ROI_yStart, 0, ROI_xSize, ROI_ySize, slices ) )

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
		if thresholdingType == 0:
			if autoLocalThresholding == 0:
				# arguments are: ImagePlus imp, String myMethod, boolean noWhite, boolean noBlack, boolean doIwhite,
				# boolean doIset, boolean doIlog , boolean doIstackHistogram
				methods = [ "Default", "Huang", "Intermodes", "IsoData",  "Li", "MaxEntropy","Mean", "MinError(I)", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag" , "Triangle", "Yen" ]
				met = methods[10]
				print 'Used thresholding method:', met
				ret = AT().exec( imp, met, True, True, True, False, False, True )
				#ret = at.exec( imp, methods[len(methods)-1], True, True, True, False, False, True )
				newImp = ret[1]
				thres = ret[0]
				# print threshold
				print 'Threshold:', thres
				#newImp.setTitle( methods[len(methods)-1] + ' threshold: ' + str(thres) )
			else:
				methods = [ "Bernsen", "Contrast", "Mean", "Median", "MidGrey", "Niblack", "Otsu", "Phansalkar", "Sauvola" ]
				met = methods[7]
				print 'Used local thresholding method:', met
				#CON().run("8-bit")
				#SCON( imp ).convertToGray8()
				#aImp = ALT().exec( imp, met, thresholdingRadius, 0, 0, True )
				#newImp = aImp[0]
				#SCON( newImp ).convertToGray16()
		elif thresholdingType == 1:
			if chosenData < 6:
				iFactor = i%50
				if iFactor == 0:
					curTRatio = 0
				else:
					curTRatio = float(iFactor-1)/float(50)
				intens = curTRatio * endI + (1-curTRatio) * startI
			else:
				curTRatio = float(abs( fixEndT - i ))/float(abs( fixEndT - fixStartT ))
				intens = curTRatio * startI + (1-curTRatio) * endI
				
			#if endT != startT:
			#	curTRatio = float(i-1)/float(abs( endT - startT ))
			#else:
			#	curTRatio = 0
			#intens = startI - curTRatio*float(abs( endI - startI ))
				
			print 'Threshold:', intens
			newImp = convertToBinaryMask( imp, int(intens) )
			#newImp = cropVoxelValuesBelowThreshold( imp, intens )
		# trained segmentation
		else:
			#zp = ZProjector(imp)
			#zp.setMethod( ZProjector.MAX_METHOD )
			#zp.doProjection()
			#imp = zp.getProjection()
			segmentator = WS( imp )
			#if i < 19:
				#segmentator.loadClassifier( classifierPath1 )
			#else:
				#segmentator.loadClassifier( classifierPath2 )
			segmentator.loadClassifier( classifierPath )
				
			newImp = segmentator.applyClassifier( imp )
			#SCON( newImp ).convertToGray8()
			slices = newImp.getNSlices()
			stack1 = newImp.getStack()
			stack2 = IMGS( newImp.getWidth(), newImp.getHeight() )
			for s in range( 0, slices ):
				ip = stack1.getProcessor(s+1)
				ip.invert()
				label = stack1.getSliceLabel(s+1)
				ip.resetMinAndMax()
				stack2.addSlice( label, ip.convertToByte(True) )
				
			newImp.setStack( stack2 )

		if resampleData == 1:
			scaler = Scale()
			newImg = Image.wrap( newImp )
			sImg = scaler.run( newImg, 1., 1., zRatio, 1., 1., Scale.LINEAR )
			newImp = sImg.imageplus()

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
			print "Saving to", output
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
	