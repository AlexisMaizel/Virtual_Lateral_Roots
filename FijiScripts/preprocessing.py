# preprocessing steps of data resulting in segmented nuclei or membrane BUT they are not labeled
# so generation of connected components etc. is done later in MATLAB
# The following steps can be done:
# - manual or automatic thresholding
# - load a classifier from the trainable WEKA segmentation plugin to segment
# - save of TIFF files and/or MIPs
# - possibility to resample data based on z resolution (but not recommended because of resulting too large data sizes)
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
import ij.plugin.GaussianBlur3D as Gauss
import ij.process.StackStatistics as STAT
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
endT = 10

applyThresholding = 1
# auto, manual or training thresholding: 0, 1, or 2
thresholdingType = 2

saveThresholding = 1
saveRawImageMIP = 0
saveThresholdImageMIP = 1

# save time-dependent stack of MIPs
saveRawMultiImageMIP = 0
saveThresholdMultiImageMIP = 0

# properties if auto local thresholding is selected
useAutoLocalThresholding = 0
thresholdingRadius = 50

# properties if manual thresholding based on percentage of total intensity range
useManualPercentageThresholding = 1
percentage = 0.4

resampleData = 0

# nuclei (=1) or membrane (=0) channel?
nucleiChannel = 1

cropImage = 0
ROI_xStart = 400
ROI_yStart = 600
ROI_xSize = 1100
ROI_ySize = 800

showRawMIP = 1
showThresholdMIP = 1

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
# 20160427
elif chosenData == 6:
	if nucleiChannel == 1:
		sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\red\\cropped_spim_TL'
		outputFolder = 'I:\\SegmentationResults\\Preprocessing\\20160427\\changed_t'
		MIPoutput = MIPoutput + '20160427\\MIP'
		classifierPath = 'I:\\SegmentationResults\\Preprocessing\\20160427\\nucleiClassifierT10-168.model'
	else:
		sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\green\\cropped_spim_TL'
		outputFolder = 'I:\\SegmentationResults\\Preprocessing\\20160427\\changed_membrane_t'
		MIPoutput = MIPoutput + '20160427\\MIP_membrane'
		classifierPath = 'I:\\SegmentationResults\\Preprocessing\\20160427\\membraneClassifierT10-168.model'
		
	appendix = '_Angle1.tif'
	zRatio = 1
	I1 = 1500
	I2 = 1500
	fixStartT = 10
	fixEndT = 168
# 20160428
elif chosenData == 7:
	sourceFolder = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\nuclei\\left\\cropped_Ch1_CamL_T00'
	#sourceFolder = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\nuclei\\right\\_Ch1_CamR_T00'
	outputFolder = 'I:\\SegmentationResults\\Preprocessing\\20160428\\changed_t'
	appendix = '.tif'
	classifierPath1 = 'I:\\SegmentationResults\\Preprocessing\\20160428\\nucleiClassifierT10-21.model'
	classifierPath2 = 'I:\\SegmentationResults\\Preprocessing\\20160428\\nucleiClassifierT22-34.model'
	split = 21
	MIPoutput = MIPoutput + '20160428\\MIP'
	zRatio = 4
	I1 = 800
	I2 = 800
	fixStartT = 10
	fixEndT = 34
# 20160426
else:
	if nucleiChannel == 1:
		sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160426\\red\\cropped_spim_TL'
		if useManualPercentageThresholding == 1 and thresholdingType == 1:
			outputFolder = 'I:\\SegmentationResults\\Preprocessing\\20160426\\manThres_t'
			MIPoutput = MIPoutput + '20160426\\MIP_manThres'
		else:
			outputFolder = 'I:\\SegmentationResults\\Preprocessing\\20160426\\changed_t'
			MIPoutput = MIPoutput + '20160426\\MIP_trainedThres_t'
		
		#classifierPath = 'I:\\SegmentationResults\\Preprocessing\\20160426\\singleNucleiClassifier.model'
		classifierPath = 'I:\\SegmentationResults\\Preprocessing\\20160426\\nucleiClassifierT3-23.model'
		#classifierPath2 = 'I:\\SegmentationResults\\Preprocessing\\20160426\\nucleiClassifierT14-23.model'
		split = 13
	else:
		sourceFolder = 'I:\\NewDatasets\\Zeiss\\20160426\\green\\cropped_spim_TL'
		outputFolder = 'I:\\SegmentationResults\\Preprocessing\\20160426\\changed_membrane_t'
		MIPoutput = MIPoutput + '20160426\\MIP_membrane'
		classifierPath = 'I:\\SegmentationResults\\WEKA\\20160426\\membraneClassifier.model'

	# for not normalized data
	appendix = '_Angle1.tif'
	#appendix = '.tif'
	
	zRatio = 7.5
	I1 = 1000
	I2 = 1000
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
			if useAutoLocalThresholding == 0:
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
			if useManualPercentageThresholding == 0:
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
			else:
				#stats = STAT( imp )
				#minIntens = stats.histMin
				#maxIntens = stats.histMax

				
				#slices = imp.getNSlices()
				#stack = imp.getStack()
				#minIntens = 65535
				#maxIntens = 0
				#for s in range( 0, slices ):
					#ip = stack.getProcessor(s+1)
					#miI = ip.getMin()
					#maI = ip.getMax()
					#print miI, maI
					#if miI < minIntens:
						#minIntens = miI

					#if maI > maxIntens:
						#maxIntens = maI

				
				#print minIntens, maxIntens
				#intens = minIntens + (maxIntens - minIntens) * percentage

				stack = imp.getStack()
				slices = imp.getNSlices()
				meanI = 0
				dim = 0
				for s in range( 0, slices ):
					ip = stack.getProcessor(s+1)
					for x in range( 0, imp.getWidth() ):
						for y in range( 0, imp.getHeight() ):
							val = ip.get( x, y )
							if val > 1000:
								meanI = meanI + val
								dim = dim + 1
					IJ.showProgress(s+1, slices)

				#dim = imp.getWidth()*imp.getHeight()*slices
				intens = float(meanI)/float(dim)
				
			#if endT != startT:
			#	curTRatio = float(i-1)/float(abs( endT - startT ))
			#else:
			#	curTRatio = 0
			#intens = startI - curTRatio*float(abs( endI - startI ))
				
			print 'Threshold:', intens
			newImp = convertToBinaryMask( imp, int(intens) )
			#newImp = cropVoxelValuesBelowThreshold( imp, int(intens) )

			# apply blurring
			#if useManualPercentageThresholding == 1:
			#	Gauss().blur( newImp, 2., 2., 2. )

		# trained segmentation
		else:
			#zp = ZProjector(imp)
			#zp.setMethod( ZProjector.MAX_METHOD )
			#zp.doProjection()
			#imp = zp.getProjection()
			segmentator = WS( imp )

			#if i <= split:
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
	