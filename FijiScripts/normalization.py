import time
import ij.IJ as IJ
import ij.process.StackStatistics as STAT
import ij.ImageStack as IMGS
from ij.io import FileSaver

startT = 3
endT = 23

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
# assuming that we have 16 bit images
maxI = 0
minI = 65535
# determine min and max intensities values
for i in range(fixStartT, fixEndT+1):		
	if i < 10:
		digit = '00' + str(i)
	elif i < 100:
		digit = '0' + str(i)
	elif i < 1000:
		digit = str(i)

	imagePath = sourceFolder + digit + appendix

	# open image file
	print 'Opening image', imagePath
	imp = IJ.openImage( imagePath )

	stats = STAT( imp )
	minIntens = stats.histMin
	maxIntens = stats.histMax

	if minIntens < minI:
		minI = minIntens

	if maxIntens > maxI:
		maxI = maxIntens

print minI, maxI

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
	slices = imp.getNSlices()
	stack1 = imp.getStack()
	stack2 = IMGS( imp.getWidth(), imp.getHeight() )
	print 'Normalizing image...'
	for s in range( 0, slices ):
		ip = stack1.getProcessor(s+1)
		for x in range( 0, imp.getWidth() ):
			for y in range( 0, imp.getHeight() ):
				oldI = ip.get( x, y )
				newI = 65535 * (oldI - minI)/(maxI - minI)
				newI = int(newI)
				ip.set( x, y, newI )
		label = stack1.getSliceLabel(s+1)
		ip.resetMinAndMax()
		stack2.addSlice( label, ip )
		IJ.showProgress(s+1, slices)
				
	imp.setStack( stack2 )
	fs = FileSaver( imp )
	print "Saving to", outputPath
	fs.saveAsTiff( outputPath )
	