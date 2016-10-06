import sys
import ij.IJ as IJ
from ij.io import FileSaver
import ij.plugin.ZProjector as ZProjector
import ij.ImagePlus as ImagePlus
sys.path.append( 'C:\\Jens\\VLRRepository\\FijiScripts' )

# data can be in [1,9]
chosenData = 8
applyNucleiCropping = 1
applyMembraneCropping = 0

if chosenData == 6:
	nucleiFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\red\spim_TL'
	membraneFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\green\\spim_TL'
	cropNucleiOutput = 'I:\\NewDatasets\\ilastikWorkshopData\\20160427\\nuclei\\division_cropped_nuclei_T'
	cropMembraneOutput = 'I:\\NewDatasets\\ilastikWorkshopData\\20160427\\membrane\\division_cropped_membrane_T'
	startT = 43
	endT = 50 #168
	#ROI_xStart = 0
	#ROI_yStart = 0
	#ROI_xSize = 640
	#ROI_ySize = 640
	# old before 05/09/2016
	#ROI_xStart = 200
	#ROI_yStart = 410
	#ROI_xSize = 200
	#ROI_ySize = 400
	# single division analysis for small region
	ROI_xStart = 240
	ROI_yStart = 410
	ROI_xSize = 40
	ROI_ySize = 70
	ROI_zStart = 220
	ROI_zSize = 55 # if zero then it will be set later if the image size is read else chosen size is taken
	# single division analysis from t in [43,50]
	zStart = 0
	zLength = 640
	appendix = '_Angle1.tif'
elif chosenData == 7:
	nucleiFolder = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\nuclei\\left\\_Ch1_CamL_T00'
	membraneFolder = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\membrane\\left\\_Ch0_CamL_T00'
	cropNucleiOutput = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\nuclei\\left\\cropped_Ch1_CamL_T00'
	cropMembraneOutput = 'I:\\NewDatasets\\2016-04-28_17.35.59_JENS\\Tiffs\\membrane\\left\\cropped_Ch0_CamL_T00'
	startT = 10
	endT = 34
	ROI_xStart = 500
	ROI_yStart = 350
	ROI_xSize = 1000
	ROI_ySize = 950
	ROI_zStart = 0
	ROI_zSize = 0 # if zero then it will be set later if the image size is read else chosen size is taken
	appendix = '.tif'
elif chosenData == 8:
	nucleiFolder = 'I:\\NewDatasets\\Zeiss\\20160426\\red\\spim_TL'
	membraneFolder = 'I:\\NewDatasets\\Zeiss\\20160426\\green\\spim_TL'
	cropNucleiOutput = 'I:\\NewDatasets\\Zeiss\\20160426\\red\\test_cropped_spim_TL'
	cropMembraneOutput = 'I:\\NewDatasets\\Zeiss\\20160426\\green\\test_cropped_spim_TL'
	startT = 3
	endT = 23
	#ROI_xStart = 400
	#ROI_yStart = 600
	#ROI_xSize = 1100
	#ROI_ySize = 800
	ROI_xStart = 700
	ROI_yStart = 500
	ROI_xSize = 400
	ROI_ySize = 700
	ROI_zStart = 0
	ROI_zSize = 0 # if zero then it will be set later if the image size is read else chosen size is taken
	appendix = '_Angle1.tif'
elif chosenData == 9:
	membraneFolder = 'I:\\NewDatasets\\20160706\\beamExp_Ch2_CamL_T00'
	cropMembraneOutput = 'I:\\NewDatasets\\20160706\\cropped_beamExp_Ch2_CamL_T00'
	startT = 0
	endT = 0
	ROI_xStart = 300
	ROI_yStart = 0
	ROI_xSize = 800
	ROI_ySize = 2048
	ROI_zStart = 0
	ROI_zSize = 0 # if zero then it will be set later if the image size is read else chosen size is taken
	appendix = '.tif'

for t in range(startT, endT+1):
	if t < 10:
		digit = '00' + str(t)
	elif t < 100:
		digit = '0' + str(t)
	elif t < 1000:
		digit = str(t)

	if chosenData < 9 and applyNucleiCropping == 1:
		imageNPath = nucleiFolder + digit + appendix
		imageNCroppedPath = cropNucleiOutput + digit + appendix
	
		print 'Opening image', imageNPath
		impN = IJ.openImage( imageNPath )
		impN.setRoi( ROI_xStart, ROI_yStart, ROI_xSize, ROI_ySize )
		slices = impN.getNSlices()
		if ROI_zSize == 0:
			ROI_zSize = slices
		stack = impN.getStack()
		# crop(int x, int y, int z, int width, int height, int depth)
		impN.setStack( stack.crop( ROI_xStart, ROI_yStart, ROI_zStart, ROI_xSize, ROI_ySize, ROI_zSize ) )
		#impN.setStack( stack.crop( ROI_xStart, ROI_yStart, 130, ROI_xSize, ROI_ySize, 220 ) )
		#impN.setStack( stack.crop( ROI_xStart, ROI_yStart, 275, ROI_xSize, ROI_ySize, 1 ) )
		fs = FileSaver( impN )
		fs.saveAsTiff( imageNCroppedPath )
	
		zp = ZProjector(impN)
		zp.setMethod( ZProjector.MAX_METHOD )
		zp.doProjection()
		MIPimpN = zp.getProjection()
		#MIPimpN.show()

	if applyMembraneCropping == 1:
		imageMPath = membraneFolder + digit + appendix
		imageMCroppedPath = cropMembraneOutput + digit + appendix
	
		print 'Opening image', imageMPath
		impM = IJ.openImage( imageMPath )
		impM.setRoi( ROI_xStart, ROI_yStart, ROI_xSize, ROI_ySize )
		slices = impM.getNSlices()
		if ROI_zSize == 0:
			ROI_zSize = slices
		stack = impM.getStack()
		# crop(int x, int y, int z, int width, int height, int depth)
		impM.setStack( stack.crop( ROI_xStart, ROI_yStart, ROI_zStart, ROI_xSize, ROI_ySize, ROI_zSize ) )
		#impM.setStack( stack.crop( ROI_xStart, ROI_yStart, 275, ROI_xSize, ROI_ySize, 1 ) )
		fs = FileSaver( impM )
		fs.saveAsTiff( imageMCroppedPath )
	
		zp = ZProjector(impM)
		zp.setMethod( ZProjector.MAX_METHOD )
		zp.doProjection()
		MIPimpM = zp.getProjection()
		#MIPimpM.show()
	