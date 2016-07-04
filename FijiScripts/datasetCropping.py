import sys
import ij.IJ as IJ
from ij.io import FileSaver
import ij.plugin.ZProjector as ZProjector
import ij.ImagePlus as ImagePlus
sys.path.append( 'C:\\Jens\\VLRRepository\\FijiScripts' )

# data can be in [1,8]
chosenData = 6

if chosenData == 6:
	nucleiFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\red\spim_TL'
	membraneFolder = 'I:\\NewDatasets\\Zeiss\\20160427\\green\\spim_TL'
	cropNucleiOutput = 'I:\\NewDatasets\\Zeiss\\20160427\\red\\cropped_spim_TL'
	cropMembraneOutput = 'I:\\NewDatasets\\Zeiss\\20160427\\green\\cropped_spim_TL_slice'
	startT = 10
	endT = 10 #168
	ROI_xStart = 1
	ROI_yStart = 1
	ROI_xSize = 640
	ROI_ySize = 640
	#ROI_xStart = 180
	#ROI_yStart = 300
	#ROI_xSize = 120
	#ROI_ySize = 150
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
	appendix = '.tif'
elif chosenData == 8:
	nucleiFolder = 'I:\\NewDatasets\\Zeiss\\20160426\\red\\spim_TL'
	membraneFolder = 'I:\\NewDatasets\\Zeiss\\20160426\\green\\spim_TL'
	cropNucleiOutput = 'I:\\NewDatasets\\Zeiss\\20160426\\red\\cropped_spim_TL'
	cropMembraneOutput = 'I:\\NewDatasets\\Zeiss\\20160426\\green\\cropped_spim_TL'
	startT = 1
	endT = 23
	ROI_xStart = 400
	ROI_yStart = 600
	ROI_xSize = 1100
	ROI_ySize = 800
	appendix = '_Angle1.tif'

for t in range(startT, endT+1):
	if t < 10:
		digit = '00' + str(t)
	elif t < 100:
		digit = '0' + str(t)
	elif t < 1000:
		digit = str(t)
	
	imageNPath = nucleiFolder + digit + appendix
	imageMPath = membraneFolder + digit + appendix
	imageNCroppedPath = cropNucleiOutput + digit + appendix
	imageMCroppedPath = cropMembraneOutput + digit + appendix

	print 'Opening image', imageNPath
	impN = IJ.openImage( imageNPath )
	impN.setRoi( ROI_xStart, ROI_yStart, ROI_xSize, ROI_ySize )
	slices = impN.getNSlices()
	stack = impN.getStack()
	# crop(int x, int y, int z, int width, int height, int depth)
	impN.setStack( stack.crop( ROI_xStart, ROI_yStart, 0, ROI_xSize, ROI_ySize, slices ) )
	#fs = FileSaver( impN )
	#fs.saveAsTiff( imageNCroppedPath )

	print 'Opening image', imageMPath
	impM = IJ.openImage( imageMPath )
	impM.setRoi( ROI_xStart, ROI_yStart, ROI_xSize, ROI_ySize )
	slices = impM.getNSlices()
	stack = impM.getStack()
	# crop(int x, int y, int z, int width, int height, int depth)
	impM.setStack( stack.crop( ROI_xStart, ROI_yStart, 275, ROI_xSize, ROI_ySize, 1 ) )
	fs = FileSaver( impM )
	fs.saveAsTiff( imageMCroppedPath )

	zp = ZProjector(impM)
	zp.setMethod( ZProjector.MAX_METHOD )
	zp.doProjection()
	MIPimpM = zp.getProjection()
	MIPimpM.show()
	