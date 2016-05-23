import sys
import ij.IJ as IJ
import ij.plugin.ZProjector as ZProjector
import ij.ImagePlus as ImagePlus
sys.path.append( 'C:\\Jens\\VLRRepository\\FijiScripts' )

t = 1
s = 187
chosenData = 7

if chosenData == 3:
	dataPath = 'I:\\SegmentationResults\\RACE\\item_0022_SliceBySliceFusionFilter\\cropped_t'
	dataPath2D = 'I:\\SegmentationResults\\RACE\\item_0019_NScaleMorphologicalWatershedFilter\\cropped_t'
	appendix = '_c1_SliceBySliceFusionFilter_Out1.tif'
elif chosenData == 6:
	dataPath = 'I:\\SegmentationResults\\RACE\\item_0022_SliceBySliceFusionFilter\\spim_TL'
	dataPath2D = 'I:\\SegmentationResults\\RACE\\item_0019_NScaleMorphologicalWatershedFilter\\spim_TL'
	appendix = '_Angle1_SliceBySliceFusionFilter_Out1.tif'
elif chosenData == 7:
	dataPath = 'I:\\SegmentationResults\\RACE\\item_0022_SliceBySliceFusionFilter\\spim_TL'
	dataPath2D = 'I:\\SegmentationResults\\RACE\\item_0019_NScaleMorphologicalWatershedFilter\\spim_TL'
	appendix = '_Angle1_SliceBySliceFusionFilter_Out1.tif'

if t < 10:
	digit = '00' + str(t)
elif t < 100:
	digit = '0' + str(t)
elif t < 1000:
	digit = str(t)

imagePath = dataPath + digit + appendix
imagePath2D = dataPath2D + digit + appendix
	
imp = IJ.openImage( imagePath )
imp2D = IJ.openImage( imagePath2D )
#imp.show()
# MIP
if imp:
	zp = ZProjector(imp)
	zp.setMethod( ZProjector.MAX_METHOD )
	zp.doProjection()
	MIPimp = zp.getProjection()
	MIPimp.show()

# specific slice
if imp2D:
	stack = imp2D.getStack()
	slices = imp2D.getNSlices()
	if s < slices:
		ip = stack.getProcessor(s+1)
		impS = ImagePlus()
		impS.setProcessor( 'MIP_Slice' + str(s), ip )
		zpS = ZProjector(impS)
		zpS.setMethod( ZProjector.MAX_METHOD )
		zpS.doProjection()
		MIPimpS = zpS.getProjection()
		MIPimpS.show()