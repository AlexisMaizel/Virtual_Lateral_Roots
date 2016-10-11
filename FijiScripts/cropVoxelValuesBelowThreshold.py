# helper methods for segmentation
import ij

# all values below threshold are set to 0, else set to 1
def convertToBinaryMask( imp, threshold ):
	stack = imp.getStack()
	slices = imp.getNSlices()
	for i in range( 0, slices ):
		ip = stack.getProcessor(i+1)
		ip.threshold( threshold )
		stack.setProcessor( ip, i+1 )
	
	return ij.ImagePlus( 'BinaryImage', stack )

# all values below threshold are set to 0, else the value keeps the same
def cropVoxelValuesBelowThreshold( imp, threshold ):
	stack = imp.getStack()
	slices = imp.getNSlices()
	width = imp.getWidth()
	height = imp.getHeight()
	for s in range( 0, slices ):
		ip = stack.getProcessor(s+1)
		pixels = ip.getPixels()
		for i in range( 0, width*height ):
			if pixels[i]&0xffff <= threshold:
				pixels[i] = 0
				
		#ip.setPixels( pixels )
		#stack.setProcessor( ip, s+1 )
	
	return ij.ImagePlus( 'ThresholdImage', stack )

	# stack = imp.getStack()
	# sliceCount = imp.getNSlices()
	# width = imp.getWidth()
	# height = imp.getHeight()
	# print( str(width) + 'x' + str(height) + 'x' + str(sliceCount) )
	# for i in range(0, sliceCount):
		# ip = stack.getProcessor(i+1)
		# for y in range(0, height):
			# for x in range(0, width):
				# value = ip.getPixelValue(x, y)

				# if value < threshold:
					# value = 0
					# ip.putPixelValue( x, y, value )
				# #else:
					# #value = 1
					
		# stack.setProcessor( ip, i+1 )

# #	if imp.getBitDepth() == 8:
	# #	ip2 = ip2.convertToByte(false)
		
	# return ij.ImagePlus( 'NewImage', stack )
