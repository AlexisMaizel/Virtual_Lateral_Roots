function generateCellShape( imageStack, cc, maxInt, cellRadius, newFileName )
% set all values to zero
imageStack(:,:,:) = 0;
height = size( imageStack, 1);
width = size( imageStack, 2);
slices = size( imageStack, 3);

% generate a new cell at the cell center
cellC = floor(cc);
for c=1:size(cc,1)
  xStart = cellC(c,1)-cellRadius;
  xEnd = cellC(c,1)+cellRadius;
  yStart = cellC(c,2)-cellRadius;
  yEnd = cellC(c,2)+cellRadius;
  zStart = cellC(c,3)-cellRadius;
  zEnd = cellC(c,3)+cellRadius;
  for x=xStart:xEnd
    for y=yStart:yEnd
      for z=zStart:zEnd
        if x > 0 && x < width && y > 0 && y < height && z > 0 && z < slices
          imageStack( y, x, z ) = maxInt(c,1);
        end
      end
    end
  end
end
writeTIFstack( imageStack, char(newFileName), 2^31 );