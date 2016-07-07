function generateCellShape( width, height, slices, cc, maxInt, newFileName )
% set all values to zero
imageStack = zeros( height, width, slices, 'uint16' );

% 0 for sphere or 1 for ellipsoid
cellRadius = 15;
shape = 1;

% generate a new cell at the cell center
cellC = floor(cc);
% ellipsoid parameter
a = cellRadius;
if shape == 1
  b = floor(cellRadius/2);
  d = floor(cellRadius/2);
else
  b = a;
  d = a;
end
for c=1:size(cc,1)
  maxIntensity = maxInt(c,1);
  minIntensity = maxInt(c,1)/2.;
  xStart = cellC(c,1)-a;
  xEnd = cellC(c,1)+a;
  yStart = cellC(c,2)-b;
  yEnd = cellC(c,2)+b;
  zStart = cellC(c,3)-d;
  zEnd = cellC(c,3)+d;
  for x=xStart:xEnd
    for y=yStart:yEnd
      for z=zStart:zEnd
        xx = (x - cellC(c,1))*(x - cellC(c,1));
        yy = (y - cellC(c,2))*(y - cellC(c,2));
        zz = (z - cellC(c,3))*(z - cellC(c,3));
        %if sqrt(xx + yy + zz) <= cellRadius
        length = xx/(a*a) + yy/(b*b) + zz/(d*d);
        % length is in [0,1]
        % 0 -> max intensity
        % 1 -> min intensity
        intensity = length * minIntensity + (1-length) * maxIntensity;
        if xx/(a*a) + yy/(b*b) + zz/(d*d) <= 1.
          if x > 0 && x < width && y > 0 && y < height && z > 0 && z < slices
            imageStack( y, x, z ) = intensity;
          end
        end
      end
    end
  end
end
writeTIFstack( imageStack, char(newFileName), 2^31 );
