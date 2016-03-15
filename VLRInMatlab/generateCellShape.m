function generateCellShape( imageStack, cc, maxInt, cellRadius, newFileName )
% set all values to zero
height = size( imageStack, 1);
width = size( imageStack, 2);
slices = size( imageStack, 3);
imageStack = zeros( height, width, slices, 'uint16' );

% generate a new cell at the cell center
cellC = floor(cc);
% ellipsoid parameter
a = cellRadius;
b = floor(cellRadius/2);
d = floor(cellRadius/2);
%s = 995;
%t = 1000;
%rng(0, 'twister')
for c=1:size(cc,1)
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
        if xx/(a*a) + yy/(b*b) + zz/(d*d) <= 1.
          if x > 0 && x < width && y > 0 && y < height && z > 0 && z < slices
            %r = (t - s).*rand(1,1) + s;
            %r = floor(r);
            imageStack( y, x, z ) = maxInt(c,1);
          end
        end
      end
    end
  end
end
writeTIFstack( imageStack, char(newFileName), 2^31 );
