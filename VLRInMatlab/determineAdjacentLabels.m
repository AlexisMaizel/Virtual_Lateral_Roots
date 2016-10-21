function labels = determineAdjacentLabels( WSImage, subs )
height = size( WSImage, 1 );
width = size( WSImage, 2 );
slices = size( WSImage, 3 );
labels = [];
nRadius = 2;

% check boundaries: 3D case
if slices > 1
  if subs(1)-nRadius < 1 || subs(2)-nRadius < 1 || subs(3)-nRadius < 1 ||...
      subs(1)+nRadius > height || subs(2)+nRadius > width || subs(3)+nRadius > slices
    return;
  end
  % create sphere/cube
  [ xs, ys, zs ] = ndgrid( -nRadius:nRadius, -nRadius:nRadius, -nRadius:nRadius );
  sphereVoxels = (xs).^2 + (ys).^2 + (zs).^2 <= nRadius.^2;
  
  % index of all non-zero elements and we only have to add the start1 and
  % start2 values to get all spherical values with center pos1 and pos2
  ind = find(sphereVoxels);
  [ xi, yi, zi ] = ind2sub( size(sphereVoxels), ind );
  
  si = size( xi, 1 );
  labels = zeros( 1, si );
  for i=1:si
    index = [ xi( i, 1 )-nRadius-1+subs(1) yi( i, 1 )-nRadius-1+subs(2) zi( i, 1 )-nRadius-1+subs(3) ];
    labels( 1, i ) = WSImage( index(1), index(2), index(3) );
  end
  % 2D case
else
  if subs(1)-nRadius < 1 || subs(2)-nRadius < 1 || subs(1)+nRadius > height ||...
      subs(2)+nRadius > width
    return;
  end
  % create circle/square
  [ xs, ys ] = ndgrid( -nRadius:nRadius, -nRadius:nRadius );
  circleVoxels = (xs).^2 + (ys).^2 <= nRadius.^2;
  
  % index of all non-zero elements and we only have to add the start1 and
  % start2 values to get all spherical values with center pos1 and pos2
  ind = find(circleVoxels);
  [ xi, yi ] = ind2sub( size(circleVoxels), ind );
  
  si = size( xi, 1 );
  labels = zeros( 1, si );
  for i=1:si
    index = [ xi( i, 1 )-nRadius-1+subs(1) yi( i, 1 )-nRadius-1+subs(2) ];
    labels( 1, i ) = WSImage( index(1), index(2), 1 );
  end
end
labels = unique( labels );
