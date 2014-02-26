function indexSet = determineColorAssignment( eigenVectors, positiveIndex, principalComponents )

% index set for the two principal components of the 2D ellipse
% the number indicates the number of 'similar' vectors related to the three
% principal components for the 3D ellipsoids. The similarity is measures by
% the cosine angle between the vectors.
indexSet = [ 0 0 ];

for i=1:2
  minAngle = 180;
  index = 1;
  for a=1:3
    ei = eigenVectors( a, : );
    pc = principalComponents( i, : );
    curAngle = acos( dot(ei,pc)/(norm(ei)*norm(pc)) );
    curAngle = curAngle * 180/pi;
    
    if curAngle < minAngle
      minAngle = curAngle;
      index = a;
    end
  end
  % set the positive or negative eigenvalue information of the closest
  % eigenvector
  indexSet( 1, i ) = positiveIndex( index, 1 );
end
