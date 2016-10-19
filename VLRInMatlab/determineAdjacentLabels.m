function labels = determineAdjacentLabels( WSImage, subs )
height = size( WSImage, 1 );
width = size( WSImage, 2 );
slices = size( WSImage, 3 );
labels = [];
if subs(1)+1 <= height
  labels = [ labels, WSImage( subs(1)+1, subs(2), subs(3) ) ];
end
if subs(1)-1 > height
  labels = [ labels, WSImage( subs(1)-1, subs(2), subs(3) ) ];
end
if subs(2)+1 <= width
  labels = [ labels, WSImage( subs(1), subs(2)+1, subs(3) ) ];
end
if subs(2)-1 > width
  labels = [ labels, WSImage( subs(1), subs(2)-1, subs(3) ) ];
end
if subs(3)+1 <= slices
  labels = [ labels, WSImage( subs(1), subs(2), subs(3)+1 ) ];
end
if subs(3)-1 > slices
  labels = [ labels, WSImage( subs(1), subs(2), subs(3)-1 ) ];
end
labels = unique( labels );

