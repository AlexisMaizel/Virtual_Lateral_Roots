function points = generateContourPoints( dataName, first, eps )
% I always choose 16 contour points for which 7 are used
% at the top and bottom while 1 is used for left and right
% 17 are chosen because the first one has to be set again
% to enclose the surface
numContourMarks = 17;
points = zeros( numContourMarks, 3 );
if strcmp( dataName, '120830_raw' )
  if first == true
    points( 1, 1:2 ) = [ 112 -68 ];
    points( 2, 1:2 ) = [ 160 -65 ];
    points( 3, 1:2 ) = [ 230 -65 ];
    points( 4, 1:2 ) = [ 300 -59 ];
    points( 5, 1:2 ) = [ 370 -65 ];
    points( 6, 1:2 ) = [ 440 -65 ];
    points( 7, 1:2 ) = [ 520 -65 ];
    points( 8, 1:2 ) = [ 520 -75 ];
    points( 9, 1:2 ) = [ 520 -100 ];
    points( 10, 1:2 ) = [ 440 -100 ];
    points( 11, 1:2 ) = [ 370 -100 ];
    points( 12, 1:2 ) = [ 300 -100 ];
    points( 13, 1:2 ) = [ 230 -100 ];
    points( 14, 1:2 ) = [ 160 -100 ];
    points( 15, 1:2 ) = [ 115 -100 ];
    points( 16, 1:2 ) = [ 113 -75 ];
    points( 17, 1:2 ) = [ 112 -68 ];
  else
    points( 1, 1:2 ) = [ 100 -42 ];
    points( 2, 1:2 ) = [ 160 -15 ];
    points( 3, 1:2 ) = [ 230  18 ];
    points( 4, 1:2 ) = [ 275  45 ];
    points( 5, 1:2 ) = [ 370  10 ];
    points( 6, 1:2 ) = [ 440 -20 ];
    points( 7, 1:2 ) = [ 520 -52 ];
    points( 8, 1:2 ) = [ 520 -65 ];
    points( 9, 1:2 ) = [ 520 -100 ];
    points( 10, 1:2 ) = [ 440 -100 ];
    points( 11, 1:2 ) = [ 370 -100 ];
    points( 12, 1:2 ) = [ 300 -100 ];
    points( 13, 1:2 ) = [ 230 -100 ];
    points( 14, 1:2 ) = [ 160 -100 ];
    points( 15, 1:2 ) = [ 115 -100 ];
    points( 16, 1:2 ) = [ 113 -75 ];
    points( 17, 1:2 ) = [ 100 -42 ];  
  end
elseif strcmp( dataName, '121211_raw' )
  if first == true
    points( 1, 1:2 ) = [ -10 -45 ];
    points( 2, 1:2 ) = [ 90 -35 ];
    points( 3, 1:2 ) = [ 190 -35 ];
    points( 4, 1:2 ) = [ 290 -35 ];
    points( 5, 1:2 ) = [ 390 -35 ];
    points( 6, 1:2 ) = [ 490 -50 ];
    points( 7, 1:2 ) = [ 600 -60 ];
    points( 8, 1:2 ) = [ 600 -75 ];
    points( 9, 1:2 ) = [ 600 -90 ];
    points( 10, 1:2 ) = [ 490 -90 ];
    points( 11, 1:2 ) = [ 390 -90 ];
    points( 12, 1:2 ) = [ 290 -90 ];
    points( 13, 1:2 ) = [ 190 -90 ];
    points( 14, 1:2 ) = [ 90 -90 ];
    points( 15, 1:2 ) = [ -10 -90 ];
    points( 16, 1:2 ) = [ -10 -70 ];
    points( 17, 1:2 ) = [ -10 -35 ];
  else
    points( 1, 1:2 ) = [ -10 -45 ];
    points( 2, 1:2 ) = [ 90 -25 ];
    points( 3, 1:2 ) = [ 190 60 ];
    points( 4, 1:2 ) = [ 290 130 ];
    points( 5, 1:2 ) = [ 390 25 ];
    points( 6, 1:2 ) = [ 490 -45 ];
    points( 7, 1:2 ) = [ 600 -60 ];
    points( 8, 1:2 ) = [ 600 -75 ];
    points( 9, 1:2 ) = [ 600 -90 ];
    points( 10, 1:2 ) = [ 490 -90 ];
    points( 11, 1:2 ) = [ 390 -90 ];
    points( 12, 1:2 ) = [ 290 -90 ];
    points( 13, 1:2 ) = [ 190 -90 ];
    points( 14, 1:2 ) = [ 90 -90 ];
    points( 15, 1:2 ) = [ -10 -90 ];
    points( 16, 1:2 ) = [ -10 -70 ];
    points( 17, 1:2 ) = [ -10 -45 ];  
  end
end

% export of initial contour positions used in the model before the
% offset is applied such that the model contour points are always
% within the increased surface realized with the eps offset
if first == true
  fileName = strcat( '/tmp/conPoints-', dataName, '.txt' );
  fileId = fopen( char(fileName), 'w' );
  fprintf( fileId, '%1d\n', numContourMarks-1 );
  for p=1:size(points,1)-1
    fprintf( fileId, '%4f %4f\n', points(p,1:2) );
  end
  fprintf( fileId, '\n' );
  fclose( fileId );
end

% apply the offset to the points according to their position
% on the boundary of the surface
% top left point
points(1,1:2) = [ points(1,1)-eps points(1,2)+eps ];
for p=2:6
  points(p,1:2) = [ points(p,1) points(p,2)+eps ];
end
% top right point
points(7,1:2) = [ points(7,1)+eps points(7,2)+eps ];
points(8,1:2) = [ points(8,1)+eps points(8,2) ];
% bottom right point
points(9,1:2) = [ points(9,1)+eps points(9,2)-eps ];
for p=10:14
  points(p,1:2) = [ points(p,1) points(p,2)-eps ];
end
% bottom left point
points(15,1:2) = [ points(15,1)-eps points(15,2)-eps ];
points(16,1:2) = [ points(16,1)-eps points(16,2) ];
points(17,1:2) = [ points(17,1)-eps points(17,2)+eps ];
