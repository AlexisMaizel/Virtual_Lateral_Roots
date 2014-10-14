function points = generateContourPoints( dataName, first, eps )
% I always choose 16 contour points for which 7 are used
% at the top and bottom while 1 is used for left and right
numContourMarks = 16;
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
  else
    points( 1, 1:2 ) = [ 100 -42 ];
    points( 2, 1:2 ) = [ 160 -15 ];
    points( 3, 1:2 ) = [ 230  18 ];
    points( 4, 1:2 ) = [ 300  65 ];
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
  end
elseif strcmp( dataName, '121204_raw_2014' )
  if first == true
    points( 1, 1:2 ) = [ 112 -220 ];
    points( 2, 1:2 ) = [ 160 -210 ];
    points( 3, 1:2 ) = [ 230 -200 ];
    points( 4, 1:2 ) = [ 300 -195 ];
    points( 5, 1:2 ) = [ 370 -194 ];
    points( 6, 1:2 ) = [ 440 -210 ];
    points( 7, 1:2 ) = [ 520 -220 ];
    points( 8, 1:2 ) = [ 520 -230 ];
    points( 9, 1:2 ) = [ 520 -250 ];
    points( 10, 1:2 ) = [ 440 -245 ];
    points( 11, 1:2 ) = [ 370 -235 ];
    points( 12, 1:2 ) = [ 300 -230 ];
    points( 13, 1:2 ) = [ 230 -230 ];
    points( 14, 1:2 ) = [ 160 -240 ];
    points( 15, 1:2 ) = [ 115 -250 ];
    points( 16, 1:2 ) = [ 113 -230 ];
  else
    points( 1, 1:2 ) = [ 112 -230 ];
    points( 2, 1:2 ) = [ 160 -210 ];
    points( 3, 1:2 ) = [ 230 -175 ];
    points( 4, 1:2 ) = [ 300 -125 ];
    points( 5, 1:2 ) = [ 370 -105 ];
    points( 6, 1:2 ) = [ 470 -160 ];
    points( 7, 1:2 ) = [ 570 -200 ];
    points( 8, 1:2 ) = [ 570 -230 ];
    points( 9, 1:2 ) = [ 570 -245 ];
    points( 10, 1:2 ) = [ 470 -260 ];
    points( 11, 1:2 ) = [ 370 -260 ];
    points( 12, 1:2 ) = [ 300 -260 ];
    points( 13, 1:2 ) = [ 230 -260 ];
    points( 14, 1:2 ) = [ 160 -260 ];
    points( 15, 1:2 ) = [ 115 -260 ];
    points( 16, 1:2 ) = [ 113 -245 ]; 
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
  end
elseif strcmp( dataName, '130508_raw' )
  if first == true
    points( 1, 1:2 ) = [ 25 260 ];
    points( 2, 1:2 ) = [ 130 270 ];
    points( 3, 1:2 ) = [ 180 270 ];
    points( 4, 1:2 ) = [ 230 270 ];
    points( 5, 1:2 ) = [ 280 270 ];
    points( 6, 1:2 ) = [ 330 270 ];
    points( 7, 1:2 ) = [ 400 260 ];
    points( 8, 1:2 ) = [ 400 245 ];
    points( 9, 1:2 ) = [ 400 230 ];
    points( 10, 1:2 ) = [ 330 220 ];
    points( 11, 1:2 ) = [ 280 220 ];
    points( 12, 1:2 ) = [ 230 220 ];
    points( 13, 1:2 ) = [ 180 220 ];
    points( 14, 1:2 ) = [ 130 220 ];
    points( 15, 1:2 ) = [ 25 230 ];
    points( 16, 1:2 ) = [ 25 245 ];
  else
    points( 1, 1:2 ) = [ 25 250 ];
    points( 2, 1:2 ) = [ 130 300 ];
    points( 3, 1:2 ) = [ 180 350 ];
    points( 4, 1:2 ) = [ 230 350 ];
    points( 5, 1:2 ) = [ 280 325 ];
    points( 6, 1:2 ) = [ 330 300 ];
    points( 7, 1:2 ) = [ 400 250 ];
    points( 8, 1:2 ) = [ 400 225 ];
    points( 9, 1:2 ) = [ 400 200 ];
    points( 10, 1:2 ) = [ 330 200 ];
    points( 11, 1:2 ) = [ 280 200 ];
    points( 12, 1:2 ) = [ 230 200 ];
    points( 13, 1:2 ) = [ 180 200 ];
    points( 14, 1:2 ) = [ 130 200 ];
    points( 15, 1:2 ) = [ 25 210 ];
    points( 16, 1:2 ) = [ 25 240 ];
  end    
elseif strcmp( dataName, '130607_raw' )
  if first == true
    points( 1, 1:2 ) = [ -60 240 ];
    points( 2, 1:2 ) = [ 50 250 ];
    points( 3, 1:2 ) = [ 110 250 ];
    points( 4, 1:2 ) = [ 170 250 ];
    points( 5, 1:2 ) = [ 230 250 ];
    points( 6, 1:2 ) = [ 290 250 ];
    points( 7, 1:2 ) = [ 350 240 ];
    points( 8, 1:2 ) = [ 350 215 ];
    points( 9, 1:2 ) = [ 350 175 ];
    points( 10, 1:2 ) = [ 290 175 ];
    points( 11, 1:2 ) = [ 230 185 ];
    points( 12, 1:2 ) = [ 170 195 ];
    points( 13, 1:2 ) = [ 110 195 ];
    points( 14, 1:2 ) = [ 50 195 ];
    points( 15, 1:2 ) = [ -60 195 ];
    points( 16, 1:2 ) = [ -60 215 ];
  else
    points( 1, 1:2 ) = [ -175 225 ];
    points( 2, 1:2 ) = [ -80 235 ];
    points( 3, 1:2 ) = [ 20 265 ];
    points( 4, 1:2 ) = [ 120 350 ];
    points( 5, 1:2 ) = [ 220 370 ];
    points( 6, 1:2 ) = [ 320 275 ];
    points( 7, 1:2 ) = [ 400 175 ];
    points( 8, 1:2 ) = [ 420 165 ];
    points( 9, 1:2 ) = [ 410 120 ];
    points( 10, 1:2 ) = [ 320 150 ];
    points( 11, 1:2 ) = [ 220 160 ];
    points( 12, 1:2 ) = [ 120 180 ];
    points( 13, 1:2 ) = [ 20 180 ];
    points( 14, 1:2 ) = [ -80 175 ];
    points( 15, 1:2 ) = [ -175 175 ];
    points( 16, 1:2 ) = [ -175 175 ]; 
  end
end

% export of initial contour positions used in the model before the
% offset is applied such that the model contour points are always
% within the increased surface realized with the eps offset
if first == true
  fileName = strcat( '/tmp/conPoints-', dataName, '.txt' );
  fileId = fopen( char(fileName), 'w' );
  fprintf( fileId, '%1d\n', numContourMarks );
  for p=1:size(points,1)
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
