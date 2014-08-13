function points = generateContourPoints( dataName, first )
numContourMarks = 17;
points = zeros( numContourMarks, 3 );
if strcmp( dataName, '120830_raw' )
  if first == true
    points( 1, 1:2 ) = [ 112 -68 ];
    points( 2, 1:2 ) = [ 160 -65 ];
    points( 3, 1:2 ) = [ 230 -59 ];
    points( 4, 1:2 ) = [ 300 -65 ];
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
end