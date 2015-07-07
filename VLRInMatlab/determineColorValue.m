function value = determineColorValue( magnitudes, min, max )

average = 0;
for m=1:size( magnitudes, 1 )
  average = average + magnitudes(m);
end

average = average / size( magnitudes, 1 );
average = convertToReal( log(average) );
value = ( average - min )/( max - min );