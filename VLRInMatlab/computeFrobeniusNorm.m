function norm = computeFrobeniusNorm( matrix )
nRows = size( matrix, 1 );
nColumns = size( matrix, 2 );
sum = 0.;

for i=1:nRows
  for j=1:nColumns
    sum = sum + matrix( i, j ) * matrix( i, j );
  end
end

norm = sqrt( sum );