function pos = computeBezierCoordinate( S, u, v )

pos = zeros( 3, 1 );
dim = size( S, 1 );

for i=0:dim-1
  for j=0:dim-1
    s = nchoosek(dim-1,i) * power(v,i) * power(1-v, dim-1-i) *...
      nchoosek(dim-1,j) * power(u,j) * power(1-u, dim-1-j);
    % x coord
    pos(1) = pos(1) + (s * S( i+1, j+1, 1 ));
    % y coord
    pos(2) = pos(2) + (s * S( i+1, j+1, 2 ));
  end
end