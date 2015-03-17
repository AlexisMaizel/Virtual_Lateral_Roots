function Q = bezierSurfaceInterp( S, steps )

u = linspace(0,1,steps);
v = u;

for i = 1:length(u)
  for j = 1:length(v)
    Q( i, j, : ) = computeBezierCoordinate( S, u(i), v(j) );
  end
end
