function [Q, curS] = updateBezierSurfaceBoundary( factor, initialS, finalS,...
  curS, dimCP, steps )
% loop over control points
for i=1:dimCP
  for j=1:dimCP
    % only interpolate boundary control points
    if i == 1 || i == dimCP || j == 1 || j == dimCP
      curS( i, j, : ) = (1-factor) * initialS( i, j, : ) +...
        factor * finalS( i, j, : );
    end
  end
end
% resolution of interpolated points in bezier surface; actually
% only required for plotting the bezier surface
Q = bezierSurfaceInterp( curS, steps );