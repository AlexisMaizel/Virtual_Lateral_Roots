function [ S, finalS, Q ] = initializeBezierSurface( dimCP, min, max,...
  offset, steps, emphasizeDomeTip )
% initial bezier surface
stepSize = zeros(2);
half = zeros(2);
for s=1:2
  stepSize(1, s) = abs(2*offset + max(1, s) - min(1, s))/(dimCP-1);
  half(1, s) = abs(max(1, s) - min(1, s))/2.;
end
for i=1:dimCP
  for j=1:dimCP
    % x coord
    S( i, j, 1 ) = min(1, 1) + (j-1)*stepSize(1,1) - offset;
    % y coord
    S( i, j, 2 ) = min(1, 2) + (i-1)*stepSize(1,2) - offset;
    % z coord
    S( i, j, 3 ) = 0.;
  end
end

% resolution of interpolated points in bezier surface; actually
% only required for plotting the bezier surface
Q = bezierSurfaceInterp( S, steps );

% final bezier surface is set manually but only for the boundary control
% points while the inner cps are given by the uniform interpolation done
% above
finalS = S;

if dimCP == 7
  % row 1, bottom
  finalS( 1, 1, : ) = [ -280 -50 0 ];
  finalS( 1, 2, : ) = [ -180 -50 0 ];
  finalS( 1, 3, : ) = [ -80 -50 0 ];
  finalS( 1, 4, : ) = [ 0 -50 0 ];
  finalS( 1, 5, : ) = [ 100 -50 0 ];
  finalS( 1, 6, : ) = [ 200 -50 0 ];
  finalS( 1, 7, : ) = [ 300 -50 0 ];
  % row 2 till 6 for column 1 and 7
  finalS( 2, 1, : ) = [ -280 -45 0 ];
  finalS( 3, 1, : ) = [ -280 -40 0 ];
  finalS( 4, 1, : ) = [ -280 -35 0 ];
  finalS( 5, 1, : ) = [ -280 -30 0 ];
  finalS( 6, 1, : ) = [ -280 -25 0 ];
  finalS( 2, 7, : ) = [ 300 -45 0 ];
  finalS( 3, 7, : ) = [ 300 -40 0 ];
  finalS( 4, 7, : ) = [ 300 -35 0 ];
  finalS( 5, 7, : ) = [ 300 -30 0 ];
  finalS( 6, 7, : ) = [ 300 -25 0 ];
  % row 7, top
  if emphasizeDomeTip == 0
    finalS( 7, 1, : ) = [ -280 -20 0 ];
    finalS( 7, 2, : ) = [ -180 15 0 ];
    finalS( 7, 3, : ) = [ -80 90 0 ];
    finalS( 7, 4, : ) = [ 0 85 0 ];
    finalS( 7, 5, : ) = [ 100 90 0 ];
    finalS( 7, 6, : ) = [ 200 15 0 ];
    finalS( 7, 7, : ) = [ 300 -20 0 ];
  else
    finalS( 7, 1, : ) = [ -280 -10 0 ];
    finalS( 7, 2, : ) = [ -180 30 0 ];
    finalS( 7, 3, : ) = [ -80 180 0 ];
    finalS( 7, 4, : ) = [ 0 170 0 ];
    finalS( 7, 5, : ) = [ 100 180 0 ];
    finalS( 7, 6, : ) = [ 200 30 0 ];
    finalS( 7, 7, : ) = [ 300 -10 0 ];
  end
else
  disp( strcat( 'Manual final surface is not supported yet for current dimCP ', num2str(dimCP) ) );
end