function [ u, v, dir, planePos, TF ] = initTransformations( dataStr, cView, registerBase )
% get stored eigenvectors for the last time step to set the same
% direction view for each time step
coeff = getNormalizedPrincipalComponents( dataStr, 1 );

% set PC depending on the viewing direction
if cView == 1
  dir = coeff(:,2);
  u = coeff(:,1);
  v = coeff(:,3);
elseif cView == 2
  dir = coeff(:,3);
  u = coeff(:,1);
  v = coeff(:,2);
  if strcmp( dataStr, '121211_raw' )
    v = -v;
  end
elseif cView == 3
  dir = -coeff(:,1);
  u = -coeff(:,3);
  v = coeff(:,2);
  if strcmp( dataStr, '121211_raw' )
    v = -v;
  end
end

% set plane position
planePos = dir * 1;

plane = [ planePos(1) planePos(2) planePos(3)...
  u(1) u(2) u(3)...
  v(1) v(2) v(3) ];
TF = createBasisTransform3d( 'g', plane );

if strcmp( dataStr, '121204_raw_2014' )
  TF = createRotationOz( degtorad( 1. ) ) * TF;
elseif strcmp( dataStr, '121211_raw' )
  TF = createRotationOz( degtorad( 1. ) ) * TF;
elseif strcmp( dataStr, '130607_raw' )
  TF = createRotationOz( degtorad( -3.5 ) ) * TF;
end

% apply manual translation of data sets for registration
TF = getDataSetTranslationMatrix( dataStr, registerBase ) * TF;

