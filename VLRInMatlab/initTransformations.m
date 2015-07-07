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
  u = coeff(:,3);
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
  TF = createRotationOz( degTorad( 2. ) ) * TF;
elseif strcmp( dataStr, '121211_raw' )
  TF = createRotationOz( degTorad( 1. ) ) * TF;
elseif strcmp( dataStr, '130508_raw' )
  TF = createRotationOz( degTorad( -1.5 ) ) * TF;
elseif strcmp( dataStr, '130607_raw' )
  TF = createRotationOz( degTorad( -4. ) ) * TF;
end

if cView == 3
  if strcmp( dataStr, '120830_raw' )
    TF = createRotationOz( degTorad( 10. ) ) * TF;
    TF = createTranslation3d( [ 10. 0 0 ] ) * TF;
  elseif strcmp( dataStr, '121204_raw_2014' )
    TF = createRotationOz( degTorad( -10. ) ) * TF;
    TF = createTranslation3d( [ -10 0 0 ] ) * TF;
  elseif strcmp( dataStr, '130508_raw' )
    TF = createTranslation3d( [ -25 0 0 ] ) * TF;
  end
end

if cView == 1
  if strcmp( dataStr, '121204_raw_2014' )
    TF = createRotationOz( degTorad( -3. ) ) * TF;
    TF = createTranslation3d( [ 0 -9 0 ] ) * TF;
  elseif strcmp( dataStr, '130508_raw' )
    TF = createRotationOz( degTorad( 2. ) ) * TF;
    TF = createTranslation3d( [ 0 -35 0 ] ) * TF;
  elseif strcmp( dataStr, '130607_raw' )
    TF = createRotationOz( degTorad( -5. ) ) * TF;
  end
end

% apply manual translation of data sets for registration
TF = getDataSetTranslationMatrix( dataStr, registerBase ) * TF;

