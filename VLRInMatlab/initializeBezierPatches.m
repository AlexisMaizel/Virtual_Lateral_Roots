function [ S, finalS, Q ] = initializeBezierPatches( numPatches, min, max,...
  offset, emphasizeDomeTip )

% each patch has 16 control points
%number of interpolated values between end control points
ni = 10;
%uniform parameterization
u = linspace(0,1,ni);
v = u;

% generate matrix of control points
if numPatches == 1
  stepSize = zeros( 2 );
  for s=1:2
    stepSize(1, s) = abs(2*offset + max(1, s) - min(1, s))/3.;
  end
 
  for i=1:4
    for j=1:4
      % x coord
      S( i, j, 1, 1 ) = min(1, 1) + (j-1)*stepSize(1,1) - offset;
      % y coord
      S( i, j, 2, 1 ) = min(1, 2) + (i-1)*stepSize(1,2) - offset;
      % z coord
      S( i, j, 3, 1 ) = 0.;
    end
  end
  
  % interpolation of patch
  Q(:,:,:,1) = bezierpatchinterp(S(:,:,:,1), u, v);
  
elseif numPatches == 4
  stepSize = zeros(2);
  half = zeros(2);
  for s=1:2
    stepSize(1, s) = abs(2*offset + max(1, s) - min(1, s))/6.;
    half(1, s) = abs(max(1, s) - min(1, s))/2.;
  end
  % loop over patches
  for pi=1:2
    for pj=1:2
      % loop over control points
      pIndex = pj + (pi-1)*2;
      for i=1:4
        for j=1:4
          % x coord
          S( i, j, 1, pIndex ) = (pj-1)*half(1,1) + min(1, 1) + (j-1)*stepSize(1,1) - offset;
          % y coord
          S( i, j, 2, pIndex ) = (pi-1)*half(1,2) + min(1, 2) + (i-1)*stepSize(1,2) - offset;
          % z coord
          S( i, j, 3, pIndex ) = 0.;
        end
      end
      % interpolation of patch
      Q(:,:,:,pIndex) = bezierpatchinterp(S(:,:,:,pIndex), u, v);
    end
  end
  
  % manually initialize the bezier surface at the final registered time step
  finalS = S;
  % loop over patches
  for p=1:numPatches
    % bottom left
    if p == 1
      finalS( 1, 1, :, p ) = [ -280 -50 0 ];
      finalS( 1, 2, :, p ) = [ -180 -50 0 ];
      finalS( 1, 3, :, p ) = [ -80 -50 0 ];
      finalS( 1, 4, :, p ) = [ 0 -50 0 ];
      finalS( 2, 1, :, p ) = [ -280 -45 0 ];
      finalS( 3, 1, :, p ) = [ -280 -40 0 ];
      finalS( 4, 1, :, p ) = [ -280 -35 0 ];
    % bottom right
    elseif p == 2
      finalS( 1, 1, :, p ) = [ 0 -50 0 ];
      finalS( 1, 2, :, p ) = [ 100 -50 0 ];
      finalS( 1, 3, :, p ) = [ 200 -50 0 ];
      finalS( 1, 4, :, p ) = [ 300 -50 0 ];
      finalS( 2, 4, :, p ) = [ 300 -45 0 ];
      finalS( 3, 4, :, p ) = [ 300 -40 0 ];
      finalS( 4, 4, :, p ) = [ 300 -35 0 ];
    % top left
    elseif p == 3
      if emphasizeDomeTip == 0
        finalS( 4, 1, :, p ) = [ -280 -20 0 ];
        finalS( 4, 2, :, p ) = [ -180 15 0 ];
        finalS( 4, 3, :, p ) = [ -80 90 0 ];
        finalS( 4, 4, :, p ) = [ 0 85 0 ];
      else
        finalS( 4, 1, :, p ) = [ -280 -10 0 ];
        finalS( 4, 2, :, p ) = [ -180 30 0 ];
        finalS( 4, 3, :, p ) = [ -80 180 0 ];
        finalS( 4, 4, :, p ) = [ 0 170 0 ];
      end
      finalS( 1, 1, :, p ) = [ -280 -35 0 ];
      finalS( 2, 1, :, p ) = [ -280 -30 0 ];
      finalS( 3, 1, :, p ) = [ -280 -25 0 ];
    % top right
    else
      if emphasizeDomeTip == 0
        finalS( 4, 1, :, p ) = [ 0 85 0 ];
        finalS( 4, 2, :, p ) = [ 100 90 0 ];
        finalS( 4, 3, :, p ) = [ 200 15 0 ];
        finalS( 4, 4, :, p ) = [ 300 -20 0 ];
      else
        finalS( 4, 1, :, p ) = [ 0 170 0 ];
        finalS( 4, 2, :, p ) = [ 100 180 0 ];
        finalS( 4, 3, :, p ) = [ 200 30 0 ];
        finalS( 4, 4, :, p ) = [ 300 -10 0 ];
      end
      finalS( 1, 4, :, p ) = [ 300 -35 0 ];
      finalS( 2, 4, :, p ) = [ 300 -30 0 ];
      finalS( 3, 4, :, p ) = [ 300 -25 0 ];
    end
    % interpolation of patch
    %Q(:,:,:,p) = bezierpatchinterp(finalS(:,:,:,p), u, v);
  end
else
  disp( strcat( 'Number of patches is not supported yet! ', num2str(numPatches) ) );
end
