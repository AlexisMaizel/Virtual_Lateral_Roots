function [Q, curS] = updateBezierSurfaceBoundary( factor, initialS, finalS, curS, numPatches )
%number of interpolated values between end control points
ni = 10;
%uniform parameterization
u = linspace(0,1,ni);
v = u;

if numPatches == 1
  for i=1:4
    for j=1:4
      % only interpolate boundary control points
      if i == 1 || i == 4 || j == 1 || j == 4
        curS( i, j, :, 1 ) = (1-factor) * initialS( i, j, :, 1 ) + factor * finalS( i, j, :, 1 );
      end
    end
  end
  
  % interpolation of patch
  Q(:,:,:,1) = bezierpatchinterp(curS(:,:,:,1), u, v);
elseif numPatches == 4
  % loop over patches
  for p=1:numPatches
    % loop over control points
    for i=1:4
      for j=1:4
        % only interpolate boundary control points
        if (i == 1 && (p == 1 || p == 2)) ||...
            (i == 4 && (p == 3 || p == 4)) ||...
            (j == 1 && (p == 1 || p == 3)) ||...
            (j == 4 && (p == 2 || p == 4))
          curS( i, j, :, p ) = (1-factor) * initialS( i, j, :, p ) + factor * finalS( i, j, :, p );
        end
      end
    end
    % interpolation of patch
    Q(:,:,:,p) = bezierpatchinterp(curS(:,:,:,p), u, v);
  end
else
  disp( strcat( 'Number of patches is not supported yet! ', num2str(numPatches) ) );
end