function tp = applyTransformations( p, planePos, u, v, TF, dataStr )
tp = projectOnPlane( p, planePos, u, v );
tp = transformPoint3d( tp, TF );
rotMat = getManualRotationMatrix( dataStr );
tp = [tp 1];
tp = rotMat * tp';
tp = [ tp(1) tp(2) tp(3) ];

% if strcmp( dataStr, '130607_raw' )
%   tp(1) = -tp(1);
% end