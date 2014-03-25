function tp = applyTransformations( p, planePos, u, v, TF, dataStr, renderMasterFile )
tp = projectOnPlane( p, planePos, u, v );
tp = transformPoint3d( tp, TF );
rotMat = getManualRotationMatrix( dataStr, renderMasterFile );
transMat = getManualTranslationMatrix( dataStr, renderMasterFile );
tp = [tp 1];
tp = transMat * tp';
tp = rotMat * tp;
tp = [ tp(1) tp(2) tp(3) ];

% if strcmp( dataStr, '130607_raw' )
%   tp(1) = -tp(1);
% end