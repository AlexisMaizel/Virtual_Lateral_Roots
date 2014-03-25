function tp = applyTransformations( p, planePos, u, v, TF, dataStr, renderMasterFile, index )
tp = projectOnPlane( p, planePos, u, v );
tp = transformPoint3d( tp, TF );
% manually choose always the transformation applied to the last index aka
% number of cells = 143
%index = 5;
rotMat = getManualRotationMatrix( dataStr, renderMasterFile, index );
transMat = getManualTranslationMatrix( dataStr, renderMasterFile, index );
tp = [tp 1];
tp = transMat * tp';
tp = rotMat * tp;
tp = [ tp(1) tp(2) tp(3) ];