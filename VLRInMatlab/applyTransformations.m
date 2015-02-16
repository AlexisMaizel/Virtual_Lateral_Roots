function tp = applyTransformations( p, planePos, u, v, TF, dataStr, renderMasterFile )
tp = projectOnPlane( p, planePos, u, v );
tp = transformPoint3d( tp, TF );
% manually choose always the transformation applied to the last index aka
% number of cells = 143
%index = 5;
%renderMasterFile = 0;
rotMat = getManualRotationMatrix( dataStr, renderMasterFile );
transMat = getManualTranslationMatrix( dataStr, renderMasterFile );
tp = [tp 1];
tp = transMat * tp';
tp = rotMat * tp;
tp = [ tp(1) tp(2) tp(3) ];