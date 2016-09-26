function [nI] = normalizeImage( I, bd )
im = im2double( I );
nI = ( im - mean( im(:) ) )/ std( im(:) );
if bd == 8
  nI = im2uint8( nI );
elseif bd == 16
  nI = im2uint16( nI );
end