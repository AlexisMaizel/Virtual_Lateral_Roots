function exportBezierSurface( curI, surface, surfaceDir, dimCP )
if curI < 10
  digit = strcat( '_00' );
elseif curI < 100
  digit = strcat( '_0' );
else
  digit = strcat( '_' );
end

% output file for complete bezier surface at curI
fileName = strcat( surfaceDir, 'bezierTensorSurface', digit, num2str(curI), '.txt' );
fileId = fopen( char(fileName), 'w' );

% loop over control points
for i=1:dimCP
  for j=1:dimCP
    fprintf( fileId, '%4f %4f %4f ', surface( i, j, : ) );
  end
  fprintf( fileId, '\n' );
end
fprintf( fileId, '\n' );
