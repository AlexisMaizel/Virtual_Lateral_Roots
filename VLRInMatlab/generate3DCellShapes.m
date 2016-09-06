function [ h, cellTracks ] = generate3DCellShapes( imageStack, txtPath, slices, first, t, cellTracks )
connectivity = 26;
height = size( imageStack, 1 );
CC = bwconncomp( imageStack, connectivity );
S = regionprops( CC, 'Centroid', 'Area', 'BoundingBox', 'PixelList', 'PixelIdxList' );
numCCs = size(S, 1);
cm = colormap( jet(numCCs) );
remainingCCs = 0;
hold on
shading interp
light
lighting phong
fileID = fopen( char(txtPath), 'w' );
cellShapes = [];
for i=1:numCCs
  disp( strcat( num2str(i), '/', num2str(numCCs) ) );
  voxels = S(i, :).PixelList;
  numVoxels = size( voxels, 1 );
  if numVoxels > 50 && numVoxels < 500000
    %k = boundary( voxels );
    k = convhull( voxels );
    if slices == 1
      h = plot( voxels(k,1), -voxels(k,2)+height );
    else
      color = cm( i, : );
      h = trimesh( k, voxels(:,1), -voxels(:,2)+height, voxels(:,3),...
        'Facecolor', color, 'FaceLighting', 'gouraud', 'LineStyle', 'none', 'FaceAlpha', 1 );
      [ rf, rv ] = reducepatch( h, 0.05 );
      fprintf( fileID,'%u\n', size( rv, 1 ) );
      for v=1:size(rv, 1)
        fprintf( fileID,'%f %f %f\n', rv(v, 1), rv(v, 2), rv(v, 3) );
      end
    end
    remainingCCs = remainingCCs + 1;
    cellShapes = [ cellShapes ; { S( i, : ).Centroid, remainingCCs, cm( remainingCCs, : ) } ];
  end
end
% if it is the first time step store for each remaining cell its
% centroid and its unique color for identifcation of cell tracks
if first == 1
  cellTracks( t ) = cellShapes;
else
  % else consider the initial set of values from the last time step and
  % find the nearest neighbor of the current extracted cell shape and
  % store it in the cellTracks map
  valueSet = values( cellTracks( t-1 ) );
  % find for each cell shape in the current time step the nearest neighbor
  % in the last time step and create a cell track between them
  rows = size( cellShapes, 1 );
  for r=1:rows
    % TODO
  end
end
fclose(fileID);
remainingCCs
hold off