function index = getNeighborTileIndex( neighbor, centerIndex, columns )

% start at top left -> first row
if neighbor > 0 && neighbor < 4
  index = centerIndex - columns + neighbor - 2;
% continue with second row
elseif neighbor == 4
  index = centerIndex - 1;
elseif neighbor == 5
  index = centerIndex + 1;
% and process the last row
else
  index = centerIndex + columns + neighbor - 7;
end