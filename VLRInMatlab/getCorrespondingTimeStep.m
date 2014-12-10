function [ curT, numNormCells ] = getCorrespondingTimeStep( curI, minI, maxI, maxT, numCellsPerTimeStep )

% only continue with this time step that is synchronized in
% number of cells for each data set
numNormCells = getNormalizedCellNumber( curI, 18, 143, minI, maxI );
curT = 1;
epsilon = 1;
found = 0;
while found == 0
  for j=1:maxT
    if numNormCells - epsilon < numCellsPerTimeStep(j,1)...
        && numNormCells + epsilon > numCellsPerTimeStep(j,1)
      curT = j;
      found = 1;
      break;
    end
  end
  % if no time step was found with the corresponding number of cells
  % then increase the search radius by setting epsilon
  if found == 0
    epsilon = epsilon + 1;
  end
end