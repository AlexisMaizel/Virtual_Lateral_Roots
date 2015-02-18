function [ curT, numNormCells, numCellExists ] = ...
  getCorrespondingTimeStep( curI, minI, maxI, maxT, numCellsPerTimeStep, numTotalCells, considerAllCells )

% only continue with this time step that is synchronized in
% number of cells for each data set
numRegisteredCellsStart = 18;
if considerAllCells == 1
  numRegisteredCellsEnd = 265;
else
  numRegisteredCellsEnd = 143;
end
numNormCells = getNormalizedCellNumber( curI, numRegisteredCellsStart, numRegisteredCellsEnd, minI, maxI );
curT = 1;
epsilon = 1;
found = 0;
while found == 0
  for j=1:maxT
    if numNormCells - epsilon < numCellsPerTimeStep(j,1)...
        && numNormCells + epsilon > numCellsPerTimeStep(j,1)
      curT = j;
      found = 1;
      numCellExists = 1;
      break;
    end
  end
  % if no time step was found with the corresponding number of cells
  % then increase the search radius by setting epsilon but only if
  % numNormCells is really smaller than the total number of cells in the
  % data set
  if found == 0
    epsilon = epsilon + 1;
  end
end

if numTotalCells < numNormCells && epsilon > 5
  numCellExists = 0;
end