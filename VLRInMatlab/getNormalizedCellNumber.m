function numCells = getNormalizedCellNumber( normIndex, minCells, maxCells, minIndex, maxIndex )

ratio = (maxCells-minCells)/(maxIndex-minIndex);

numCells = normIndex * ratio;
numCells = numCells + minCells - ratio;
numCells = round(numCells);