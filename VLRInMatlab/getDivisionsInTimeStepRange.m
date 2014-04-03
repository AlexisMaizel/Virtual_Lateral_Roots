function divs = getDivisionsInTimeStepRange( startT, endT, onlyMasterCellFile, cellFileMap, divisionProperties )
numDivisions = size( divisionProperties, 1 );
divs = [];
for d=1:numDivisions
  % first check time step range
  time = divisionProperties( d, 2 );
  if time >= startT && time <= endT
    % then check if cell division occurs in the master cell file
    if onlyMasterCellFile == 1 &&...
      0 ~= cellFileMap( divisionProperties( d, 1 ) )
      continue;
    end
    divs = [ divs; divisionProperties(d, :) ];
  end
end