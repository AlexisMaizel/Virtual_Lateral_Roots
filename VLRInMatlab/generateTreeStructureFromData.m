function generateTreeStructureFromData( cellData, numTimeSteps, startT, endT )
forestMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
CellIDToNodeIDMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
cIDtoNIDMap = cell( numTimeSteps, 1 );
nodeIndexListOfNewRoots = [];
for t=1:numTimeSteps
  cIDtoNIDMap{t,1} = CellIDToNodeIDMap;
end
dimSize = size( cellData, 1 );
for i=1:dimSize
  cellID = cellData{i, 1};
  time = cellData{i, 5};
  precursorID = cellData{i, 8};
  lin = cellData{i, 6};
  
  % check if the current time step is later than endT
  if time > endT
    continue;
  end
  
  % generate new lineage tree with root node
  if precursorID == -1
    tr = tree( [ cellData{ i, : } ] );
    forestMap(lin) = tr;
    cIDtoNIDMap{time, 1}(cellID) = 1;
    % else add node to existing lineage
  else
    tr = forestMap(lin);
    addToNodeID = cIDtoNIDMap{time, 1}(precursorID);
    [ tr, newNodeID ] = tr.addnode( addToNodeID, [ cellData{ i, : } ] );
    forestMap(lin) = tr;
    cIDtoNIDMap{time, 1}(cellID) = newNodeID;
  end
  
  % if the start time is greater than 1 we store the index nodeID of the
  % new root nodes for which later the set of new subtrees are generated
  if time == startT && startT > 1
    if precursorID == -1
      nodeIndexListOfNewRoots = [ nodeIndexListOfNewRoots ; lin 1 ];
    else
      nodeIndexListOfNewRoots = [ nodeIndexListOfNewRoots ; lin newNodeID ];
    end
  end
end

numLineages = 1;
allNodes = [];
if startT > 1
  croppedForestMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
  for l=1:size( nodeIndexListOfNewRoots, 1)
    lin = nodeIndexListOfNewRoots( l, 1 );
    nodeID = nodeIndexListOfNewRoots( l, 2 );
    croppedForestMap(numLineages) = forestMap(lin).subtree(nodeID);
    if l ~= size( nodeIndexListOfNewRoots, 1)
      numLineages = numLineages + 1;
    end
  end
  allKeys = cell2mat( keys(croppedForestMap) );
  for k=1:size( allKeys, 2 )
    allNodes = [ allNodes croppedForestMap( allKeys(1, k) ).nnodes() ];
  end
else
  allKeys = cell2mat( keys(forestMap) );
  numLineages = size( allKeys, 2 );
  
  for k=1:size( allKeys, 2 )
    allNodes = [ allNodes forestMap( allKeys(1, k) ).nnodes() ];
  end
end
disp( strcat( 'Lineages:', num2str(numLineages) ) );
