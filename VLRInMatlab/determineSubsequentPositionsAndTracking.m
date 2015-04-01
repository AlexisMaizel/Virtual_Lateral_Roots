function [ cPos, nPos ] = determineSubsequentPositionsAndTracking( dimData, curT, nextT,...
  cellDatas, dataStr, centerPosPerTimeStep, excludeOutliers, planePos, u, v,...
  TF)
% required vectors
curPos = [];
cellIDs = [];
cPos = [];
nPos = [];
for j=1:dimData
  if cellDatas{j, 5} == curT
    if cellDatas{j, 7} ~= 0
      continue;
    end
    % this is a special case for the data set 130508 for which we
    % ignore the two cells that arise in the master cell file with
    % lineage ID 5
    if excludeOutliers == 1
      if strcmp( dataStr, '130508_raw' ) && cellDatas{j, 6} == 5
        continue;
      end
    end
    
    % get position of current cell
    p = [ cellDatas{j, 2} cellDatas{j, 3} cellDatas{j, 4} ];
    p = p - centerPosPerTimeStep(curT,:);
    p = applyTransformations( p, planePos, u, v, TF, dataStr, 1 );
    curPos = [ curPos ; p ];
    cellIDs = [ cellIDs ; cellDatas{j, 1} ];
  end
end
% go over the whole loop again to match successor cells based on
% lineage information
% the intention is to create two vectors cPos and nPos that share the same
% dimension in such a way that for each i a cell moved from cPos(i) to
% nPos(i)
nextPos = [];
for j=1:dimData
  % also check next time step
  if cellDatas{j, 5} == nextT
    if cellDatas{j, 7} ~= 0
      continue;
    end
    % this is a special case for the data set 130508 for which we
    % ignore the two cells that arise in the master cell file with
    % lineage ID 5
    if excludeOutliers == 1
      if strcmp( dataStr, '130508_raw' ) && cellDatas{j, 6} == 5
        continue;
      end
    end
    
    % get position of next cell
    p = [ cellDatas{ j, 2 } cellDatas{j, 3} cellDatas{j, 4} ];
    p = p - centerPosPerTimeStep(nextT,:);
    p = applyTransformations( p, planePos, u, v, TF, dataStr, 1 );
    
    % ID of current cell
    cellID = cellDatas{j, 1};
    % get precursor list
    [ precurIDList, numEntries ] = getPrecursorIDList( cellDatas{j, 8} );
    
    % if the cell is a daughter cell of a non-division cell
    % then the pos of the next time step is located at the same index
    % position as in the previous time step
    [row,col] = find( cellIDs == cellID );
    if size( row, 1 ) > 0
      cPos = [ cPos ; curPos(row, :) ];
      nPos = [ nPos ; p ];
      % else there is a division and we store both positions of
      % the two daughter cells to later generate a cell that is associated with
      % the cell at the previous time step
    else
      pre = numEntries;
      while pre > 0
        precurID = precurIDList( pre, 1 );
        [row,col] = find( cellIDs == precurID );
        if size( row, 1 ) > 0
          cPos = [ cPos ; curPos(row, :) ];
          nPos = [ nPos ; p ];
          break;
        else
          pre = pre - 1;
        end
      end
    end
  end
end