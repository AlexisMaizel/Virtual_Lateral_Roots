function T = computeTopologicalLink( p1, p2, triC, triN, matPosC, matPosN,...
  cellIdsC, cellIdsN, objectIdC, objectIdN, objectLinksC, deltaT,...
  appearedLinksPerCell, disappearedLinksPerCell,...
  numAveragedLinksPerCell, triangulationType,...
  centerPosPerTimeStep, curTC, curTN )
  % topological term
  T1 = zeros(3);
  T2 = zeros(3);
  
  % number of added links between two time steps
  numAppearedLinksPerCell = size( appearedLinksPerCell, 2 );
  
% loop over all appeared linked neighbors determining the first part of T
for a=1:numAppearedLinksPerCell
  % current object id of neighbor
  neighborId = appearedLinksPerCell( 1, a );
  
  % link at time step t + deltaT
  pos1 = getCellPosition( neighborId, triN, cellIdsN,...
    triangulationType, matPosN );
  pos1 = pos1 - centerPosPerTimeStep(curTN,:);
  
  pos2 = getCellPosition( objectIdN, triN, cellIdsN,...
    triangulationType, matPosN );
  pos2 = pos2 - centerPosPerTimeStep(curTN,:);
  
  l = pos1 - pos2;
    
  % compute link matrix
  m = getLinkMatrix( l, l );
  
  % update matrix T
  T1 = T1 + m;
end

% NEW: also add the link between the cell in time step t + deltaT
% and the precursor which does not have the same id
if isConserved( objectIdN, objectLinksC ) == 0
  % link between position at time step t and t + deltaT
  l = p2 - p1;
  
  % compute link matrix
  m = getLinkMatrix( l, l );
  
  T1 = T1 + m;
  
  numAppearedLinksPerCell = numAppearedLinksPerCell + 1;
end

% after processing each neighbor, divide each entry by number
% of appeared neighbors -> averaging
if numAppearedLinksPerCell > 0
  T1 = T1./numAppearedLinksPerCell;
end

% multiply the factor of deltaN_a/N_tot (see paper in Appendix C1)
if numAppearedLinksPerCell > 0
  T1 = T1.*( numAppearedLinksPerCell / numAveragedLinksPerCell);
end

% divide by deltaT
T1 = T1./deltaT;

% number of added links between two time steps
numDisappearedLinksPerCell = size( disappearedLinksPerCell, 2 );

% loop over all disappeared linked neighbors determining the second part of T
for d=1:numDisappearedLinksPerCell
  % current object id of neighbor
  neighborId = disappearedLinksPerCell( 1, d );
  
  % link at time step t
  pos1 = getCellPosition( neighborId, triC, cellIdsC,...
    triangulationType, matPosC );
  pos1 = pos1 - centerPosPerTimeStep(curTC,:);
  
  pos2 = getCellPosition( objectIdC, triC, cellIdsC,...
    triangulationType, matPosC );
  pos2 = pos2 - centerPosPerTimeStep(curTC,:);
  
  l = pos1 - pos2;
  
  % compute link matrix
  m = getLinkMatrix( l, l );
  
  % update matrix T
  T2 = T2 + m;
end

% after processing each neighbor, divide each entry by number
% of disappeared neighbors -> averaging
if numDisappearedLinksPerCell > 0
  T2 = T2./numDisappearedLinksPerCell;
end

% multiply the factor of deltaN_a/N_tot (see paper in Appendix C1)
if numDisappearedLinksPerCell > 0
  T2 = T2.*( numDisappearedLinksPerCell / numAveragedLinksPerCell);
end

% divide by deltaT
T2 = T2./deltaT;

T = T1 - T2;