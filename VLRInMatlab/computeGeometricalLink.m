function B = computeGeometricalLink( p1, p2, triC, triN, matPosC, matPosN,...
  cellIdsC, cellIdsN, deltaT, conservedLinksPerCell, numConservedLinksPerCell,...
  numAveragedLinksPerCell, triangulationType,...
  centerPosPerTimeStep, curTC, curTN )
B = zeros(3);
% loop over all linked and conserved neighbors
% determining B
for n=1:numConservedLinksPerCell
  % current object id of neighbor
  neighborId = conservedLinksPerCell( 1, n );
  
  % link at time step t
  l1 = getCellPosition( neighborId, triC, cellIdsC,...
    triangulationType, matPosC ) - p2;
  % link at time step t + deltaT
  l2 = getCellPosition( neighborId, triN, cellIdsN,...
    triangulationType, matPosN ) - p1;
  
  % translate points to center
  l1 = l1 - centerPosPerTimeStep(curTC,:);
  l2 = l2 - centerPosPerTimeStep(curTN,:);
  
  % compute average of the two links
  la = ( l1 + l2 )/2.;
  
  % compute difference between the two links
  ld = l2 - l1;
  
  % compute link matrix which is the matrix C in the paper (Appendix C1)
  C = getLinkMatrix( la, ld/deltaT );
  
  % sum of C and the transposed one
  B = B + C + C';
end

% after processing each neighbor, divide each entry by number
% of conserved neighbors -> averaging
if numConservedLinksPerCell > 0
  B = B./numConservedLinksPerCell;
end

% multiply the factor of N_c/N_tot (see paper in Appendix C1)
if numAveragedLinksPerCell > 0
  B = B.*( numConservedLinksPerCell / numAveragedLinksPerCell);
end