function [lineColorIndex, linePos, minMaxEigenValueIndex,...
  positiveEigenvalueVector, minMaxSemiAxisVector, centerEllipse,...
  timePositions, indexColorSet, contributions, magnitudes ]...
  = computeTimeEvolution( uniqueEdgesC, uniqueEdgesN, cellIdsC, cellIdsN,...
  numCellsN, triC, triN, matPosC, matPosN, cellPrecursorsN, triangulationType,...
  termTypeStr, dataStr, planePos, u, v, TF, deltaT, centerPosPerTimeStep,...
  curTC, curTN, renderMasterFile, cellFileMap )
% scaling of ellipses
scaling = 5;
% number of subdivisions for ellipsoids
nEllip = 10;

% first get mapping of vertex ids of delaunay triangulation to object
% ids of raw data set and store the results as a dimx2 matrix
% including the link information between two cells
objectLinksC = generateLinksBetweenObjectsId( uniqueEdgesC, cellIdsC );
objectLinksN = generateLinksBetweenObjectsId( uniqueEdgesN, cellIdsN );

% determine the number of conserved links for all cells
conservedLinks = intersect( objectLinksC, objectLinksN, 'rows' );

% disappeared links
disappearedLinks = setxor( objectLinksC, conservedLinks, 'rows' );

% appeared links
appearedLinks = setxor( objectLinksN, conservedLinks, 'rows' );

% iterating over the loop in order to determine the number of considered
% cells
numConsideredCells = 0;
for c=1:numCellsN
  % get objectId of current cell
  objectIdN = cellIdsN( c );
  % only render the master cell file if required
  if renderMasterFile == 1 &&...
      0 ~= cellFileMap( objectIdN )
    continue;
  end
  
  numConsideredCells = numConsideredCells + 1;
end

lineColorIndex = zeros( numConsideredCells, 9 );
linePos = zeros( numConsideredCells, 18 );
minMaxSemiAxisVector = zeros( numConsideredCells, 6 );
centerEllipse = zeros( numConsideredCells, 3 );
timePositions = zeros( numConsideredCells, 6 );
minMaxEigenValueIndex = zeros( numConsideredCells, 3 );
indexColorSet = zeros( numConsideredCells, 2 );
positiveEigenvalueVector = zeros( numConsideredCells, 3 );

% contributions for B and T term realted to their sum B+T
% first entry is |B|/(|B|+|T|) and second one is |T|/(|B|+|T|)
contributions = zeros( numConsideredCells, 2 );

% magnitudes of longest elongations of ellipses
% for geometrical and topological parts
% first entry is L_B/(L_B+L_T) and second one is L_T/(L_B+L_T)
magnitudes = zeros( numConsideredCells, 2 );

% consider each cell between two time steps
% and render the time evolution as an ellipsoid located at the
% position of the cell at time step t + deltaT
% Note that we consider the TIME STEP t + deltaT and look back at
% time step t which cell is related to the second time step
nc = 1;
for c=1:numCellsN
  % get objectId of current cell
  objectIdN = cellIdsN( c );
  % only render the master cell file if required
  if renderMasterFile == 1 &&...
      0 ~= cellFileMap( objectIdN )
    continue;
  end
  
  % get position of current cell
  if triangulationType == 1
    p1 = [ triN.Points( c, 1 ) triN.Points( c, 2 ) triN.Points( c, 3 ) ];
  else
    p1 = [ matPosN( c, 1 ) matPosN( c, 2 ) matPosN( c, 3 ) ];
  end
  
  p1 = p1 - centerPosPerTimeStep(curTN,:);
  
  % check if the current cell already existed in the last time
  % step; if not then back traverse its precursors until the
  % corresponding object id is found
  if isConserved( objectIdN, objectLinksC ) == 0
    objectIdC = getPrecursorID( objectIdN, cellPrecursorsN{ c }, cellIdsC );
  else
    objectIdC = objectIdN;
  end
  
  % get the position of the previous cell in the previous time step
  p2 = getCellPosition( objectIdC, triC, cellIdsC,...
    triangulationType, matPosC );
  
  p2 = p2 - centerPosPerTimeStep(curTC,:);
  
  % get all neighbors for the two found cells
  nVecC = getNeighborsOfObjectId( objectIdC, objectLinksC );
  nVecN = getNeighborsOfObjectId( objectIdN, objectLinksN );
  
  % vector of objectIds of cells which are conserved in the next time
  % step
  conservedLinksPerCell = getNeighborsOfObjectId( objectIdN, conservedLinks );
  % number of conserved links between two time steps
  numConservedLinksPerCell = size( conservedLinksPerCell, 2 );
  
  % get vector of objectsIds of cell which are added in the next time
  % step
  appearedLinksPerCell = getNeighborsOfObjectId( objectIdN, appearedLinks );
  
  % get vector of objectsIds of cell which disappeared in the next time
  % step
  disappearedLinksPerCell = getNeighborsOfObjectId( objectIdN, disappearedLinks );
  
  % number of links of current and next cell
  numTotalLinksPerCellC = size( nVecC, 2 );
  numTotalLinksPerCellN = size( nVecN, 2 );
  
  if numTotalLinksPerCellC == 0 && numTotalLinksPerCellN == 0
    continue;
  end
  
  % number of averaged links
  numAveragedLinksPerCell = ( numTotalLinksPerCellC + numTotalLinksPerCellN )/2.;
  
  % geometrical term
  if strcmp( termTypeStr, 'B' ) || strcmp( termTypeStr, 'All' )
    B = computeGeometricalLink( p1, p2, triC, triN, matPosC, matPosN,...
      cellIdsC, cellIdsN, deltaT, conservedLinksPerCell, numConservedLinksPerCell,...
      numAveragedLinksPerCell, triangulationType,...
      centerPosPerTimeStep, curTC, curTN );
  end
  
  % topological term
  if strcmp( termTypeStr, 'T' ) || strcmp( termTypeStr, 'All' )
    T = computeTopologicalLink( p1, p2, triC, triN, matPosC, matPosN,...
      cellIdsC, cellIdsN, objectIdC, objectIdN, objectLinksC, deltaT,...
      appearedLinksPerCell, disappearedLinksPerCell,...
      numAveragedLinksPerCell, triangulationType,...
      centerPosPerTimeStep, curTC, curTN );
  end
  
  % compute the final texture for the current ellipsoid
  if strcmp( termTypeStr, 'B' )
    M = B;
  elseif strcmp( termTypeStr, 'T' )
    M = T;
  else
    M = B + T;
  end
  
  % check if the matrix is zero then draw no ellipse
  if all( M == 0 )
    continue;
  end
  
  if strcmp( termTypeStr, 'All' )
    % compute the contributions of the single terms
    BContr = computeFrobeniusNorm(B);
    TContr = computeFrobeniusNorm(T);
    sumContr = BContr + TContr;%computeFrobeniusNorm(B+T);
    contributions( nc, 1 ) = BContr;%/sumContr;
    contributions( nc, 2 ) = TContr;%/sumContr;
    
    % compute the magnitudes of the single terms
    magnitudes( nc, 1 ) = determineMagnitude( B );
    magnitudes( nc, 2 ) = determineMagnitude( T );
    %lengthBT = lengthB + lengthT;%determineMagnitude( B+T )
  end
  
  % compute the eigenvectors and eigenvalues of matrix M
  % The columns of Q are the eigenvectors and the diagonal
  % elements of D are the eigenvalues
  [Q,D] = eig(M);
  
  % check if the eigenvalues are smaller than zero; if so, then do
  % not draw a line and consider the absolute value of it -> TODO
  positiveEigenvalue = [ 1 ; 1 ; 1 ];
  radii = diag(D);
  for e=1:3
    if radii(e) < 0.
      radii(e) = -radii(e);
      positiveEigenvalue( e, 1 ) = 0;
    end
  end
  
  % store the order of increasing eigen values
  [ ~, index ] = sort( radii );
  minMaxEigenValueIndex(nc, :) = [ index(1) index(2) index(3) ];
  
  positiveEigenvalueVector(nc, :) = [ positiveEigenvalue(1)...
    positiveEigenvalue(2) positiveEigenvalue(3) ];
  
  % line color by default black for each line
  lineColor = [ 0 0 0 0 0 0 0 0 0 ];
  % set line colors depending on the computed term and
  % eigenvalue
  if strcmp( termTypeStr, 'B' )
    for e=1:3
      % if the eigenvalue is positive then use a red color
      if positiveEigenvalue( e, 1 ) == 1
        lineColor( 1, 1 +3*e-3 ) = 1;
        lineColor( 1, 2 +3*e-3 ) = 0;
        lineColor( 1, 3 +3*e-3 ) = 0;
        % negative eigenvalue
      elseif positiveEigenvalue( e, 1 ) == 0
        lineColor( 1, 1 +3*e-3 ) = 0;
        lineColor( 1, 2 +3*e-3 ) = 0;
        lineColor( 1, 3 +3*e-3 ) = 1;
      end
    end
  end
  
  % store the coloring of lines depending on the sign of the
  % eigenvalues
  lineColorIndex(nc, :) = lineColor;
  
  % radii of the ellipsoid
  radii = sqrt( radii );
  
  % scaling such that even with small deformation changes
  % the ellipsoids can be identified
  radii = radii.*scaling;
  
  xEigVec = Q(:, 1);
  yEigVec = Q(:, 2);
  zEigVec = Q(:, 3);
  
  % draw the single cell as ellipsoid
  [ x, y, z ] = ellipsoid( 0, 0, 0, radii(1)/2., radii(2)/2., radii(3)/2., nEllip );
  %ellipPos = (p1+p2)/2.;
  ellipPos = p1;
  X = ellipPos(1) + x*xEigVec(1) + y*yEigVec(1) + z*zEigVec(1);
  Y = ellipPos(2) + x*xEigVec(2) + y*yEigVec(2) + z*zEigVec(2);
  Z = ellipPos(3) + x*xEigVec(3) + y*yEigVec(3) + z*zEigVec(3);
  
  semiLines = zeros( 1, 18 );
  % draw the three major axes in the origin which are then
  % rotated according to the eigenvectors
  for l=1:3
    if l == 1
      sX = [ -radii(1)/2., radii(1)/2. ];
      sY = [ 0, 0 ];
      sZ = [ 0, 0 ];
    elseif l == 2
      sX = [ 0, 0 ];
      sY = [ -radii(2)/2., radii(2)/2. ];
      sZ = [ 0, 0 ];
    else
      sX = [ 0, 0 ];
      sY = [ 0, 0 ];
      sZ = [ -radii(3)/2., radii(3)/2. ];
    end
    lineX = ellipPos(1) + sX*xEigVec(1) + sY*yEigVec(1) + sZ*zEigVec(1);
    lineY = ellipPos(2) + sX*xEigVec(2) + sY*yEigVec(2) + sZ*zEigVec(2);
    lineZ = ellipPos(3) + sX*xEigVec(3) + sY*yEigVec(3) + sZ*zEigVec(3);
    
    % TODO
    projLine1 = applyTransformations( [ lineX(1) lineY(1) lineZ(1) ], planePos, u, v, TF, dataStr, renderMasterFile, 1 );
    projLine2 = applyTransformations( [ lineX(2) lineY(2) lineZ(2) ], planePos, u, v, TF, dataStr, renderMasterFile, 1 );
    
    % and store the start/end points of the lines in linePos
    if l == 1
      semiLines(1, 1:6) = [ projLine1 projLine2 ];
    elseif l == 2
      semiLines(1, 7:12) = [ projLine1 projLine2 ];
    else
      semiLines(1, 13:18) = [ projLine1 projLine2 ];
    end
  end
  
  indexColorSet(nc, :) = determineColorAssignment(...
    [ xEigVec(1) xEigVec(2) xEigVec(3) ;...
    yEigVec(1) yEigVec(2) yEigVec(3) ;...
    zEigVec(1) zEigVec(2) zEigVec(3) ], positiveEigenvalue,...
    [ projLine1 ; projLine2 ] );
  
  % update semi axes in 3D
  linePos(nc, :) = semiLines;
  
  % project each vertex of the ellipsoid onto the plane
  dimP = size( X, 1 );
  for q=1:dimP
    for p=1:dimP
      curPos = applyTransformations( [ X(p,q) Y(p,q) Z(p,q) ], planePos, u, v, TF, dataStr, renderMasterFile, 1 );
      X(p,q) = curPos(1);
      Y(p,q) = curPos(2);
      Z(p,q) = curPos(3);
    end
  end
  
  ellipPos = applyTransformations( ellipPos, planePos, u, v, TF, dataStr, renderMasterFile, 1 );
  pStart = applyTransformations( p2, planePos, u, v, TF, dataStr, renderMasterFile, 1 );
  pEnd = applyTransformations( p1, planePos, u, v, TF, dataStr, renderMasterFile, 1 );
  timePositions(nc, :) = [ pStart pEnd ];
  centerEllipse(nc, :) = ellipPos;
  % the direction is now the normal of the x-y plane
  minMaxS = determineAxes( X, Y, Z, ellipPos, [ 0 0 1 ] );
  minMaxSemiAxisVector(nc, :) = minMaxS;
  nc = nc + 1;
end