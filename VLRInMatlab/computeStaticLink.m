function M = computeStaticLink( p1, nVec, centerData, tri, matPos, triangulationType )
% number of links of current cell
numLinks = size( nVec, 2 );

% initialize texture with zeros for each cell
M = zeros(3);

% loop over all linked neighbors
for n=1:numLinks
  % vertex ID
  verID = nVec(1,n);
  
  % get the corresponding neighbor position
  if triangulationType == 1
    p2 = [ tri.Points( verID, 1 ) tri.Points( verID, 2 ) tri.Points( verID, 3 ) ];
  else
    p2 = [ matPos( verID, 1 ) matPos( verID, 2 ) matPos( verID, 3 ) ];
  end
  
  p2 = p2 - centerData;
  
  % compute link matrix
  lMat = getLinkMatrix( p2-p1, p2-p1 );
  
  % update texture matrix
  M = M + lMat;
end

% after processing each neighbor, divide each entry by number
% of neighbors -> averaging
M = M./numLinks;
