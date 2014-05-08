function length = determineMagnitude( M )
% check if the matrix is zero then return a length of zero
if all( M == 0 )
  length = 0.;
  return;
end

% compute the eigenvectors and eigenvalues of matrix M
% The columns of Q are the eigenvectors and the diagonal
% elements of D are the eigenvalues
[~,D] = eig(M);

% get eigenvalues
radii = diag(D);

% store the order of increasing eigen values
[ ~, index ] = sort( radii );

if -radii(index(1)) > radii(index(3))
  length = -radii(index(1));
  length = sqrt( length );
  length = -length;
else
  length = radii(index(3));
  length = sqrt( length );
end

