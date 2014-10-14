% output format of values
format longG
% read distance matrix given as vector
path = strcat( '/tmp/vec.txt' );
fileID = fopen( char(path) );
dist = dlmread( path );
%dim = size( dist, 1 )
%entries = 502*501/2
fclose(fileID);
% read data
path = strcat( '/tmp/data.txt' );
fileID = fopen( char(path) );
data = dlmread( path, ' ' );
fclose(fileID);

linkageStr = { 'average' 'centroid' 'complete' 'median' 'single' 'ward', 'weighted' };

lnTypeS = 1;
lnTypeE = 7;

v = zeros(7, 1);

for l=lnTypeS:lnTypeE
  %X = pdist( data );
  if size( dist, 1 ) == 1
    X = dist;
  else
    X = dist';
  end
  
  %Y = randn( 550, 1 );
  %X = pdist( Y );
  Z = linkage( X, linkageStr( 1, l ) );
  %dendrogram(Z);
  % cophenetic correlation coefficient
  v(l,1) = cophenet( Z, X );
end

v