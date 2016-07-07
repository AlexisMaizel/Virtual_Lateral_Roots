function ccLocalMaximaMap = determineNumberOfLocalMaximaOfCC( CC, S, minVoxelCount, connectivity )
ccLocalMaximaMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );

% init adjacency matrix for each cc with zeros
for c=1:CC.NumObjects
  areaVal = S(c, :).Area;
  % ignore cc with pixeal area smaller than some threshold
  if areaVal < minVoxelCount
    continue;
  end
  
  %diary('matlabOutput.txt')
  %[center, radii, evecs, pars ] = ellipsoid_fit( S(c,:).PixelList );
  %radii( 3, 1 ) = radii( 3, 1 )*anisotropyZ;
  %r = int16( mean( radii( 1:2, 1 ) ) )
  %r = int16(max( radii ))
  r = 5;
  
  list = S(c,:).PixelList;
  numPixels = size( list, 1 );
  bb = S(c,:).BoundingBox;
  
  xmin = min( list(:, 1) );
  ymin = min( list(:, 2) );
  zmin = min( list(:, 3) );
  
  % for each pixel object set an 0 else background is 1
  mat = ones( bb(1, 4), bb(1, 5), bb(1, 6) );
  for p=1:numPixels
    pos = list( p, : );
    pos = pos - [ xmin-1 ymin-1 zmin-1 ];
    mat( pos(1,1), pos(1,2), pos(1,3) ) = 0;
  end
  
  % compute the Euclidean distance between each zero and its nearest
  % non-zero element
  D = bwdist( mat );
%   Dmax = minmaxfilt( D, r, 'max', 'same' );
%   ip = find( Dmax == D );
%   [ x, y, z ] = ind2sub( size(D), ip );
%   idxkeep = find( x>1 & x<size(D,1) & y>1 & y<size(D,2) & z>1 & z<size(D,3));
%   x=x(idxkeep);
%   y=y(idxkeep);
%   z=z(idxkeep);
%   ip=ip(idxkeep);
  
  % find to local maxima in the distance matrix
  msk = true(r, r, r);
  mid = int32((r+1)/2);
  msk(mid, mid, mid) = false;
  
%   DMIP = zeros( size(D,1), size(D,2) );
%   for d1=1:size(D,1)
%     for d2=1:size(D,2)
%       for d3=1:size(D,3)
%         if D( d1, d2, d3 ) > DMIP( d1, d2 )
%           DMIP( d1, d2 ) = D( d1, d2, d3 );
%         end
%       end
%     end
%   end

  %diary('matlabOutput.txt')
  areaVal = S(c, :).Area;
  centroid = S(c, :).Centroid;
  
  %area = S(c, :).Area
  D_dil = imdilate(D, msk);
  DCC = zeros( size(D,1), size(D,2), size(D,3) );
  for d1=1:size(D,1)
    for d2=1:size(D,2)
      for d3=1:size(D,3)
        if D( d1, d2, d3 ) > 2 && D( d1, d2, d3 ) >= D_dil( d1, d2, d3 )
          DCC( d1, d2, d3 ) = 1;
        end
      end
    end
  end
  
%   DCC = zeros( size(D,1), size(D,2), size(D,3) );
%   for l=1:size(x,1)
%     DCC( x(l,1), y(l,1), z(l,1) ) = 1;
%   end

  NCC = bwconncomp( DCC, connectivity );
  NS = regionprops( NCC, 'Centroid', 'Area' );
  numCCs = size(NS, 1)
  
  % old version
%   locMax = [];
%   % compute distances between each centroid of cc
%   distances = zeros( size(NS, 1), size(NS, 1) );
%   mergedCC = cell(1, size(NS, 1));
%   for i=1:size(NS, 1)
%     p1 = NS(i,:).Centroid + [ xmin-1 ymin-1 zmin-1 ];
%     mergedCC{1, i} = [ i ];
%     for j=1:size(NS, 1)
%       if j > i
%         p2 = NS(j,:).Centroid + [ xmin-1 ymin-1 zmin-1 ];
%         distances( i, j ) = sqrt( (p1(1,1)-p2(1,1))*(p1(1,1)-p2(1,1)) + (p1(1,2)-p2(1,2))*(p1(1,2)-p2(1,2)) + (p1(1,3)-p2(1,3))*(p1(1,3)-p2(1,3)) );
%       end
%     end
%   end
  
%   
%   for i=1:size(NS, 1)
%     for j=1:size(NS, 1)
%       if j>i && distances( i, j ) < r
%         mergedCC{1, i} = [ mergedCC{1, i} mergedCC{1, j} ];
%         mergedCC{1, j} = [];
%       end
%     end
%   end
  
  X = [];
  for i=1:size(NS, 1)
    X = [ X ; NS(i,:).Centroid + [ xmin-1 ymin-1 zmin-1 ] ];
  end
  %[ idx, distances ] = rangesearch( X, X, r )
  Y = pdist(X);
  if size(NS, 1) > 1
    Z = linkage( Y, 'centroid' );
    T = cluster( Z, 'cutoff', r, 'criterion', 'distance' );
    k = max(T);
  elseif size(NS, 1) == 1
    k = 1;
  else
    k = 0;
  end
  
  %numCCs = size(NS, 1)
%   for i=1:length( mergedCC )
%     %mergedCC{ 1, i }
%     centroid = [ 0 0 0 ];
%     si = length( mergedCC{ 1, i } );
%     for j=1:si
%       ind = mergedCC{ 1, i }(j);
%       centroid = centroid + NS(ind, :).Centroid;
%     end
%     if si > 0
%       centroid = centroid/double(si);
%       centroid = centroid + [ xmin-1 ymin-1 zmin-1 ];
%       locMax = [ locMax ; centroid ];
%     end
%   end
  %cells = size( locMax, 1 )
  %diary off
  
  ccLocalMaximaMap(c) = k;
  %break;
end
