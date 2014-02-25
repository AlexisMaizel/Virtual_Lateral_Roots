function planeCoord = projectOnPlane( x, r, u, v )
  % x is the coord in 3D space which should be projected
  % r is the support vector of the plane
  % u is the first spanning vector
  % v is the seccond spanning vector
  pos = ( dot( x, u )/dot( u, u ) ) * u + ( dot( x, v )/dot( v, v ) ) * v;
  
  % after finding the point on the plane tranlate it by the magnitude of r
  T = [ 1 0 0 r(1) ; 0 1 0 r(2) ; 0 0 1 r(3) ; 0 0 0 1 ];
  planeCoord = T * [ pos ; 1 ];
  planeCoord = [ planeCoord(1) planeCoord(2) planeCoord(3) ];
end