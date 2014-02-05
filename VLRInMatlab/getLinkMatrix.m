function lMat = getLinkMatrix( p1, p2 )
l = p2 - p1;

lMat = [ l(1)*l(1) l(1)*l(2) l(1)*l(3) ; l(2)*l(1) l(2)*l(2) l(2)*l(3) ; l(3)*l(1) l(3)*l(2) l(3)*l(3)];
