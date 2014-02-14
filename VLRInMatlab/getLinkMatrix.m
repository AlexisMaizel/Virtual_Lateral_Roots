function lMat = getLinkMatrix( p, q )
lMat = [ p(1)*q(1) p(1)*q(2) p(1)*q(3) ; p(2)*q(1) p(2)*q(2) p(2)*q(3) ; p(3)*q(1) p(3)*q(2) p(3)*q(3)];
