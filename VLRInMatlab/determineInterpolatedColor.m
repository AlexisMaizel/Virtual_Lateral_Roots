function color = determineInterpolatedColor( value, min, max )
% colors for contribution of B and T terms
colors = [ 222./255., 235./255., 247./255. ;
  %158./255., 202./255., 225./255. ;
  49./255., 130./255., 189./255. ];

t = (max-value)/(max-min);
color = [ t * colors(1, 1) + (1-t) * colors(2, 1)
          t * colors(1, 2) + (1-t) * colors(2, 2)
          t * colors(1, 3) + (1-t) * colors(2, 3) ];