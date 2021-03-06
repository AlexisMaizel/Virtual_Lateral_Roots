function color = determineInterpolatedColor( value, min, max, colorType )
if colorType == 0
% colors for contribution of B and T terms -> green gradient
colors = [ 229./255., 245./255., 224./255. ;
  49./255., 163./255., 84./255. ];
% colors = [ 255./255., 255./255., 255./255. ;
%   0./255., 255./255., 0./255. ];
% positive magnitude -> red gradient
elseif colorType == 1
  colors = [ 254./255., 224./255., 210./255. ;
  222./255., 45./255., 38./255. ];
%   colors = [ 255./255., 255./255., 255./255. ;
%   255./255., 0./255., 0./255. ];
% negative magnitude -> blue gradient
elseif colorType == 2
  colors = [ 222./255., 235./255., 247./255. ;
  49./255., 130./255., 189./255. ];
% colors = [ 255./255., 255./255., 255./255. ;
%   0./255., 0./255., 255./255. ];
else
  color = [ 0 0 0 ];
  return;
end

if colorType ~= 2
  t = (max-value)/(max-min);
  color = [ t * colors(1, 1) + (1-t) * colors(2, 1)
    t * colors(1, 2) + (1-t) * colors(2, 2)
    t * colors(1, 3) + (1-t) * colors(2, 3) ];
else
  t = (max-value)/(max-min);
  color = [ (1-t) * colors(1, 1) + t * colors(2, 1)
    (1-t) * colors(1, 2) + t * colors(2, 2)
    (1-t) * colors(1, 3) + t * colors(2, 3) ];
end