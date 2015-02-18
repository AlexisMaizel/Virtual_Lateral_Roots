function [min, max] = updateMinMax( pos, min, max )

dim = size( pos, 2 );
for i=1:dim
  if pos(1, i) < min(i)
    min(i) = pos(1, i);
  end
  
  if pos(1, i) > max(i)
    max(i) = pos(1, i);
  end
end