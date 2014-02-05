function lineageId = getLineageIdFromPos( x, y, z, t, data )

% get dimension aka number of lines
col = size(data{1});
numLines = col(1,1);

curPos = [ x y z ];

for j=1:numLines
  if data{5}(j) == t
    pos = [ data{2}(j) data{3}(j) data{4}(j) ];
    if pos == curPos
      lineageId = getLineageId( data{1}(j), data{7}(j) );
      break;
    end
  end
end
