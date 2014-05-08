function [ Min, Max ] = determineMinMax( array )

colSize = size( array, 2 );
if colSize == 1
Min = min(array);
Max = max(array);
elseif colSize == 2
  Min = 100000000;
  Max = -100000000;
  minM1 = min( array( :, 1 ) );
  minM2 = min( array( :, 2 ) );
  maxM1 = max( array( :, 1 ) );
  maxM2 = max( array( :, 2 ) );
  if minM1 < minM2
    minM = minM1;
  else
    minM = minM2;
  end
  if maxM1 > maxM2
    maxM = maxM1;
  else
    maxM = maxM2;
  end
  
  if minM <= Min
    Min = minM;
  end
  if maxM > Max
    Max = maxM;
  end
end

