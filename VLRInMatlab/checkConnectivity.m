function [ connList ] = checkConnectivity( p, mat, max, connectivity )

connList = [];
if size(p, 2) == 3
  if connectivity == 6
    if p(1)+1 <= max(1)
      connList = [ connList mat( p(1)+1, p(2), p(3) ) ];
    end
    if p(1)-1 > 0
      connList = [ connList mat( p(1)-1, p(2), p(3) ) ];
    end
    if p(2)+1 <= max(2)
      connList = [ connList mat( p(1), p(2)+1, p(3) ) ];
    end
    if p(2)-1 > 0
      connList = [ connList mat( p(1), p(2)-1, p(3) ) ];
    end
    if p(3)+1 <= max(3)
      connList = [ connList mat( p(1), p(2), p(3)+1 ) ];
    end
    if p(3)-1 > 0
      connList = [ connList mat( p(1), p(2), p(3)-1 ) ];
    end
  elseif connectivity == 18
    if p(1)+1 <= max(1)
      connList = [ connList mat( p(1)+1, p(2), p(3) ) ];
    end
    if p(1)-1 > 0
      connList = [ connList mat( p(1)-1, p(2), p(3) ) ];
    end
    if p(2)+1 <= max(2)
      connList = [ connList mat( p(1), p(2)+1, p(3) ) ];
    end
    if p(2)-1 > 0
      connList = [ connList mat( p(1), p(2)-1, p(3) ) ];
    end
    if p(3)+1 <= max(3)
      connList = [ connList mat( p(1), p(2), p(3)+1 ) ];
    end
    if p(3)-1 > 0
      connList = [ connList mat( p(1), p(2), p(3)-1 ) ];
    end
    
    if p(1)+1 <= max(1) && p(2)+1 <= max(2)
      connList = [ connList mat( p(1)+1, p(2)+1, p(3) ) ];
    end
    if p(1)-1 > 0 && p(2)-1 > 0
      connList = [ connList mat( p(1)-1, p(2)-1, p(3) ) ];
    end
    if p(1)+1 <= max(1) && p(2)-1 > 0
      connList = [ connList mat( p(1)+1, p(2)-1, p(3) ) ];
    end
    if p(1)-1 > 0 && p(2)+1 <= max(2)
      connList = [ connList mat( p(1)-1, p(2)+1, p(3) ) ];
    end
    
    if p(1)+1 <= max(1) && p(3)+1 <= max(3)
      connList = [ connList mat( p(1)+1, p(2), p(3)+1 ) ];
    end
    if p(1)-1 > 0 && p(3)-1 > 0
      connList = [ connList mat( p(1)-1, p(2), p(3)-1 ) ];
    end
    if p(1)+1 <= max(1) && p(3)-1 > 0
      connList = [ connList mat( p(1)+1, p(2), p(3)-1 ) ];
    end
    if p(1)-1 > 0 && p(3)+1 <= max(3)
      connList = [ connList mat( p(1)-1, p(2), p(3)+1 ) ];
    end
    
    if p(2)+1 <= max(2) && p(3)+1 <= max(3)
      connList = [ connList mat( p(1), p(2)+1, p(3)+1 ) ];
    end
    if p(2)-1 > 0 && p(3)-1 > 0
      connList = [ connList mat( p(1), p(2)-1, p(3)-1 ) ];
    end
    if p(2)+1 <= max(2) && p(3)-1 > 0
      connList = [ connList mat( p(1), p(2)+1, p(3)-1 ) ];
    end
    if p(2)-1 > 0 && p(3)+1 <= max(3)
      connList = [ connList mat( p(1), p(2)-1, p(3)+1 ) ];
    end
  elseif connectivity == 26
    % TODO
  else
    disp( 'Connectivity is not supported!' )
  end
else
  disp( 'Not supported dimension of pos!' )
end
  