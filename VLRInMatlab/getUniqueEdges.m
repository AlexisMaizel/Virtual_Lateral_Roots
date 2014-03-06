function uniqueEdges = getUniqueEdges( tri )

  dim = size( tri, 1 );
  edgeMatrix = [ dim*6:2 ];
  
  for i=1:dim
    count = 1;
    for j=1:3
      for k=j+1:4
        if tri( i, j ) < tri( i, k )
          edgeMatrix( (i-1)*6 + count, 1 ) = tri( i, j );
          edgeMatrix( (i-1)*6 + count, 2 ) = tri( i, k );
        else
          edgeMatrix( (i-1)*6 + count, 2 ) = tri( i, j );
          edgeMatrix( (i-1)*6 + count, 1 ) = tri( i, k );
        end
        count = count + 1;
      end
    end
  end
  % at last store only the unique edges and sort them
  uniqueEdges = unique( edgeMatrix, 'rows' );