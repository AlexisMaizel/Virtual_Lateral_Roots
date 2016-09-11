function [] = drawStackedBar( minNC, maxNC, numStacks, numBins, data,...
  pl, type )
hold off
% generate histogram for each attribute (cell file or cell layer)
hi = cell( numStacks, 1 );
if type == 2
  for cf=1:numStacks
    hi{ cf, 1 } = hist( data{ cf, 1 }, numBins );
  end
else
  for cf=1:numStacks
    % get complete data from map container
    if isKey( data, cf ) == 1
      vals = cell2mat( data( cf ) );
      hi{ cf, 1 } = hist( vals, numBins );
    else
      hi{ cf, 1 } = hist( 0, numBins );
    end
  end
end
hold on
% 
x = linspace( minNC, maxNC, numBins );
y = zeros( numBins, numStacks );
for cf=1:numStacks
  for b=1:numBins
    y( b, cf ) = y( b, cf ) + hi{ cf, 1 }( 1, b );
  end
end
xplot = 1:numel(x);
ba = bar( xplot, y, 'stacked' );
NumTicks = 3;
L = get(gca,'XLim');
%xpos = linspace( L(1), L(2), NumTicks );
xpos = linspace( L(1), numBins, NumTicks );
xIndex = linspace( 1, numBins, NumTicks );
labels = cell( 1, NumTicks );

if pl == 2 || pl == 4
  prec = '% .4f';
else
  prec = '% .2f';
end

for t=1:NumTicks
  labels{ 1, t } = num2str(x( round(xIndex(1, t) ) ), prec );
end
set( gca, 'XTick', xpos, 'XTickLabel', cellstr(labels) )
% cell files
if type == 2
  ba(1).FaceColor = 'black'; % -4
  ba(2).FaceColor = 'black'; % -3
  ba(3).FaceColor = 'red'; % -2
  ba(4).FaceColor = 'green'; % -1
  ba(5).FaceColor = 'blue'; % 0
  ba(6).FaceColor = 'yellow'; % 1
  ba(7).FaceColor = 'cyan'; % 2
  ba(8).FaceColor = 'magenta'; % 3
  T = [ 0.5, 0, 0.5 ; 0, 0, 0 ; 1, 0, 0 ; 0, 1, 0 ; 0, 0, 1 ; 1, 1, 0 ; 0, 1, 1 ; 1, 0, 1 ];
  xsp = [ 0 ; 1/7 ; 2/7 ; 3/7 ; 4/7 ; 5/7 ; 6/7 ; 7/7 ];
  custMap = interp1( xsp, T, linspace(0,1,8));
  colormap( custMap )
else
  % cell layers
  cm = generateLayerColorMap();
  colormap( cm );
  for c=1:numStacks
     ba(c).FaceColor = cm( mod( c-1, size(cm, 1) )+1, : );
  end
end
colorbar
%caxis( [minCF, maxCF] )