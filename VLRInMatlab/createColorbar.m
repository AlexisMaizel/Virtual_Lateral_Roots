function COLORBAR = createColorbar( colorType )
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
else
  colors = [ 222./255., 235./255., 247./255. ;
  49./255., 130./255., 189./255. ];
% colors = [ 255./255., 255./255., 255./255. ;
%   0./255., 0./255., 255./255. ];
end
n = 256; % size of new color map
m = size(colors,1);
t0 = linspace(0,1,m)';
t = linspace(0,1,n)';
r = interp1(t0,colors(:,1),t);
g = interp1(t0,colors(:,2),t);
b = interp1(t0,colors(:,3),t);
cmap = [r,g,b];
cm = colormap( cmap );
%   numTicks = 9;
%   steps = (maxContr-minContr)/numTicks;
%   for l=1:numTicks+1
%     labels{l} = [ num2str(minContr+(l-1)*steps) ];
%   end
COLORBAR = colorbar( 'location', 'eastoutside' );
set( COLORBAR, 'YTick', 0:1:2,  'YTickLabel', [ 5 10 ] );