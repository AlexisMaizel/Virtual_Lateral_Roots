function texthandle = drawTriangleLabels( tri )
hold on
ic = incenter(tri);
numtri = size( tri,1 );
z = zeros( numtri, 1 );
trilabels = arrayfun(@(x) {sprintf('%d', x)}, (0:numtri-1)');
texthandle = text(ic(:,1), ic(:,2), z, trilabels, 'FontSize', 10, 'FontWeight', 'bold', ...
      'HorizontalAlignment', 'center', 'Color', 'blue');
hold off
end

