setWorkingPathProperties()

chosenData = 1;
dataStr = { 'BD' 'Random' };

simulationSteps = 1;
startT = 500;
endT = 500;

radEllip = 8;
lineWidth = 1;

% drawing parameters
drawGraph = 0;
drawDistributions = 1;
drawVoronoi = 0;

% graph coloring: 1 -> cell layering, 2 -> graph properties
graphColoring = 2;

numBins = 20;

% number of generated graph properties each displayed in a different
% subplot
numGraphProperties = 4;
graphPropString = { 'Degree' 'Closeness' 'Betweenness' 'PageRank' };

% output format of values
format longG

% figure properties
f = figure( 'Name', 'Graph Topology Analysis', 'Position', [10 10 1400 800] );
% activate orbit rotation by default
cameratoolbar( 'SetMode', 'orbit' );
% activate none coord system by default for not resetting the camera up
% vector when rotating
cameratoolbar( 'SetCoordSys', 'none' );
% show camera toolbar by default
cameratoolbar( 'Show' );

% gca is the current axes handle
set( gca,'nextplot','replacechildren' );
% gcf is the current figure handle
%lighting phong
set( gcf, 'Renderer', 'zbuffer' );
lighting gouraud
set( gcf, 'Renderer', 'OpenGL' );
set( gcf,'nextplot','replacechildren' );
set( gcf, 'color', [ 1 1 1 ] );

% surface instance
ELLIP = [];
ELLIPPATCH = [];
LINES = [];

imageOutputPath = strcat( 'I:\GraphTopologyAnalysis\', dataStr( 1, chosenData ), '\' );
mkdir( char(imageOutputPath) );

tic
% read simulation data
[ cellData, connMat ] = readSimulationData( dataStr( 1, chosenData ),...
  simulationSteps, endT );
toc

curS = 1;
nucleiCounter = 1;
tic
% loop over all time steps
for curT=startT:endT
  if curT < 10
    digit = '00';
  elseif curT < 100
    digit = '0';
  else
    digit = '';
  end
  
  hideHandle( ELLIP );
  hideHandle( ELLIPPATCH );
  hideHandle( LINES );
  
  for pl=1:numGraphProperties
    subplot( numGraphProperties/2, numGraphProperties/2, pl, 'replace' )
    if drawGraph == 1
      axis( [ -300 300 -50 150 ] );
    end
    axis on
    grid off;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title( strcat( dataStr( 1, chosenData ), ' Time Step ', num2str(curT), '\_', char( graphPropString( 1, pl ) ) ) );
    if drawGraph == 1
      daspect( [ 1 1 1 ] );
    end
  end
  
  % adjacency matrix of current time step
  adjaMat = connMat{ curS, 1 }{ curT, 1 };
  % and cell data information
  data = cellData{ curS, 1 }{ curT, 1 };
  
  % number of cells in the current time step
  numCells = size( adjaMat, 1 );
  
  % compute graph properties
  G = graph( adjaMat );
  
  % max number of cell layers
  numCellLayers = max( data( :, 7 )+1 );
  
  for pl=1:numGraphProperties
    subplot( numGraphProperties/2, numGraphProperties/2, pl )
    
    if pl == 1
      gp = centrality( G, 'degree' );
    elseif pl == 2
      gp = centrality( G, 'closeness' );
    elseif pl == 3
      gp = centrality( G, 'betweenness' );
    elseif pl == 4
      gp = centrality( G, 'pagerank' );
    end
    
    if graphColoring == 1
      cm = generateLayerColorMap();
    else
      cm = jet( numCells );
    end
    colormap( cm );
    minNC = min(gp);
    maxNC = max(gp);
    if minNC ~= maxNC
      cellIndex = floor( ((gp-minNC)/range( gp )) * (numCells-1) ) + 1;
    else
      cellIndex = gp;
    end
    
    cellLayerToGP = cell( numCellLayers, 1 );
    
    % draw an ellipsoid for each cell
    hold on
    for c=1:numCells
      % get position of current cell
      p = [ data( c, 4 ) data( c, 5 ) 10 ];
      layer = data( c, 7 )+1;
      
      % separate the gp values into each cell layer
      cellLayerToGP{ layer, 1 } = [ cellLayerToGP{ layer, 1 } gp( c, 1 ) ];
      
      if drawGraph == 1
        if graphColoring == 1
          color = cm( layer, : );
        else
          color = cm( cellIndex(c, 1), : );
        end
        [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
          drawEllipse3d( p(1), p(2), p(3), radEllip, radEllip, 0, 0 );
        set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
        set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
        nucleiCounter = nucleiCounter+1;
      end
    end
    
    % draw lines
    if drawGraph == 1
      gplot2( adjaMat, data( :, 4:5 ), '-k', 'LineWidth', lineWidth );
      colorbar
      caxis( [minNC, maxNC] )
    end
    
    if drawDistributions == 1
      drawStackedBar( minNC, maxNC, numCellLayers, numBins, cellLayerToGP,...
        pl, 1 );
    end
  end
  
  if drawGraph == 1
    suffixStr = '_Data';
  else
    suffixStr = '_Distribution';
  end
  
  filePath = strcat( imageOutputPath, dataStr( 1, chosenData ), '_T', digit, num2str(curT), suffixStr, '.png' );
  export_fig( gcf, char(filePath), '-m2', '-png' );
end
toc
