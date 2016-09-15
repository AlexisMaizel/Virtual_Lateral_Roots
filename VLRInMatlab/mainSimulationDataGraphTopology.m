setWorkingPathProperties()

chosenData = 2;
dataStr = { 'BD' 'Random' };

startT = 1;
stepT = 20;
endT = 500;
% startS should always start with 1!!!!
startS = 1;
endS = 100;

radEllip = 8;
lineWidth = 1;

% resolution of grid
resGrid = 20;

% drawing parameters
drawAverageGraph = 1;
drawAverageDistribution = 0;
draw2DGraph = 0;
drawGrid = 0;
drawVoronoi = 0;

storeImages = 1;

% determine global min and max values
determineMinMaxValues = 0;

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
GRIDTILES = [];

imageOutputPath = strcat( 'I:\GraphTopologyAnalysis\', dataStr( 1, chosenData ), '\' );
mkdir( char(imageOutputPath) );

disp( 'Reading simlulation data' )
tic
maxS = endS;
maxT = 500;
% read simulation data
[ cellData, connMat ] = readSimulationData( dataStr( 1, chosenData ),...
  maxS, maxT );
toc

totalMinAxes = [ -300 -50 -10 ];
totalMaxAxes = [ 300 150 10 ];
nucleiCounter = 1;
                       
if determineMinMaxValues == 1
  globalMinNC = ones( numGraphProperties, 1 ).*realmax;
  globalMaxNC = zeros( numGraphProperties, 1 );
else
  % if the graph is drawn then the min and max values are different because
  % of computing the mean value which is not the case for the distibution
  % plot
  % MEAN
%   if drawAverageGraph == 1
%     % averaged min and max values for each graph property
%     globalMinNC = [ 1 ; 0.00299954965288543 ; 0.00531914893617021 ; 0.00977880469131011 ];
%     globalMaxNC = [ 6.52 ; 1 ; 186.936937678646 ; 0.5 ];
%   else
%     % averaged min and max values for each graph property
%     globalMinNC = [ 1 ; 0.0022271714922049 ; 0.333333333333333 ; 0.00767674555703311 ];
%     globalMaxNC = [ 10 ; 1 ; 511.429512096124 ; 0.5 ];
%   end
  % VARIANCE
  if drawAverageGraph == 1
    % averaged min and max values for each graph property
    globalMinNC = [ 0.0114926987560844 ; 7.52316384526264e-37 ; 5.91645678915759e-31 ; 7.70371977754894e-34 ];
    globalMaxNC = [ 2.88888888888889 ; 0.00114237035833957 ; 16929.3290544878 ; 0.000972376745430563 ];
  else
    % averaged min and max values for each graph property
    globalMinNC = [ 1 ; 0.0022271714922049 ; 0.333333333333333 ; 0.00767674555703311 ];
    globalMaxNC = [ 10 ; 1 ; 511.429512096124 ; 0.5 ];
  end
  % BD
%   globalMinNC =
% 
%         0.0114926987560844
%       7.52316384526264e-37
%       5.91645678915759e-31
%       7.70371977754894e-34
% 
% 
% globalMaxNC =
% 
%           1.51246537396122
%       2.12039263552961e-05
%           7140.95067435156
%       0.000256099284644677
% Random
% globalMinNC =
% 
%         0.0131555555555556
%       2.30377511592669e-10
%        0.00263128112267995
%        6.1862096868431e-08
% 
% 
% globalMaxNC =
% 
%           2.88888888888889
%        0.00114237035833957
%           16929.3290544878
%       0.000972376745430563
end

disp( 'Traversing time steps' )
tic
% loop over all time steps
for curT=startT:stepT:endT
  tileGrid = cell( numGraphProperties, 1 );
  cellLayers = cell( numGraphProperties, 1 );
  hideHandle( ELLIP );
  hideHandle( ELLIPPATCH );
  hideHandle( LINES );
  hideHandle( GRIDTILES );
  
  for pl=1:numGraphProperties
    subplot( numGraphProperties/2, numGraphProperties/2, pl, 'replace' )
    if drawGrid == 1
      % initialize 2D grid
      [ rows, columns ] =...
        generate2DGrid(...
        [totalMinAxes(1) totalMinAxes(2)],...
        [totalMaxAxes(1) totalMaxAxes(2)], resGrid );
    end
    
    tileGrid{ pl, 1 } = cell( rows*columns, 1 );
    cellLayers{ pl, 1 } = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
    % remove all entries at the beginning of each time step
    %remove( cellLayers{ pl, 1 }, keys( cellLayers{ pl, 1 } ) );
    if drawAverageGraph == 1 || draw2DGraph == 1
      axis( [ -300 300 -50 150 ] );
    end
    axis on
    grid off;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title( strcat( dataStr( 1, chosenData ), ' Time Step ', num2str(curT), '\_', char( graphPropString( 1, pl ) ) ) );
    if drawAverageGraph == 1 || draw2DGraph == 1
      daspect( [ 1 1 1 ] );
    end
  end
  
  % loop over all simulation steps
  for curS=startS:endS    
    % adjacency matrix of current time step
    adjaMat = connMat{ curS, 1 }{ curT, 1 };
    % and cell data information
    data = cellData{ curS, 1 }{ curT, 1 };
    
    % number of cells in the current time step
    numCells = size( adjaMat, 1 );
    
    % compute graph properties
    G = graph( adjaMat );
    
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
      
      if draw2DGraph == 1
        if graphColoring == 1
          cm = generateLayerColorMap();
        else
          numColors = 20;
          cm = jet( numColors );
        end
        colormap( cm );
        minNC = min(gp);
        maxNC = max(gp);
        if minNC ~= maxNC
          colorIndex = round( ((gp-minNC)/range( gp )) * (numColors-1) ) + 1;
        else
          colorIndex = ones( size( gp ) );
        end
      end
      
      % draw an ellipsoid for each cell
      hold on
      for c=1:numCells
        % get position of current cell
        p = [ data( c, 4 ) data( c, 5 ) 10 ];
        layer = data( c, 7 )+1;
        % separate the gp values into each cell layer
        if isKey( cellLayers{ pl, 1 }, layer ) == 1
          tempCell = cellLayers{ pl, 1 }( layer );
          dim = size( tempCell, 1 );
          tempCell{ dim+1, 1 } = gp( c, 1 );
          cellLayers{ pl, 1 }( layer ) = tempCell;
        else
          tempCell = cell( 1, 1 );
          tempCell{ 1, 1 } = gp( c, 1 );
          cellLayers{ pl, 1 }( layer ) = tempCell;
        end
        
        % locate the cell within the grid and assign its gp value to the
        % grid in order to average over the data in each tile
        tileIndex = getTileIndex( p, [totalMinAxes(1) totalMinAxes(2)],...
          [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns );
        tileGrid{ pl, 1 }{ tileIndex } = [ tileGrid{ pl, 1 }{ tileIndex }; gp( c, 1 ) ];
        
        if draw2DGraph == 1
          if graphColoring == 1
            color = cm( mod( layer-1, size(cm, 1) )+1, : );
          else
            color = cm( colorIndex( c, 1 ), : );
          end
          [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
            drawEllipse3d( p(1), p(2), p(3), radEllip, radEllip, 0, 0 );
          set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
          set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
          nucleiCounter = nucleiCounter+1;
        end
      end
      
      % draw lines
      if draw2DGraph == 1
        gplot2( adjaMat, data( :, 4:5 ), '-k', 'LineWidth', lineWidth );
        lineCounter = lineCounter + 1;
        colorbar
        if minNC ~= maxNC
          caxis( [minNC, maxNC] )
        end
      end
    end
  end
  
  if determineMinMaxValues == 1
    for pl=1:numGraphProperties
      for gt=1:rows*columns
        numValues = size( tileGrid{ pl, 1 }{ gt }, 1 );
        if numValues ~= 0
          %averageVal = mean( tileGrid{ pl, 1 }{ gt } );
          averageVal = var( tileGrid{ pl, 1 }{ gt }, 1 );
          
          if averageVal > 0
            if averageVal < globalMinNC( pl, 1 )
              globalMinNC( pl, 1 ) = averageVal;
            end
          end
          
          if averageVal >= globalMaxNC( pl, 1 )
            globalMaxNC( pl, 1 ) = averageVal;
          end
        end
      end
    end
  end
  
  if drawAverageGraph == 1
    % combine all simulation results for the current time step
    % get overall min and max values of the graph properties
    numColors = 20;
    cm = jet( numColors );
    colormap( cm );
    hold on
    tileCounter = 1;
    for pl=1:numGraphProperties
      subplot( numGraphProperties/2, numGraphProperties/2, pl )
      %tg = cell2mat( tileGrid{ pl, 1 } );
      %totalMinNC = min( tg, [], 1 );
      %totalMaxNC = max( tg, [], 1 );
      totalMinNC = log( globalMinNC( pl, 1 ) );
      totalMaxNC = log( globalMaxNC( pl, 1 ) );
      colorbar
      if totalMinNC ~= totalMaxNC
        caxis( [totalMinNC, totalMaxNC] )
      end
      for gt=1:rows*columns
        numValues = size( tileGrid{ pl, 1 }{ gt }, 1 );
        if numValues ~= 0
          %averageVal = mean( tileGrid{ pl, 1 }{ gt } );
          averageVal = var( tileGrid{ pl, 1 }{ gt }, 1 );
           if averageVal == 0
             averageVal = globalMinNC( pl, 1 );
           end
           averageVal = log( averageVal );
          if totalMinNC ~= totalMaxNC
            colorIndex = round( ((averageVal-totalMinNC)/(totalMaxNC-totalMinNC)) * (numColors-1) ) + 1;
          else
            colorIndex = 1;
          end
          color = cm( colorIndex, : );
          GRIDTILES( tileCounter ) = drawColoredTile( gt,...
            [totalMinAxes(1) totalMinAxes(2)],...
            [totalMaxAxes(1) totalMaxAxes(2)], resGrid, rows, columns, color );
          tileCounter = tileCounter + 1;
        end
      end
    end
  end
  
  if drawAverageDistribution == 1
    for pl=1:numGraphProperties
      subplot( numGraphProperties/2, numGraphProperties/2, pl )
      %tg = cell2mat( tileGrid{ pl, 1 } );
      %totalMinNC = min( tg, [], 1 );
      %totalMaxNC = max( tg, [], 1 );
      if pl ~= 1
        totalMinNC = log( globalMinNC( pl, 1 ) );
        totalMaxNC = log( globalMaxNC( pl, 1 ) );
      else
        totalMinNC = globalMinNC( pl, 1 );
        totalMaxNC = globalMaxNC( pl, 1 );
      end
      drawStackedBar( totalMinNC, totalMaxNC, size( cellLayers{ pl, 1 }, 1 ),...
        numBins, cellLayers{ pl, 1 }, pl, 1 );
    end
  end
  if storeImages == 1
    if drawAverageGraph == 1
      suffixStr = '_Data';
    else
      suffixStr = '_Distribution';
    end
    filePath = strcat( imageOutputPath, dataStr( 1, chosenData ), num2str(curT),...
      '_T', num2str(startT), '-', num2str(endT),...
      '_S', num2str(startS), '-', num2str(endS), suffixStr, '.png' );
    export_fig( gcf, char(filePath), '-m2', '-png' );
  end
end

toc

if determineMinMaxValues == 1
  globalMinNC
  globalMaxNC
end
