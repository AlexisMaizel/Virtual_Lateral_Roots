function [ cellDatas, dimData, maxT, numCellsPerTimeStep, centerPosPerTimeStep,...
  totalMinAxes, totalMaxAxes, cellFileMap ] =...
  prepareData( dataStr, startData, endData, numData, visualizationType, renderSingleCellFile, cView )
cpuT = cputime;
% cellData is the main array with all relevant information
% for the further analysis:
% ObjectID X Y Z Timepoint LineageID TrackGroup
cellDatas = cell( numData );
cellFileMap = cell( numData );
maxT = zeros( numData, 1 );
minAxes = zeros( numData, 3 );
maxAxes = zeros( numData, 3 );
totalMinAxes = [ 5000 5000 5000 ];
totalMaxAxes = [ -5000 -5000 -5000 ];
numCellsPerTimeStep = cell( numData, 1 );
dimData = zeros( numData, 1 );
centerPosPerTimeStep = cell( numData, 1 );
% preprocessing of data sets to store and determine different properties
for dataIndex=startData:endData
  % reading raw data
  path = strcat( '../FinalVLRForMatlab/', dataStr( 1, dataIndex ), '.csv' );
  fileID = fopen( char(path) );
  % format of data sets:
  % ObjectID X Y Z Timepoint Radius Precursors Color Lineage TrackID TrackColor TrackGroup Layer DivisionType
  formatSpec = '%d %f %f %f %d %d %q %q %d %q %q %d %d %q';
  if strcmp( dataStr( 1, dataIndex ), '131203_raw' )
    formatSpec = '%d %f %f %f %d %d %q %q %d %d %d %q';
  end
  % read data and ignore the first four header lines
  data = textscan( fileID, formatSpec, 'HeaderLines', 4, 'Delimiter', ';' );
  fclose(fileID);
  
  % set PC depending on the viewing direction
  coeff = getPrincipalComponents( dataStr( 1, dataIndex ), renderSingleCellFile );
  if cView == 1
    dir = coeff(:,2);
    u = coeff(:,1);
    v = coeff(:,3);
  elseif cView == 2
    dir = coeff(:,3);
    u = coeff(:,1);
    v = coeff(:,2);
    if strcmp( dataStr( 1, dataIndex ), '121211_raw' )
      v = -v;
    end
    if strcmp( dataStr( 1, dataIndex ), '120830_raw' ) &&...
        renderSingleCellFile == 0
      v = -v;
    end
  elseif cView == 3
    dir = -coeff(:,1);
    u = -coeff(:,3);
    v = coeff(:,2);
    if strcmp( dataStr( 1, dataIndex ), '121211_raw' )
      v = -v;
    end
  end
  
  % set plane position
  planePos = dir * 1;
  plane = [ planePos(1) planePos(2) planePos(3)...
    u(1) u(2) u(3)...
    v(1) v(2) v(3) ];
  TF = createBasisTransform3d( 'g', plane );
  
  % get dimension aka number of lines
  col = size(data{1});
  numLines = col(1,1);
  % store relevant columns in corresponding data structures
  % and ignore the others
  % Object Id
  IdCol = data{1};
  % X coord
  XCol = data{2};
  % Y coord
  YCol = data{3};
  % z coord; Note that due to resampling, the z value is multiplied with 2
  if strcmp( dataStr( 1, dataIndex ), '121204_raw_2014' )
    ZCol = -2 * data{4};
  else
    ZCol = 2 * data{4};
  end
  % Time Step
  TCol = data{5};
  % string of precursors
  PCol = data{7};
  % Lineage Tree Id
  LCol = data{9};
  % temporal for data 131203
  if strcmp( dataStr( 1, dataIndex ), '131203_raw' )
    % Cell File Id
    CFCol = data{10};
    % Cell Layer Id
    %CLCol = data{11};
    % Divison Type
    %DCol = data{12};
  else
    % Cell File Id
    CFCol = data{12};
    % Cell Layer Id
    %CLCol = data{13};
    % Divison Type
    %DCol = data{14};
  end
  
  % store min and max of cell files
  %minCF = min( CFCol );
  %maxCF = max( CFCol );
  
  l = 1;
  % first determine dimension of cellData
  % in order to initialize the size of cellData
  % -> huge performance speed up
  while (l < numLines+1)
    firstCellId = IdCol(l);
    secondCellId = IdCol(l+1);
    if firstCellId == secondCellId
      firstTS = TCol(l);
      secondTS = TCol(l+1);
      deltaTS = secondTS - firstTS;
      dimData( dataIndex ) = dimData( dataIndex ) + deltaTS + 1;
      l = l+2;
    else
      dimData( dataIndex ) = dimData( dataIndex )+1;
      l = l+1;
    end
  end
  
  cellData = cell( dimData( dataIndex ), 8 );
  cellFileMap{ dataIndex } = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
  
  l = 1;
  nl = 1;
  % interpolate the missing positions in between
  while (l < numLines+1)
    firstCellId = IdCol(l);
    secondCellId = IdCol(l+1);
    % insert first line
    cellData( nl, : ) = {firstCellId XCol(l) YCol(l) ZCol(l) TCol(l) LCol(l) CFCol(l) PCol(l)};
    nl = nl+1;
    
    % interpolate between cell positions
    if firstCellId == secondCellId
      firstTS = TCol(l);
      secondTS = TCol(l+1);
      deltaTS = secondTS - firstTS;
      steps = double( 1./double(deltaTS) );
      
      for s=1:deltaTS-1
        t = firstTS + s;
        k = double( double(s)*steps );
        x = double(1-k) * XCol(l) + double(k) * XCol(l+1);
        y = double(1-k) * YCol(l) + double(k) * YCol(l+1);
        z = double(1-k) * ZCol(l) + double(k) * ZCol(l+1);
        % insert all relevant data into main data structure
        cellData( nl, : ) = {firstCellId x y z t LCol(l) CFCol(l) PCol(l)};
        nl = nl+1;
      end
      
      % insert last line
      cellData( nl, : ) = {firstCellId XCol(l+1) YCol(l+1) ZCol(l+1) TCol(l+1) LCol(l+1) CFCol(l+1) PCol(l+1)};
      nl = nl+1;
      % update cell file map
      cellFileMap{ dataIndex }( firstCellId ) = CFCol(l);
      
      % increment loop index
      l = l+2;
      % else cell exists only for one time step
    else
      % update cell file map
      cellFileMap{ dataIndex }( firstCellId ) = CFCol(l);
      
      l = l+1;
    end
  end
  
  cellDatas{ dataIndex } = cellData;
  
  % initialize number of rows depending on the number of time steps of
  % current data set
  maxT( dataIndex ) = max( TCol );
  numCellsPerTimeStep{dataIndex} = zeros(maxT(dataIndex), 1);
  centerPosPerTimeStep{dataIndex} = zeros(maxT(dataIndex), 3);
  MIN = [ 5000 5000 5000 ];
  MAX = [ -5000 -5000 -5000 ];

  % get maximum number of cells for each data set and time step
  % as well as determine center of each position set at a time step
  for j=1:dimData( dataIndex )
    timeStep = cellDatas{dataIndex}{j, 5};
    numCellsPerTimeStep{dataIndex}(timeStep, 1) = numCellsPerTimeStep{dataIndex}(timeStep, 1) + 1;
    pos = [ cellDatas{dataIndex}{j, 2} cellDatas{dataIndex}{j, 3} cellDatas{dataIndex}{j, 4} ];
    centerPosPerTimeStep{dataIndex}(timeStep, :) = centerPosPerTimeStep{dataIndex}(timeStep, :) + pos;
  end
  
  for t=1:maxT( dataIndex )
    centerPosPerTimeStep{dataIndex}(t, :) = centerPosPerTimeStep{dataIndex}(t, :)./numCellsPerTimeStep{dataIndex}(t, 1);
  end
  
  for j=1:dimData( dataIndex )
    timeStep = cellDatas{dataIndex}{j, 5};
    pos = [ cellDatas{dataIndex}{j, 2} cellDatas{dataIndex}{j, 3} cellDatas{dataIndex}{j, 4} ];
    pos = pos - centerPosPerTimeStep{dataIndex}(timeStep, :);
    pos = applyTransformations( pos, planePos, u, v, TF, dataStr( 1, dataIndex ) );
    
    % determine min and max values for each time step
    for m=1:3
      if pos(m) < MIN(m)
        MIN(m) = pos(m);
      end
      if pos(m) >= MAX(m)
        MAX(m) = pos(m);
      end
    end
  end

  if strcmp( visualizationType, 'Ellipses' ) || strcmp( visualizationType, 'Contour' )
    % set axes properties
    minAxes(dataIndex, :) = MIN;
    maxAxes(dataIndex, :) = MAX;
  elseif strcmp( visualizationType, 'Ellipsoids' )
    % get bounding box by determining
    % min and max of x/y/z values
    minX = min( XCol );
    minY = min( YCol );
    minZ = min( ZCol );
    maxX = max( XCol );
    maxY = max( YCol );
    maxZ = max( ZCol );
    minAxes(dataIndex, :) = [ minX minY minZ ];
    maxAxes(dataIndex, :) = [ maxX maxY maxZ ];
  end
  
  for mm=1:3
    if minAxes(dataIndex, mm) < totalMinAxes(mm)
      totalMinAxes(mm) = minAxes(dataIndex, mm);
    end
    if maxAxes(dataIndex, mm) >= totalMaxAxes(mm)
      totalMaxAxes(mm) = maxAxes(dataIndex, mm);
    end
  end
end
ElapsedTimeDataPreparation = cputime - cpuT