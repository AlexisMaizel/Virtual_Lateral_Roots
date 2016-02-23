function [ cellData, dimData ] = readRawData( dataStr )

storeOnlyLastPrecursorInfo = 0;

% reading raw data
path = strcat( '../FinalVLRForMatlab/', dataStr, '.csv' );
fileID = fopen( char(path) );
% format of data sets:
% ObjectID X Y Z Timepoint Radius Precursors Color Lineage TrackID TrackColor TrackGroup Layer DivisionType
formatSpec = '%d %f %f %f %d %d %q %q %d %q %q %d %d %q';
if strcmp( dataStr, '131203_raw' )
  formatSpec = '%d %f %f %f %d %d %q %q %d %d %d %q';
end
% read data and ignore the first four header lines
data = textscan( fileID, formatSpec, 'HeaderLines', 4, 'Delimiter', ';' );
fclose(fileID);

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
% z coord
ZCol = data{4};
% Time Step
TCol = data{5};
% string of precursors
PCol = data{7};
% Lineage Tree Id
LCol = data{9};
% temporal for data 131203
if strcmp( dataStr, '131203_raw' )
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
l = 1;
dimData = 0;
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
    dimData = dimData + deltaTS + 1;
    l = l+2;
  else
    dimData = dimData+1;
    l = l+1;
  end
end

cellData = cell( dimData, 8 );
cellFileMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'int32' );
divisionMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );

% data if daughter cells of a division
childrenData = [];

l = 1;
nl = 1;
% interpolate the missing positions in between and store the results in a
% tree structure as well as a cellData array
while (l < numLines+1)
  firstCellId = IdCol(l);
  secondCellId = IdCol(l+1);
  
  % check last precursor to detect daughter cells and the id of
  % their last parent
  % if the id is -1 than the cell starts at the root which we use as an
  % identification of non daughter cells
  % we set this value for all points interpolated in between
  lastPrecur = getLastPrecursorID( PCol(l) );
  
  % insert first line
  if storeOnlyLastPrecursorInfo == 1
    cellData( nl, : ) = {firstCellId XCol(l) YCol(l) ZCol(l) TCol(l) LCol(l) CFCol(l) lastPrecur};
  else
    cellData( nl, : ) = {firstCellId XCol(l) YCol(l) ZCol(l) TCol(l) LCol(l) CFCol(l) PCol(l)};
  end
  nl = nl+1;
  
  if lastPrecur ~= -1
    childrenData = [ childrenData ; double(lastPrecur) double(TCol(l)) XCol(l) YCol(l) ZCol(l) ];
  end
  
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
      if storeOnlyLastPrecursorInfo == 1
        cellData( nl, : ) = {firstCellId x y z t LCol(l) CFCol(l) -1};
      else
        cellData( nl, : ) = {firstCellId x y z t LCol(l) CFCol(l) PCol(l)};
      end
      nl = nl+1;
    end
    
    % insert last line
    if storeOnlyLastPrecursorInfo == 1
      cellData( nl, : ) = {firstCellId XCol(l+1) YCol(l+1) ZCol(l+1) TCol(l+1) LCol(l+1) CFCol(l+1) -1};
    else
      cellData( nl, : ) = {firstCellId XCol(l+1) YCol(l+1) ZCol(l+1) TCol(l+1) LCol(l+1) CFCol(l+1) PCol(l+1)};
    end
    nl = nl+1;
    % update cell file map
    cellFileMap( firstCellId ) = CFCol(l);
    
    % update division information
    divisionMap( firstCellId ) = [ XCol(l+1) YCol(l+1) ZCol(l+1) ];
    
    % increment loop index
    l = l+2;
    % else cell exists only for one time step
  else
    % update cell file map
    cellFileMap( firstCellId ) = CFCol(l);
    
    % update division information
    divisionMap( firstCellId ) = [ XCol(l) YCol(l) ZCol(l) ];
    
    l = l+1;
  end
end

% number of divisions
numDivisions = size( childrenData, 1 )/2.;

% store the object id, time step, the position, the division orientation and the
% division type for later processing
divisionData = zeros( numDivisions, 12 );

d = 1;
c = 1;
while d<2*numDivisions
  % object id of the division cell
  divisionData( c, 1 ) = childrenData( d, 1 );
  % time step of division cell
  divisionData( c, 2 ) = childrenData( d, 2 ) - 1;
  % position of division node
  divisionData( c, 3:5 ) = divisionMap(childrenData( d, 1 ));
  % orientation of division
  divisionData( c, 6:8 ) = [ childrenData( d, 3 ) childrenData( d, 4 ) childrenData( d, 5 ) ];
  divisionData( c, 9:11 ) = [ childrenData( d+1, 3 ) childrenData( d+1, 4 ) childrenData( d+1, 5 ) ];
  % TODO: division type
  %divisionData( c, 12 )
  d = d + 2;
  c = c + 1;
end
