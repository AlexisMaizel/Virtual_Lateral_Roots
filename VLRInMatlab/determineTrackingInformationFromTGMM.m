% execute compileMex.m before to generate mex files

TGMMPath = strcat( pwd, '/readTGMM_XMLoutput' );
addpath( TGMMPath );

% output format of values
format longG

% path to image output
imageDir = strcat( 'images/Segmentation/' );
mkdir( char(imageDir) );

resultPath = 'I:\SegmentationResults\TGMM\';
result = 'GMEMtracking3D_2016_2_18_10_1_56';
totalPath = strcat( resultPath, result, '\XML_finalResult_lht\GMEMfinalResult_frame');
startT = 1;
endT = 10;
radEllip = 10;
lineWidth = 1;
numColors = 20;

% svIdxCell:		cell array of length N, where N is the total number of objects tracked. svIdxCell{i} contains the indexes of the supervoxels belonging to the i-th object in trackingMatrix array. This index is necessary to rescue the segmentation from the .svb files output by TGMM software. The index starts in 0 following C convention. 
% trackingMatrix:		numericall array of size Nx10, where N is the number of points tracked bt TGMM over time. Each of the columns contains the following information:
% 
% 1.	Unique Id from the database to identify the point ( a large integer number)
% 2.	Cell type (represented by an integer). It is 0 if no cell type has been set for this object.
% 3.	x location of the nucleus centroid in world coordinates. Use the variable stackRes to convert from world coordinates to pixel unites.
% 4.	Same as 3 but for y location.
% 5.	Same as 3 but for z location.
% 6.	Estimated radius of the nucleus. It is 0 if this parameter was not estimated.
% 7.	Id of the cell in the previous time point. It is -1 if there is no linkage. Otherwise it has the unique id of the parent from column 1, so you can reconstruct the lineage.
% 8.	Time point of the nucleus.
% 9.	Confidence level in the tracking result. Value of 3 indicates high confidence that the object was correctly tracked. Value of 0 indicates low confidence.
% 10.	Skeleton id. All cells belonging to the same lineage have the same unique skeleton id.
[trackingMatrix, svIdxCell] = parseMixtureGaussiansXml2trackingMatrixCATMAIDformat( totalPath, startT, endT );

cmap = colorcube(numColors);

f = figure( 'Name', 'Segmentation' );
ELLIP = [];
ELLIPPATCH = [];
nucleiCounter = 1;

for t=startT:endT
  clf(f)
  xMinMax = [ -50 750 ];
  yMinMax = [ -450 -50 ];
  % axis([xmin xmax ymin ymax zmin zmax cmin cmax])
  axis( [ xMinMax(1) xMinMax(2) yMinMax(1) yMinMax(2) -10000 10000 0 1 ] );
  axis on
  daspect( [ 1 1 1 ] );
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  camproj( 'orthographic' );
%   hideHandle( ELLIP );
%   hideHandle( ELLIPPATCH );
  hold on
  for i=1:size(trackingMatrix,1)
    timeStep = trackingMatrix( i, 8 );
    if timeStep == t
      lineage = trackingMatrix( i, 10 );
      cen = trackingMatrix( i, 3:5 );
      color = cmap( mod( lineage, numColors )+1, : );
      
      [ ELLIP(nucleiCounter), ELLIPPATCH(nucleiCounter) ] =...
        drawEllipse3d( cen(1), -cen(2), cen(3), radEllip, radEllip, 0, 0 );
      set( ELLIP(nucleiCounter), 'color', color, 'LineWidth', lineWidth );
      set( ELLIPPATCH(nucleiCounter), 'FaceColor', color, 'FaceLighting', 'none' );
      nucleiCounter = nucleiCounter+1;
    end
  end
  
  tit = strcat( 'TimeStep', {' '}, num2str(t) );
  title( tit );
  
  % image output options
  if t < 10
    digit = strcat( 'TimeStep', '_00' );
  elseif t < 100
    digit = strcat( 'TimeStep', '_0' );
  else
    digit = strcat( 'TimeStep', '_' );
  end
  
  filePath = strcat( imageDir, digit, num2str(t), '.png' );
  export_fig( gcf, char(filePath), '-m2', '-png' );
end
