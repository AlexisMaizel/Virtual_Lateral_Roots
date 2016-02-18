% execute compileMex.m before to generate mex files

TGMMPath = strcat( pwd, '/readTGMM_XMLoutput' );
addpath( TGMMPath );

resultPath = 'I:\SegmentationResults\TGMM\';
result = 'GMEMtracking3D_2016_2_17_15_3_21';
totalPath = strcat( resultPath, result, '\XML_finalResult_lht\GMEMfinalResult_frame');
startT = 1;
endT = 10;
consideredT = 1;

[trackingMatrix, svIdxCell] = parseMixtureGaussiansXml2trackingMatrixCATMAIDformat( totalPath, startT, endT );

figure
axis equal
hold on
for i=1:size(trackingMatrix,1)
  timeStep = trackingMatrix( i, 8 );
  if timeStep == consideredT
    cen = trackingMatrix( i, 3:5 );
    [x,y,z] = ellipsoid( cen(1), -cen(2), cen(3), 10, 10, 10 );
    surf(x,y,z)
    shading interp
  end
end
colorbar