function determineMSERRegions()
fileName = 'I:\SegmentationResults\Preprocessing\20160426\MIPraw_T3.jpg';
I = imread( char(fileName) );
regions = detectMSERFeatures(I);
figure;
restoredefaultpath
imshow(I);
setWorkingPathProperties()
hold on;
plot(regions); % by default, plot displays ellipses and centroids