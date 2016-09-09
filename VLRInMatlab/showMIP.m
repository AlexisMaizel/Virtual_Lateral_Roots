function h = showMIP( imageStack )
mipStack = max( imageStack, [], 3 );
restoredefaultpath
h = imshow(mipStack, []);
setWorkingPathProperties()