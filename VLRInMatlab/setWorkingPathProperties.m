function setWorkingPathProperties()
% surpress warnings of using function names that already exist (TODO)
warning('off','all')
cd('C:\Jens\VLRRepository\VLRInMatlab')
addpath( genpath( 'geom3d/' ) );
addpath( genpath( 'export_fig/' ) );
addpath( genpath( '@tree/' ) );
addpath( genpath( 'BGLGraph/' ) );
addpath( genpath( 'ellipsoid_fit/' ) );
addpath( genpath( 'readTGMM_XMLoutput/' ) );
addpath( genpath( 'MinMaxFilterFolder/' ) );
addpath( genpath( 'MaximaMinima3D/' ) );
addpath( genpath( 'RidlerCalvardThresholding/' ) );
addpath( genpath( '3DSkeleton/' ) );
addpath( genpath( 'BuggenthinEtAl2013/' ) );