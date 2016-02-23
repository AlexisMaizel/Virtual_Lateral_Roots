% obj = readXMLmixtureGaussians(filenameXML);
% INPUT:
% 
% filenameXML:		string containing the path to the XML file generated by TGMM code. For exmaple, if the TGMM oputput is located at E:\TGMMruns\GMEMtracking3D_2014_6_3_15_39_21\XML_finalResult_lht and you want to read time point 100, then filenameXML = 'E:\TGMMruns\GMEMtracking3D_2014_6_3_15_39_21\XML_finalResult_lht\GMEMfinalResult_frame0100.xml'
% 
% OUTPUT:
% 
% obj:			Matlab structure of length N, where N is the number of objects detected in the time point specified by filenameXML. For each object, the structure obj(i) contains the following fields:
% 
% -id [integer]: 			unique id of the object in this particular time point. Following C convention, the index of the first element is 0.
% -lineage [integer]: 		unique id of the cell lineage the object belongs to.
% -parent [integer]: 		id of the linked object at the previous time point. Following the chain of �parent� objects reconstructs the track. A value of -1 indicates the birth of a track. Following C convention, the index of the first element is 0.
% -splitScore [float]: 		confidence level for the correct tracking of this particular object. A value of 0 indicates very low confidence and a value of 5 indicates very high confidence. Sorting elements by confidence level can guide the user in the data curation process and facilitate more effective editing of the TGMM results (see main text and Fig. 4).
% -scale [double[3]]: 		voxel scaling factors along the x-, y- and z-axis.
% -nu, beta, alpha [float]:	value of the hyper-parameters for the Bayesian GMM.
% -m [double[3]]: 		mean of the Gaussian Mixture (object centroid, in pixels).
% -W [double[3][3]]: 		precision matrix of the Gaussian Mixture (object shape).
% -*Prior: 			same as before, but for prior values obtained from the previous time point. These values are used during the inference procedure.
% -svIdx [integer[]]: 		list of indices of the super-voxels clustered by this Gaussian. Together with the �.svb� file, this information can be used to obtain precise segmentation regions for each object. Following C convention, the index of the first element is 0.
% -dims [integer]:		number of dimensions of the image in each time point. This value should be 2 or 3.

function readXMLmixtureGaussians()