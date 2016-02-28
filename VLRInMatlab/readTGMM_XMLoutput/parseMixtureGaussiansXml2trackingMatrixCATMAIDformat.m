% INPUT:
% 
% basename:		string containing the path to the XML files generated by TGMM code. For exmaple, if they are located at E:\TGMMruns\GMEMtracking3D_2014_6_3_15_39_21\XML_finalResult_lht, then basename = 'E:\TGMMruns\GMEMtracking3D_2014_6_3_15_39_21\XML_finalResult_lht\GMEMfinalResult_frame'
% frameIni:		integer with the initial time point to import (typically it is 0)
% frameEnd:		integer with the final time point to import
% 
% OUTPUT:
% 
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

function [trackingMatrix, svIdxCell, error]= parseMixtureGaussiansXml2trackingMatrixCATMAIDformat(basename,frameIni,frameEnd)

numTM = frameEnd - frameIni + 1;
trackingMatrix = zeros( 15000 * numTM, 10 );%preallocate memory
svIdxCell = cell( size(trackingMatrix,1), 1);%stores supervoxel idx for each element in trackingMatrix
trackingMatrixN = 0;
skeletonN = 0;
error = 0;
for frame=frameIni:frameEnd
    
   try
    obj=readXMLmixtureGaussians([basename num2str(frame,'%.4d') '.xml']);
   catch
     disp( 'An error occurred while reading the xml file' );
     error = 1;
     break;
   end
    
    %process each object
    mapId = zeros(length(obj),1);
    for ii=1:length(obj)
        blob=obj(ii);
        
        if( blob.m(1) < -1e31 )%dead cell
            continue;
        end
        
        if(isempty(blob.splitScore)) 
            blob.splitScore=-1e32;
        end;
        
        trackingMatrixN = trackingMatrixN + 1;
        mapId(ii) = trackingMatrixN;
        
        parId = blob.parent + 1;%matlab indexing        
        if( parId <= 0 )
            skeletonN = skeletonN + 1;
            skeletonId = skeletonN;           
            parNodeId = -1;
        else
            parId = mapIdOld( parId );
            skeletonId = trackingMatrix(parId,10);
            parNodeId = trackingMatrix(parId,1);
        end
        
        trackingMatrix(trackingMatrixN,:) = [trackingMatrixN -1 blob.m blob.betaPrior  parNodeId frame  blob.splitScore skeletonId];%used as a wildcard for different measuemrents: betaPrior has probBAckground
        svIdxCell{ trackingMatrixN } = blob.svIdx;                
    end              
    mapIdOld = mapId;
end

if error == 0
  trackingMatrix = trackingMatrix( 1:trackingMatrixN, : );
  svIdxCell = svIdxCell(1:trackingMatrixN);
end
