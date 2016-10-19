Virtual_Lateral_Roots Repository
================================
Date of writing: 09/26/2016
Author: Jens Fangerau

This is a repository storing the VLR datasets and all corresponding work for
modelling (in VV/C++), visual analysis (in scifer/C++) segmenting (in Matlab/Fiji)
and tracking (in Matlab) virtual lateral root datasets.

The content of each subfolder is explained below:


-------------------------------
FijiScripts
-------------------------------

This folder contains some scripts written in Jython that are used in Fiji to apply some
image processing algos to VLR datasets. The name of the files and therefore methods are
self-explanatory and paths of the datasets are used with respect to working on the
Wolverine workstation. The overall idea was to apply a preprocessing step for further
segmenting and tracking cells over time. Preprocessing refers to steps like automatic or
manual thresholding, creating of binary masks, etc. . But also MIP generation of single or
time series data was implemented. More details are given in the corresponding files.


-------------------------------
FinalVLRForMatlab
-------------------------------

This folder contains all five VLR datasets manually segmented/tracked by Daniel von Wangenheim
given by their time of aquisition. In contrast to the folder LRPAnalysisResults which contains
the original data results from Daniel, here, the data was processed in scifer to automatically
classify division types as well as cell files and these results are also stored in the datasets
sheets. The difference to the datasets given in the folder FinalVLRForScifer is that here a
delimiter is given by a ";" which is a " " in the datasets in folder FinalVLRForScifer.


-------------------------------
FinalVLRForScifer
-------------------------------

See description of FinalVLRForMatlab.


-------------------------------
LRPAnalysisResults
-------------------------------

See description of FinalVLRForMatlab. Additionally, the data is already interpolated between
division events as well as synchronized based on number of cells, for example.


-------------------------------
Modelling
-------------------------------

Contains all code of the simulation/modelling of a VLR in 2D realized in C++/VV. VV (vertex-vertex)
is a programming language provided by Richard Smith to generate a growing model of any 2D
structure given an initial and a last surface structure. The surface/tissue is realized using a 2D
Bezier surface. Together with the initial and last surface type plus the definition when and how
a cell should divide, a complete 2D VLR model can be established and run over time. The
installation of VV (Vlab on Linux and maybe LStudio on Windows) is not that trivial and requires
adittionally libs like Qt, libGLEW etc. with a minimum version number. For an updated installation
giude I refer to Richard Smith who provided us all necessary files and software.


-------------------------------
VLRInMatlab
-------------------------------

This folder contains all analysis and visualization of the VLR realized in Matlab. The subfolders
are external libraries or methods used to simplify own analysis methods. *.m files that do NOT start
with "main" are helper methods that are used within the main Matlab files to simplfy code reading.
All files are more or less documented and work by others are cited within the file. All methods are either
applied to Daniel's datasets (120830, 121204, 121211, 130508, 130607) or to the "new" datasets
generated recently (20160427, 20160428, 20160426). Below, a description of all *.m files
starting with "main" are given:


mainAverageDeformation
**********************

This method is inspired by the idea of Graner et al (2008): Discrete rearranging disordered
patterns, part I: Robust statistical tools in two or three dimenions to compute the deformation
of cells over time based on their underlying (delaunay) triangulation depicting the cell neighborhood.
This method was used to create the visual average deformation results in Wangenheim et al (2016):
Rules and Self-Organizing Properties of Postembryonic Plant Organ Cell Division Pattern.

mainExportTriangulation
**********************

Old method to manually generate contour points of the bezier surfaces that are used as surfaces
for the 2D modelling. The initial idea was to approximate the contour of the five VLR datasets
but we ended up using only an averaged and "optimal" 2D surface structure for the model for
comparison with the real data.

mainGraphNetworkTopology
**********************

This method was implemented to analysis and visualize graph properties of underlying triangulations
of the five VLR datasets. Four graphs properties were computed: Degree, Closeness, Betweenness and
PageRank (Random Walk). The contribution of graph property values are analyzed for each tile of a
grid that was overlayed on the viewing type of the data (side, top or radial). This is similar to
measurement of the contributions done in mainAverageDeformation.

mainMembraneSegmentationMultiProcessor
**********************

This algo applies a membrane segmentation inspired by the idea of Kawase et al (2015):
A direction-selective local-thresholding method, DSLT, in combination with a dye-based method
for automated three-dimensional segmentation of cells and airspaces in developing leaves.
Additionally, a cell shape extractor is included that enables the detection and generation of
cell shapes based on the segmented membrane information.

mainNormalizedStatic
**********************

This method realizes the analysis of a STATIC deformation for a single time step described in
Graner et al (2008): Discrete rearranging disordered patterns, part I: Robust statistical tools
in two or three dimenions. Each time step is normalized over all datasets. See also
mainStaticTexture.

mainNormalizedTime
**********************

This method realizes the analysis of a time-based deformation for several time steps described in
Graner et al (2008): Discrete rearranging disordered patterns, part I: Robust statistical tools
in two or three dimenions. It was mainly used as a playground to try different visual analysis
methods and which is also used as a basis for the method in mainAverageDeformation. Each time
step is normalized over all datasets. See also mainTimeEvolution.

mainNucleiSegmentation
**********************

Complete nuclei segmentation step from, thresholding, different methods for identifying cell objects
as for example using the membrane channel, nearest neighbor search or clustering of regional maxima
and saving the result as TIFF.

mainRegisteredBezierAndGrowthDisplacementModel
**********************

Method that computes the growth displacements (change of cell positions), register them among
all five datasets based one the number of cells and create a bezier surface for each time
step that can be used for the 2D modelling.

mainRegisteredBezierAndGrowthTensorModel
**********************

Same as mainRegisteredBezierAndGrowthDisplacementModel but instead of incorporating the growth
displacement information, the tensor values of the B (geometrical) and T (topological) terms
are included (from mainAverageDeformation).

mainRegisteredGridApproachModelAveraging
**********************

Using the grid/contribution approach of mainAverageDeformation to quantify the deformation of cells
within tiles in a grid among all five datasets.

mainRegisteredKNNApproachModelAveraging
**********************

Same as mainRegisteredGridApproachModelAveraging but instead of the grid approach, a KNN search
is applied.

mainSegmentationParameterOptimization
**********************

Initial nuclei segmentation step with numerical comparison between manual and automatic result. The
automatic approach was mainly generated with the Fiji scripts and then postprocessed in Matlab to
identify connected components and to apply the comparison.

mainSimulationDataGraphTopology
**********************

Same as mainGraphNetworkTopology but applied to the simulation/modelling data for which 100 runs
were performed.

mainStaticTexture
**********************

The main idea is the same as in mainNormalizedStatic BUT there is one main difference:
mainStaticTexture -> for each dataset, traverse over all time steps and compute the static deformation
mainNormalizedStatic -> for each time step, traverse over all datasets and compute the static deformation

mainTimeEvolution
**********************

The main idea is the same as in mainNormalizedTime BUT there is one main difference:
mainStaticTexture -> for each dataset, traverse over all time steps and compute the temporal deformation
mainNormalizedTime -> for each time step, traverse over all datasets and compute the temporal deformation

mainTrackingInformationFromTGMM
**********************

This method allows a batch parameter setting using the TGMM software to segment nuclei. It also
features a visual comparison to the real datasets as well as a simple lineage visualization.


-------------------------------
Experience_with_Segmentation_Tools.txt
-------------------------------

Text file that gives a list of properties, PROs, CONs about other segmentation and tracking tools based on my
experience and applied to our datasets.

-------------------------------
SciferCode
-------------------------------

C++ code fragments taken from scifer that were used for the analysis and visualization of the VLR datasets. An
additional README file within the folder gives more details about the corresponding source and header files.

-------------------------------
NewDatasetsProperties.pdf
-------------------------------

A small sheet giving some properties about the recently acquired datasets: 20160426, 20160427, 20160428.
