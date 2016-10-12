Fragments of Scifer code (software that is used in my former working group of visualization)
related to the visualization of the Virtual Lateral Root
========================================================
Date of writing: 10/12/2016
Author: Jens Fangerau

Here you can find some information about the content of the scifer code written in C++ and for
what purpose it was written for. Note that this code is not compilable, that dependencies could be missing and it is just copied from the main software scifer. It should just help in understanding how I realized different approaches for the analysis of the VLR. The concept of scifer is that it has a GUI and that different tasks are given as algorithms. All my work with respect to the analysis and visualization of the VLR is done in an algorithm called VLRBrowser. Within scifer, if this algo is called, a menu is shown with entires, buttons, combo boxes, etc. in which the user can choose several parameters and therefore decide the kind of visualization/analysis.

-------------------------------
LoadArabidopsisLineages
-------------------------------

This is a code fragment which was used to load/import the cell and lineage information of Daniel's data stored in the folder "FinalVLRForScifer". The das is then converted into the binary tree structure used in sicfer for faster pastprocessing and traversal of the lineage information.

-------------------------------
VLRBrowser
-------------------------------

This is the main algorithm where all analysis and visualization of the VLR starts. All the user-chosen parameters are set in the class method "SVLRBrowser::initProfile()" that might help to get an idea what the purpose of a parameter is. If a setting is chosen the corresponding code is executed that may refer to the execution of other classes given in "AnalysisClasses".

-------------------------------
AnalysisClasses
-------------------------------

Here are several classes that are used within the algo VLRBrowser (if the user has chosen them):

SCorrelationAnalysis
*******************************
All kind of information was extracted from the data and a visual correlation analysis was done using a scatterplot matrix to find correlations between several features.

SCellLayers
*******************************
This class manages the generation of cell layers that were used to automatically classify the division types (anticlinal, periclinal, radial) during the development of the VLR.

SLayerInformationIO
*******************************
All kind of import and export methods to integrate the information of layers and division schemes into the data files (also used for editing etc.) stored in the folder "FinalVLRForScifer".

SSimilarityMeasureGraphics
*******************************
Graphics class that features a set of several graphical parts (realized in OpenSceneGraph/OpenGL) that are used in all kind of visualizations realted to the VLR.

SSimilarityMeasureHeader
*******************************
A header storing several typedefs, etc. used all over the code.

SVLRDataAnalysisIO
*******************************
Import and export of intermediate data that was generated withtin the algorithm.

SCellTreeVisualizer
*******************************
2D color-coded Lineage tree visualizer of all loaded lineage data.

SCellLayerValidation
*******************************
Provides several tools to visually validate the generated cell layers like clipping, MIPs, bounding boxes.

SSimilarityMeasureUtil
*******************************
Namespace having lots of methods to compute and quantify several aspects of moving cells like cell cycle lengths, division direction, spherical coordinates, division angles.

SConvexHull
*******************************
Class for generating the convex hull of data points realized with CGAL.

SAlphaShape
*******************************
Class for generating the alpha shape of data points realized with CGAL.

SDivisionAnalysis
*******************************
Class that focuses on different approaches for analyzing and visualizing division properties.

SCurvatureAnalysis
*******************************
Class for analyzing and visualizing the curvature information of VLR data given as the convex hull or alpha shape of the data points.

SCellAnalysis
*******************************
Class for analyzing cell growth and movement directions over time.

SVoronoiDiagram
SVoronoiCell
SRenderVoronoiCells
*******************************
Class for generation, handling, and visualizing a voronoi diagram of cell data.
