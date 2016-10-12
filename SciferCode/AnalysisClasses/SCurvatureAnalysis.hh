#ifndef SCURVATUREANALYSIS_HH
#define SCURVATUREANALYSIS_HH

/**
  @file   SCurvatureAnalysis.hh
  @brief  Contains class for analyzing the curvature of VLR data
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SSimilarityMeasureGraphics.hh"

class SCurvatureAnalysis
{
public:
  SCurvatureAnalysis( osg::ref_ptr<SSliderCallback> sliderCallback,
                      boost::shared_ptr<cellTimeVector> cTS,
                      const curvatureType cType,
                      std::vector<std::pair<double,double> > &minMaxCurvatures,
                      const int shapeType,
                      osg::ref_ptr<osg::MatrixTransform> addToThisGroup );
};

#endif // SCURVATUREANALYSIS_HH
