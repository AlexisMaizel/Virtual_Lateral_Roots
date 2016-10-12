#ifndef SCELLLAYERVALIDATION_HH
#define SCELLLAYERVALIDATION_HH

/**
  @file   SCellLayerValidation.hh
  @brief  Contains class for validating layer results
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SClippingPlane.hh"
#include "SMIPValidation.hh"

#include "guiQt/SCellLayerWindow.hh"

class SCellLayerValidation
{
public:
  SCellLayerValidation( const std::string &MIPFilename,
                        const int projDir,
                        const bool useSlicing,
                        boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                        osg::ref_ptr<SSliderCallback> sliderCallback,
                        osg::ref_ptr<osg::MatrixTransform> addToThisGroup );

  ~SCellLayerValidation();

  void initializeClippingPlane( const int projDir );

  void renderClippingPlane( const int projDir );

  void renderMIPs( const int projDir );

  void update( const int index );

  void removeProjections();

  void generateBBox();

public:

  boost::shared_ptr<SAdvancedLineageTreeCollection> _lineages;

  /// Direction of the MIP projection
  int _projDir;

  /// filename to MIP raw data
  std::string _MIPFilename;

  /// use slicing
  bool _useSlicing;

  osg::ref_ptr<osg::MatrixTransform> _vlrGroup;

  boost::shared_ptr<SMIPValidation> _MIPValidation;

  /// bounding box information
  osg::BoundingBox *_bbox;
};

#endif // SCELLLAYERVALIDATION_HH
