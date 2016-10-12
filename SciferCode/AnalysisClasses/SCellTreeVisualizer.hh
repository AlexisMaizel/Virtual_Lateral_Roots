#ifndef SCELLTREEVISUALIZER_HH
#define SCELLTREEVISUALIZER_HH

#include "SSimilarityMeasureGraphics.hh"

#include "dataSet/SAdvancedLineageTreeCollection.hh"

#include "kernel/SKernel.hh"

class SCellTreeVisualizer
{
public:
  SCellTreeVisualizer( const std::vector< boost::shared_ptr<SAdvancedLineageTreeCollection> > &lineages,
                       const std::vector<std::size_t> &numTotalTimesteps,
                       const std::vector<NodeFeatureInfo> &divisionScheme,
                       const std::vector<NodeFeatureInfo> &layerValues,
                       const std::vector< std::map<int,int> > &cellFiles,
                       const std::vector<osg::Matrix> &rotMatrices,
                       const std::vector< std::vector<osg::Vec3> > &centresOfMass,
                       const osg::ref_ptr<osg::Vec4Array> layerColors,
                       const bool onlyMasterFiles,
                       const int renderLineageLineType,
                       const int renderLineageNodeType );

  void generateLineageTrees( const bool registerTrees,
                             const std::size_t divStop,
                             osg::ref_ptr<osg::Group> addToThisGroup,
                             const SNodeInfo &info,
                             const int windowId );

  void generateTreeDivisionSequence( const std::size_t sequenceAxisType,
                                     const bool loadModel,
                                     const SNodeInfo &info,
                                     const int windowId );

  void generateAveragedTreeDivisionSequence( const std::string &modelDataDirectory,
                                             const SNodeInfo &info,
                                             const int windowId );

private:

  std::vector< boost::shared_ptr<SAdvancedLineageTreeCollection> > _lineages;
  std::vector<std::size_t> _numTotalTimesteps;
  std::vector< std::map<int,int> > _cellFiles;
  std::vector<NodeFeatureInfo> _divisionScheme;
  std::vector<NodeFeatureInfo> _layerValues;
  std::vector<osg::Matrix> _rotMatrices;
  std::vector< std::vector<osg::Vec3> > _centresOfMass;
  osg::ref_ptr<osg::Vec4Array> _layerColors;
  bool _onlyMasterFiles;
  int _renderLineageLineType;
  int _renderLineageNodeType;
};

#endif // SCELLTREEVISUALIZER_HH
