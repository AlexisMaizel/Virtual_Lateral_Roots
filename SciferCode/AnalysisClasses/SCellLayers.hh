#ifndef SCELLLAYERS_HH
#define SCELLLAYERS_HH

/**
  @file   SCellLayers.hh
  @brief  Contains class for generating cell layers of arabidopsis data
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SAlphaShape.hh"
#include "SConvexHull.hh"
#include "SInteractiveNormals.hh"
#include "SInteractiveSurface.hh"
#include "SSliderCallback.hh"
#include "SSimilarityMeasureHeader.hh"

#include <osg/MatrixTransform>

class SCellLayers
{
public:

  /**
  * constructor for cell layer class
  *
  * @param lineages
  * @param changeLayerThreshold
  * @param newLayerThreshold
  * @param inverseMatrix
  * @param layerColors
  * @param arrowParam1
  * @param arrowParam2
  * @param wireframe
  *
  */
  SCellLayers( const boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
               const double changeLayerThreshold,
               const osg::Matrix &inverseMatrix,
               const osg::ref_ptr<osg::Vec4Array> layerColors,
               const float arrowParam1,
               const float arrowParam2,
               const bool wireframe,
               const int shapeType );

  /**
  * constructor for cell layer class when only the cell layer information
  * is required but nothing should be drawn
  *
  * @param lineages
  * @param changeLayerThreshold
  *
  */
  SCellLayers( const boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
               const double changeLayerThreshold );

  /**
  * Initialization all requires data types for cell layers
  *
  * Is called once in determineCellsPerTimeStep()
  */
  void initializeDataTypes();

  /**
  * Initial filling of all cells per time step
  *
  * Is called once in constructor
  */
  void determineCellsPerTimeStep();

  void determineLayers( osg::ref_ptr<osg::MatrixTransform> addToGroup,
                        osg::ref_ptr<SSliderCallback> sliderCallback,
                        const bool layerInfoLoaded,
                        const bool renderDelaunay = true );

  void determineLayersWithoutGeometry( const bool layerInfoLoaded );

  void updateLayers( const unsigned int timeStep,
                     SInteractiveSurface *is,
                     SInteractiveNormals *in );

  void generateCellLayers( const unsigned int startTimeStep,
                           SInteractiveSurface *is,
                           SInteractiveNormals *in,
                           const bool update,
                           const bool layerInfoLoaded );

  void generateLayerGeometry( const unsigned int timeStep,
                              SInteractiveSurface *is,
                              const cellSet cS,
                              std::vector<unsigned int> &cellsInLayer,
                              std::vector<boost::shared_ptr<SConvexHull> > &currentDTVector,
                              cellHistoryMap &currentCellHistoryMap,
                              const bool layerInfoLoaded );

  void generateLayerGeometry( const unsigned int timeStep,
                              SInteractiveSurface *is,
                              const cellSet cS,
                              std::vector<unsigned int> &cellsInLayer,
                              std::vector<boost::shared_ptr<SAlphaShape> > &currentDTVector,
                              cellHistoryMap &currentCellHistoryMap,
                              const bool layerInfoLoaded );

  void assignNewLayer( cellHistoryMap &cellHistories,
                       const SLineageTree *cell,
                       const int childChoice,
                       const int prevLayer );

  bool assignNewLayerValue( const unsigned int layerValue,
                            const int cellId,
                            const std::size_t timeStep,
                            unsigned int &oldLayerValue );

  int changeLayerCheck( const SLineageTree *parent,
                        const osg::Vec3 &pNormal );

  void updateCellHistory( cellHistoryMap &cellHistories,
                          const SLineageTree *cell,
                          const unsigned int newLayer );

  void checkLayerSizes( const unsigned int timeStep );

  osg::Vec3 getCenterNormal( const osg::Vec3 &vPos );

  void clearLayerInformation( const unsigned int startTimeStep );

  unsigned int determineLayerValue( const SLineageTree *cell );

  unsigned int determineLayerValueFromPair( const std::pair<int, int> &IDandT );

  osg::Vec3 getNearestNeighborNormal( const std::vector<boost::shared_ptr<SConvexHull> > &currentDTVector,
                                      const unsigned int layerValue,
                                      const unsigned int timeStep,
                                      const osg::Vec3 &cellPos );

  osg::Vec3 getNearestNeighborNormal( const std::vector<boost::shared_ptr<SAlphaShape> > &currentDTVector,
                                      const unsigned int layerValue,
                                      const unsigned int timeStep,
                                      const osg::Vec3 &cellPos );

  unsigned int getMaxNumLayers() const
  { return _maxNumLayers; }

  void setMaxNumLayers( const unsigned int maxNumberLayers )
  { _maxNumLayers = maxNumberLayers; }

  unsigned int getCurrentNumLayers() const
  { return _currentNumLayers; }

  boost::shared_ptr<cellLayerVector> getCellLayersPerTimeStep() const
  { return _cellLayersPerTimeStep; }

  boost::shared_ptr<NodeFeatureInfo> getCellsPerTimeStepInLayer() const
  { return _cellsPerTimeStepInLayer; }

  boost::shared_ptr<cellTimeVector> getCellsPerTimeStep() const
  { return _cellsPerTimeStep; }

  boost::shared_ptr< std::vector<cellHistoryMap> > getCellHistories() const
  { return _cellHistories; }

  boost::shared_ptr< std::map<const SLineageTree*, osg::Vec3> > getNormalDivisionMap() const
  { return _normalDivisionMap; }

  void printTimeComplexities() const
  {
    std::cout << "Geometry generation: " << _timerForGeometryGeneration << std::endl;
    std::cout << "Division Type generation: " << _timerForLayerGeneration << std::endl;
  }

  void transformLayerData();

private:

  /// current number of layers during generation
  unsigned int _currentNumLayers;

  /// maximal number of layers after complete
  /// automatic layer generation
  unsigned int _maxNumLayers;

  /// The loaded advanced lineage tree collection
  boost::shared_ptr<SAdvancedLineageTreeCollection> _lineages;

  /// This vector stores all lineage ids which are to be displayed.
  std::vector<int> _displayedLineageIds;

  /// threshold of a cell for changing its layer or not
  /// in other words: anticlinal/radial or periclinal division
  double _changeLayerThreshold;

  /// inverse matrix of rotation applied to data set
  osg::Matrix _inverseMatrix;

  /// unique colors for the cells
  osg::ref_ptr<osg::Vec4Array> _layerColors;

  /// parameters for arrows or spheres illustrating divisions
  float _arrowParam1;
  float _arrowParam2;

  /// render layer surface in wireframe mode
  bool _wireframe;

  /// this variable is set to true if only the layer info
  /// is required but not the geometry nodes (e.g. delaunay etc. )
  bool _onlyLayerInfo;

  /// render convex hull or alpha shape
  /// 0 -> convex hull
  /// 1 -> alpha shape
  int _shapeType;

  /// data structure to store for each time step and cell id
  /// the corresponding cell layer value
  boost::shared_ptr<NodeFeatureInfo> _cellsPerTimeStepInLayer;

  /// per time step, per layer, per cell
  /// will be updated and generated every time step
  boost::shared_ptr<cellLayerVector> _cellLayersPerTimeStep;

  /// filled and initialized at the beginning with cells
  /// per time step
  /// always const and will not be updated
  boost::shared_ptr<cellTimeVector> _cellsPerTimeStep;

  /// history of a cell in which layers this cell
  /// was in the past
  /// will be updated and generated every time step
  boost::shared_ptr< std::vector<cellHistoryMap> > _cellHistories;

  /// map for storing all normals of a division node
  boost::shared_ptr< std::map<const SLineageTree*, osg::Vec3> > _normalDivisionMap;

  /// storage of volume values for each time step and layer
  std::vector< std::vector<double> > _layerVolumes;

  double _timerForGeometryGeneration;
  double _timerForLayerGeneration;
};

#endif // SCELLLAYERS_HH
