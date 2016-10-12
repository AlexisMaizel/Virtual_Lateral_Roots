#ifndef SVLRBrowser_hh__
#define SVLRBrowser_hh__

/**
  @file   SVLRBrowser.hh
  @brief  Contains class for creating a feature browser of VLR data
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "../similarityMeasure/SCellLayers.hh"
#include "../similarityMeasure/SInteractiveCell.hh"
#include "../similarityMeasure/SGeometryNodeMaps.hh"
#include "../similarityMeasure/SLineageTreeLayers.hh"
#include "../similarityMeasure/SSliderCallback.hh"
#include "../similarityMeasure/SLineageCallback.hh"
#include "../similarityMeasure/SVoronoiDiagram.hh"
#include "../similarityMeasure/SCellLayerValidation.hh"

#include "guiQt/SCellLayerWindow.hh"
#include "guiQt/SCellLayerWindowCreator.hh"

#include "kernel/SAlgorithm.hh"
#include "kernel/SDatasetFilter.hh"
#include "kernel/SNodeInfo.hh"

#include "graphics/SSelectableEllipse.hh"
#include "graphics/SSelectableSphere.hh"

#include <boost/signals2/connection.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>

#include <osg/MatrixTransform>

class SVLRBrowser : public SAlgorithm
{
public:
  SVLRBrowser();
  ~SVLRBrowser();
  bool canUndo();
  void undo();
  std::string getMenuEntry() const;
  std::string getName() const;
  static std::string algoName();
  void execute( void* );

  void generateTreeDivisionSequence();

  void generateAveragedTreeDivisionSequence();

  void generateLineageTrees();

  void generateVoronoiCells();

  void lengthAnalysis();

  void performDivisionAnalysis();

  void performCorrelationAnalysis();

  void performDivisionSequenceAnalysis();

  void performDivisionAngleAnalysis();

  void curvatureAnalysis();

  void developmentAnalysis();

  void renderSingleTimeSteps();

  void generateDivisionArrows();

  void selectNode( osg::Node *node,
                   const bool select,
                   const NodeType::NodeType nT );

  void highlightCell( osg::Node* node, const NodeType::NodeType nT );

  void assignLayer( osg::Node* node, const int layer,
                    const NodeType::NodeType nT );

  void updateLayerGroups( osg::Node* node,
                          const int newLayer,
                          const int oldLayer );

  void setDivisionType( osg::Node* node, const int divisionType,
                        const NodeType::NodeType nT );

  void bindSlider( const std::string &sliderName );

  void updateLayerColorsForCells();

  void saveLayers();

  void apply();

  void generateFrames( osg::ref_ptr<osg::Group> renderGroup, const bool exportAsSVG );

  void setCellFeatures( const int cellId,
                        const int timeStep,
                        const int treeId,
                        const int layer,
                        const cellHistory cellHist,
                        const double xValue,
                        const double yValue,
                        const double zValue,
                        const int cellFile,
                        const std::size_t cellDivisionType );

  SLineageTree* getTreeIdOfSelectedCell( int &dataid, int &layer );

  void createLayerWindow();

  void initializeSlider();

  void exportCellDivisions( const std::string &fileName );

  double computeAngle( const SLineageTree *cell );

  void setColorMaps();

  void loadAndGenerateCellShapes();

  boost::signals2::connection onHighlightCell( const slotHighlight& slot );
  boost::signals2::connection onLayerAssignment( const slotLayerAssignment& slot );

  typedef boost::signals2::signal<void (std::vector<int>)> signalHighlightUpdate;
  typedef signalHighlightUpdate::slot_type slotHighlightUpdate;
  boost::signals2::connection onHighlightChanged( const slotHighlightUpdate& slot );

protected:
  void initProfile();

  /// Triggers the appropriate changes in visualization if _tfoc has changed
  void onSliderUpdate( int val );

private:
  /// Filter that accepts only SAdvancedLineageTreeCollection datasets
  SAcceptDatasetsOfClassFilter<SAdvancedLineageTreeCollection> _lineageFilter;

  /// Lineage dataset info
  std::vector<int> _lineageDataId;

  /// window id for the 2D lineage visualization
  int _lineageWindowId;

  /// node info for the 2D lineage visualization
  SNodeInfo _lineageInfo;

  /// node info for colorbar
  SNodeInfo _colorbarInfo;

  /// window id for the 2D division visualization
  int _divisionWindowId;

  /// node info for the 2D division visualization
  SNodeInfo _divisionInfo;

  /// window id for the 2D averaged division visualization
  int _averagedDivisionWindowId;

  /// node info for the 2D averaged division visualization
  SNodeInfo _averagedDivisionInfo;

  int _correlationAnalysisWindowId;

  SNodeInfo _correlationAnalysisInfo;

  /// id of window for precise analysis of 2D division information
  int _divisionAnalysisWindowId;

  /// precise analysis of 2D division information
  SNodeInfo _divisionAnalysisInfo;

  /// pointer to colorbar instance
  SColorMap *_colorbar;

  /// use of multiple lineage data
  std::size_t _numData;

  SNodeInfo _infoBar;

  /**
   * Variables concernig the visualization of the temporal development.
   * _tmin and _tmax refer to minimum and maximum time step found in the data.
   * _tfoc is the time step that the user's interest is focused on.
   * _twin is the time window specified by the user that (from _tfoc onwards) is
   * shown with all details (e.g. full opacity).
   * So we distinguish between the following three time intervals:
   * [_tmin,_tfoc) -> past; perhaps increasing opacity or completely pruned
   * [_tfoc,_tfoc+_twin) -> present and near future; all details
   * [_tfoc+_twin,_tmax] -> far future; increasing transparency
   */
  int _tmin, _tmax, _tfoc, _twin, _ccmax;

  /// vector of unique names for the slider
  std::vector<std::string> _sliderLabels;

  bool _firstClick;

  /// selected data id for cell layer visualization
  unsigned int _sDataId;

  /// The loaded advanced lineage tree collection
  std::vector< boost::shared_ptr<SAdvancedLineageTreeCollection> > _lineages;

  std::vector<std::size_t> _numTotalTimesteps;

  std::vector<std::size_t> _numCellCycles;

  std::vector< std::vector< std::set<osg::Vec3> > > _divisionsPerCellCycle;

  /// node infor for 3D cells
  std::vector<SNodeInfo> _vlrInfo;

  /// window id for each data set
  std::vector<int> _vlrWindowId;

  /// whole group information stored in a matrix transformation
  /// to apply initial rotation information
  osg::ref_ptr<osg::MatrixTransform> _vlrRoot;

  /// group consisting of all lineages
  osg::ref_ptr<osg::Group> _lineageLayerGroup;

  /// cell layer information stored in a matrix transformation
  /// to apply inital rotation information
  osg::ref_ptr<osg::MatrixTransform> _rootLayers;

  /// directory to voronoi data
  std::string _voronoiDirectory;

  /// group for division analysis
  osg::ref_ptr<osg::Group> _divGroup;

  std::vector<osg::Matrix> _rotMatrices;

  std::vector<osg::Matrix> _rotInverseMatrices;

  /// draw labels
  bool _divisionLabels;

  /// render only master cell files
  bool _onlyMasterFiles;

  std::size_t _gridResolution;

  /// render lineage trees until divStop-th division
  std::size_t _divStop;

  /// threshold of a cell for changing its layer or not
  /// in other words: anticlinal/radial or periclinal division
  double _changeLayerThreshold;

  /// threshold for checking if the current division is an anticlinal
  /// or radial division
  double _radialAngleThreshold;

  /// render layer surface in wireframe mode
  bool _wireframe;

  /// enable slicing of layers and cells
  bool _enableSlicing;

  /// draw 3D cell visualization
  bool _draw3DCells;

  /// render lineages or not
  bool _drawLineages;

  /// perform curvature analysis
  bool _renderCurvature;

  /// curvature type
  /// 0 -> mean curvature
  /// 1 -> gaussian curvature
  int _curvatureType;

  /// string of curvature types
  std::list<std::string> _curvatureTypes;

  /// colors for layers
  osg::ref_ptr<osg::Vec4Array> _layerColors;

  /// colors for cell files
  std::map<int, osg::Vec4> _cellFileColors;

  /// colors for cell lineages
  std::map<int, osg::Vec4> _cellLineageColors;

  std::map<std::size_t, std::size_t> _treeIdMap;

  /// colors for cell lineages
  SColorMap _areaColorMap;

  /// colors for life duration
  SColorMap _lifeDurationColorMap;

  /// cell layers instance
  boost::shared_ptr<SCellLayers> _cellLayers;

  /// class instance storing all maps from an lineage tree node
  /// to a selectable ellipse or sphere and vice versa
  boost::shared_ptr<SGeometryNodeMaps> _geomNodeMaps;

  signalHighlight _signalHighlightCell;
  signalLayerAssignment _signalLayerAssignment;
  signalHighlightUpdate _signalHighlightUpdate;

  /// filename of MIP raw images
  std::string _MIPfilename;

  /// Direction of the raw projection
  int _projDir;

  /// string types of projections: x, y, z, all
  std::list<std::string> _projections;

  boost::shared_ptr<SCellLayerValidation> _validation;

  /// slider type for either the cell cycles or time steps
  /// 0 -> time slider
  /// 1 -> cycle slider
  /// 2 -> registered step slider
  int _sliderType;
  std::list<std::string> _sliderTypes;

  /// shape type of primordium
  /// 0 -> convex hull
  /// 1 -> alpha shape
  int _shapeType;
  std::list<std::string> _shapeTypes;

  /// lineage color type
  /// 0 -> layer and division type
  /// 1 -> area value for tissue of cell
  int _lineageColorType;
  std::list<std::string> _lineageColorTypes;

  /// axis type for periclinal division sequences
  /// 0 -> time steps
  /// 1 -> cells
  int _sequenceAxisType;
  std::list<std::string> _sequenceAxisTypes;

  /// angle types for division analysis
  /// 0 -> Centre of mass
  /// 1 -> Primordium centre
  /// 2 -> Surface normal
  int _angleTypeDivAnalysis;
  std::list<std::string> _angleTypesDivAnalysis;

  /// render type of 3D VLR
  /// 0 -> Layering: spheres are colored based on layers plus delaunay triangulation
  /// 1 -> Cell Files: spheres are colored based on cell files assignment
  int _render3DType;
  std::list<std::string> _render3DTypes;

  int _renderLineageLineType;
  std::list<std::string> _renderLineageLineTypes;

  int _renderLineageNodeType;
  std::list<std::string> _renderLineageNodeTypes;

  /// color type of divisions
  /// 0 -> Type
  /// 1 -> Layer
  int _divisionColorType;
  std::list<std::string> _divisionColorTypes;

  std::vector< std::pair<SInteractiveCell*,const SLineageTree*> > _interactiveCellVector;

  osg::Vec4 _selectionColor;
  osg::Vec4 _lastClickedColor;
  SSelectableSphere *_lastClicked3DNode;
  SSelectableEllipse *_lastClicked2DNode;

  int _clickedCellId;
  std::size_t _clickedTimeStep;

  std::size_t _regSteps;

  /// radius of cell nodes
  double _cellRadius;

  /// render lineage tracks or not
  bool _renderTracks;

  /// load layer information from file or not
  bool _loadLayerInformation;

  /// load division type information from file or not
  bool _loadDivisionTypeInformation;

  bool _renderArrows;

  /// pointer to layer window creator handling the communication
  /// between the layer window widget and this algorithm
  SCellLayerWindowCreator* _layerWindowCreator;

  /// pointer to layer window which enables the access to specific
  /// methods of the widget
  SCellLayerWindow* _layerWindow;

  /// status if the signals are already connected or not
  bool _signalsConnected;

  /// voronoi diagram and rendering of cells
  boost::shared_ptr<SVoronoiDiagram> _voronoi;

  /// slider callback for interactive rendering
  osg::ref_ptr<SSliderCallback> _sliderCallback;

  /// triangulation of complete primordium or single layers
  bool _renderCompletePrimordium;

  bool _renderCorrelationAnalysis;

  /// cell file information
  std::vector< std::map<int,int> > _cellFiles;

  /// division scheme information
  std::vector<NodeFeatureInfo> _divisionScheme;

  /// layer value information
  std::vector<NodeFeatureInfo> _layerValues;

  /// max layer value for all loaded datasets
  unsigned int _maxLayers;

  /// cells per time step for each dataset
  std::vector< boost::shared_ptr<cellTimeVector> > _cellsPerTimeStep;

  std::vector< std::vector<osg::Vec3> > _centresOfMass;

  bool _renderTreeDivisionSequence;

  /// register trees based on the number of cells
  bool _registerTrees;

  bool _renderDivisionAnalysis;

  /// render layer surfaces
  bool _renderLayerSurfaces;

  /// load model data or not; if so then some features are disabled
  bool _loadModel;

  /// render cell wall which is atm only applied when model data is loaded
  bool _renderCellWalls;

  /// width of cell walls
  double _cellWallWidth;

  /// path to file with cell wall information
  std::string _cellWallDirectory;

  /// path to file in which model data properties are stored for several sample runs
  std::string _modelDataDirectory;
};

#endif
