#ifndef SDIVISIONANALYSIS_HH
#define SDIVISIONANALYSIS_HH

/**
  @file   SDivisionAnalysis.hh
  @brief  Class for performing a division analysis
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SSimilarityMeasureGraphics.hh"

struct Angles
{
  double VLRangle;
  double primAngle;
  double surfaceAngle;
  osg::Vec3 divPos;
  std::size_t cellCycle;
  mutable std::size_t numberOfCells;
  std::size_t timeStep;
  std::size_t id;
  mutable std::size_t lineage;
  mutable int cellFile;
  mutable std::string layerInfo;
  mutable std::size_t divType;
  mutable double theta; // first angle between division orientations and PCA
  mutable double phi; // second angle between division orientations and PCA
  mutable double tau; // angle to height axis of the primordium determined by PCA from Alex
};

struct setComp
{
  bool operator()(const Angles &a1, const Angles &a2) const
  {
    return ( a1.VLRangle <= a2.VLRangle );
  }
};

class SDivisionAnalysis
{
public:
  SDivisionAnalysis( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                     const std::size_t numTotalTimeSteps,
                     const std::size_t numCellCycles,
                     const std::size_t numRegSteps,
                     const bool onlyMasterCellFile,
                     osg::ref_ptr<SSliderCallback> sliderCallback,
                     std::vector<osg::ref_ptr<osg::MatrixTransform> > addToThisGroup,
                     const std::size_t sliderType,
                     const int angleType );

  void determineDivisionsPerCellCycle();

  void registerBasedOnNumberOfCells();

  void performAngleAnalysis();

  void performCylindricalAngleAnalysis();

  void renderSurfaceContours();

//  void renderDivisions( const int angleType );

  void renderContour( const std::vector<osg::Vec3> &boundary,
                      const std::size_t step );

  void renderClosestPointAndLine( const osg::Vec3 &p1,
                                  const osg::Vec3 &p2,
                                  const std::size_t t );

  void checkMinMax( const double value, const int type );

  void checkMinMaxBoundary( const double value,
                            const std::size_t dim );

  void updateBoundaryValues( const std::vector<osg::Vec3> &boundary );

  void initRegDivisionProp();

  std::size_t getRegIndex( const std::size_t numCells );

  std::pair<double, double> getMinMaxValues()
  { return _minMaxAngleValues.at(_angleType); }

  std::vector< std::vector< std::pair<double,double> > > getMinMaxBoundary()
  { return _minMaxBoundary; }

  std::vector< std::set<Angles, setComp> > getDivisionProperties()
  { return _divisionProperties; }

  std::vector< std::set<Angles, setComp> > getRegisteredDivisionProperties();

  void readAngleData();

private:

  std::vector< std::vector<std::pair<double,double> > > _minMaxBoundary;

  std::map< std::pair< std::size_t, std::size_t >, Angles > _angleValues;

  std::vector< std::set<Angles, setComp> > _divisionProperties;

  std::map< std::size_t, std::set<Angles, setComp> > _registeredDivisionProperties;

  std::vector<int> _registeredShapes;

  std::vector< std::pair<double, double> > _minMaxAngleValues;

  std::vector< cellSet > _cellsPerTimeStep;

  std::vector<osg::Vec3> _centresOfMass;

  std::size_t _numTotalTimeSteps;

  std::size_t _numCellCycles;

  std::size_t _numRegSteps;

  std::vector<int> _averagedCycleShapeTimes;

  std::size_t _sliderType;

  boost::shared_ptr<SAdvancedLineageTreeCollection> _lineages;

  osg::Matrix _rotMatrix;
  osg::Matrix _inverseRotMatrix;

  bool _onlyMasterCellFile;

  osg::ref_ptr<osg::Group> _renderGroup;

  std::vector< std::vector<osg::ref_ptr<osg::Group> > > _surfaceGroup;

  std::vector<SColorMap*> _colorMaps;

  int _angleType;

  std::size_t _currentView;

  std::vector<boost::shared_ptr<SAlphaShape> > _alphaShapes;
};

#endif // SDIVISIONANALYSIS_HH
