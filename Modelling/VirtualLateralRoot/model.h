#ifndef Model_HH
#define Model_HH

#include "ModelHeader.h"
#include "ModelUtils.h"
#include "InitSurface.h"
#include "DivisionSetting.h"

class MyModel : public Model 
{
  Q_OBJECT
public: 
  
  MyModel(QObject *parent);
  
  ~MyModel();
  
  void readParms();
    
  void modifiedFiles( const std::set<std::string>& filenames );
    
  void renderImage();
  
  void restartModel();
  
  void setMovieParameters();
  
  double determineTotalArea();
  
  void setStatus();
  
  void updateFromOld( const cell& cl, const cell& cr, const cell& c,
                      const MyTissue::division_data& ddata, MyTissue& );
  
  void checkLayerAppearance( const cell &c );
  
  void setCellProperties( const cell &c, const cell &parentCell );
  
  void setLayerValues( const cell& cl, const cell& cr, const cell& c,
                       const DivisionType::type divType, const bool updateBothLayers );
  
  void setCellFileSequence( const cell& cl, const cell& cr,
                            const cell& c, const DivisionType::type divType, const bool updateBothFiles );
  
  void step();
  
  void step_divisions();

  void applyDivisionRule( const cell &c, const std::size_t type );
  
  void applyDivision( const cell &c, MyTissue::division_data &ddata );
  
  void updateLayerCount();
  
  void initDraw( Viewer* viewer );
  
  void preDraw();
  
  void postDraw()
  {T.postDraw();}
  
  void draw( Viewer* viewer );
  
  Colorf cellColor(const cell& c);
  
  // color for contour of cells
  Colorf contourColor(const cell& c)
  {return _palette.getColor(T.contourColor);}

  // Method needed by the tissue
  Point3d position(const cell& c) const
  {return c->getPos();}
  
  // Method needed by the tissue
  Point3d position(const junction& c) const
  {return c->getPos();}
  
  void setPosition(const cell& c, const Point3d& p);
  
  void setPosition(const junction& j, const Point3d& p);
  
  // Method needed by the tissue
  void setPositionHint(const junction&, const junction&, const junction&, double){}
  
  // Method needed by the tissue
  Point3d normal(const junction& ) const
  { return Point3d(0,0,1); }
  
  // Method needed by the tissue
  Point3d normal(const cell& ) const
  { return Point3d(0,0,1);}
  
  void exportDataProperties( const std::string &filename );
  
  void readDataProperties( const std::string &filename );
  
  util::Parms _parms;
  Surface _VLRBezierSurface;
  RealSurface _VLRDataPointSurface;
  DivisionSetting *_divSetting;
  util::Palette _palette;
  MyTissue T;
#ifndef TimeAgainstCellAnalysis
  double dt;
#endif
  bool _useAreaRatio;
  double _equalAreaRatio;
  bool _useCombinedAreaRatio;
  bool _useWallRatio;
  double _divisionArea;
  double _divisionAreaRatio;
  double _divisionWallRatio;
  int bgColor;
  int stepPerView;
  std::size_t _initialCellNumber;
  std::string _realDataName;
  double _probabilityOfDecussationDivision;
  double _angleThreshold;
  bool _bezierGrowthSurface;
  double _surfaceScale;
  int _surfaceType;
  bool _useAutomaticContourPoints;
  std::string _lineageFileName;
  std::string _cellWallsFileName;
  std::string _divisionFileName;
  //std::string _timeAgainstCellsFileName;
  std::vector<std::size_t> _layerColorIndex;
  std::vector<std::size_t> _cellFileColorIndex;
  double _cellPinch;
  double _cellMaxPinch;
  std::size_t _cellColoringType;
  std::pair<std::size_t, std::size_t> _divOccurrences;
  InitSurface _initSurface;
  std::size_t _lod;
  bool _lastStep;
  std::size_t _initialSituationType;
  Point3d _initialBoundary;
  double _firstDivisionsAreaRatio;
  double _secondDivisionsAreaRatio;
  std::size_t _timeDelay;
  bool _accurateCenterOfMass;
  double _LODThreshold;
  double _avoidTrianglesThreshold;
  bool _useAlternativeDT;
  std::string _divisionType;
  layerMap _firstLayerAppearances;
  bool _useLoop;
  std::size_t _amountLoops;
  std::size_t _loopCounter;
  std::map<std::string, std::size_t> _totalLayerCount;
  bool _drawControlPoints;
  bool _drawBezierSurface;
  bool _interpolateBezierSurfaces;
  bool _drawSpheres;
  bool _drawLayerCurves;
  bool _drawCenter;
  bool _drawCellMesh;
  bool _drawPCLine;
  bool _drawCells;
  bool _drawJunctions;
  bool _loadLastModel;
  bool _onlyGrowthInHeight;
  std::size_t _highOrderPattern;
  GLUquadricObj *quadratic;
  std::size_t _maxTimeSteps;
  
  std::map< std::string, std::vector<std::size_t> > _divisionSequences;
  
  std::string _divisioAngleFileName;
  
  double _initialArea;
  double _curArea;
  double _totalArea;
  double _areaGrowthFactor;
  
  bool _renderMovies;
  unsigned int _renderMoviesIndex;
  QString _imagesFilename;
  
  std::size_t _timeFourCellStage;
  std::size_t _timeSixCellStage;
  
  std::vector<std::vector<double> > _probValues;
  std::vector<std::vector<double> > _lengths;
  std::vector<std::size_t> _choices;
  std::vector<std::pair<Point3d, Point3d> > _pcLines;
  
};

#include "model.moc"
DEFINE_MODEL(MyModel);

#endif