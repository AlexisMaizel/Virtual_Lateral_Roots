#include "ModelHeader.h"
#include "ModelExporter.h"
#include "ModelUtils.h"
#include "GraphicsClass.h"
#include "InitSurface.h"
#include "DivisionSetting.h"

//#define TimeAgainstCellAnalysis

const double start = 0.05;
const double steps = 0.0005;
#ifdef TimeAgainstCellAnalysis
static double dt = start;
#endif

class MyModel : public Model 
{
  Q_OBJECT
public: 
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
  bool _drawPCLine;
  bool _drawCells;
  bool _drawJunctions;
  bool _loadLastModel;
  bool _onlyGrowthInHeight;
  std::size_t _highOrderPattern;
  GLUquadricObj *quadratic;
  std::vector<std::size_t> _randomChoices;
  std::size_t _choiceCounter;
  std::size_t _maxTimeSteps;
  
  std::map< std::string, std::vector<std::size_t> > _divisionSequences;
  
  double _equalAreaRatio;
  
  bool _renderMovies;
  unsigned int _renderMoviesIndex;
  QString _imagesFilename;
  
  std::vector<std::vector<double> > _probValues;
  std::vector<std::vector<double> > _lengths;
  std::vector<std::size_t> _choices;
  std::vector<std::pair<Point3d, Point3d> > _pcLines;
  
  //----------------------------------------------------------------
  
  void readParms()
  {
    // read the parameters here
#ifndef TimeAgainstCellAnalysis
    _parms("Main", "Dt", dt);
#endif
    _parms("Main", "InitialCellNumber", _initialCellNumber);
    _parms("Main", "RealDataName", _realDataName);
    _parms("Main", "SubDivisionLevelOfCells", _lod);
    _parms("Main", "BezierGrowthSurface", _bezierGrowthSurface);
    _parms("Main", "SurfaceScale", _surfaceScale);
    _parms("Main", "SurfaceType", _surfaceType);
    _parms("Main", "UseAutomaticContourPoints", _useAutomaticContourPoints );
    _parms("Main", "InitialSituationType", _initialSituationType );
    _parms("Main", "CenterOfMassBasedOnTriangleFan", _accurateCenterOfMass );
    _parms("Main", "Loop", _useLoop );
    _parms("Main", "AvoidTrianglesThreshold", _avoidTrianglesThreshold );
    _parms("Main", "LoadLastModel", _loadLastModel );
    _parms("Main", "OnlyGrowthInHeight", _onlyGrowthInHeight );
    _parms("Main", "RenderMovies", _renderMovies );
    _parms("Main", "HighOrderPattern", _highOrderPattern );
    
    _parms("View", "StepPerView", stepPerView);
    _parms("View", "BackgroundColor", bgColor);
    _parms("View", "RenderCells", _drawCells );
    _parms("View", "RenderJunctions", _drawJunctions );
    _parms("View", "RenderSpheres", _drawSpheres );
    _parms("View", "RenderLayerCurves", _drawLayerCurves );
    _parms("View", "RenderControlPoints", _drawControlPoints );
    _parms("View", "RenderBezierSurface", _drawBezierSurface );
    _parms("View", "RenderCellCenter", _drawCenter );
    _parms("View", "RenderPCLine", _drawPCLine );

    _parms( "Division", "DivisionArea", _divisionArea);
    _parms( "Division", "DivisionAreaRatio", _divisionAreaRatio);
    _parms( "Division", "EqualAreaRatio", _equalAreaRatio);
    _parms( "Division", "UseAreaRatio", _useAreaRatio);
    _parms( "Division", "UseCombinedAreaRatio", _useCombinedAreaRatio);
    _parms( "Division", "UseWallRatio", _useWallRatio);
    _parms( "Division", "DivisionWallRatio", _divisionWallRatio);
    _parms( "Division", "UseAlternativeDivisionType", _useAlternativeDT );
    _parms( "Division", "DivisionType", _divisionType );
    _parms( "Division", "ProbabilityOfDecussationDivision", _probabilityOfDecussationDivision );
    _parms( "Division", "DivisionAngleThreshold", _angleThreshold );
    _parms( "Division", "CellColoringType", _cellColoringType );
    
    _parms( "Division", "FirstDivisionsAreaRatio", _firstDivisionsAreaRatio);
    _parms( "Division", "SecondDivisionsAreaRatio", _secondDivisionsAreaRatio);
    _parms( "Division", "TimeDelay", _timeDelay);

    _parms( "Tissue", "CellPinch", _cellPinch );
    _parms( "Tissue", "CellMaxPinch", _cellMaxPinch );
    
    T.readParms(_parms, "Tissue");
    T.readViewParms(_parms, "TissueView");
  }

  //----------------------------------------------------------------
  
  // Here, reread the files when they are modified
  void modifiedFiles( const std::set<std::string>& filenames )
  {
    forall(const std::string& fn, filenames)
    {
      if(fn == "pal.map")
        _palette.reread();
      else if(fn == "view.v")
        readParms();
    }
  }

  //----------------------------------------------------------------
  
  MyModel(QObject *parent) : Model(parent), _parms("view.v"),
    _VLRDataPointSurface( _parms, "Surface" ),
    _palette("pal.map"), T(_palette, this),
    _divOccurrences( std::make_pair( 0, 0 ) ),
    _lastStep( false ),
    _loopCounter( 1 ),
    _amountLoops( 100 ),
    _interpolateBezierSurfaces( true ),
    _choiceCounter( 0 ),
    _renderMoviesIndex( 0 ),
    _imagesFilename( "/tmp/images/m-" )
  {
    quadratic = gluNewQuadric();
    gluQuadricNormals( quadratic, GLU_SMOOTH );
    
    readParms();
    // Registering the configuration files
    registerFile("pal.map");
    registerFile("view.v");
    
    // update LOD threshold depending on subdivision level
    _LODThreshold = 30./_lod - 0.1;
    
    _randomChoices.clear();
    
    // read the properties of the saved data model
    if( _loadLastModel )
    {
      std::string fileName = "/tmp/modelDataProp";
      ModelUtils::appendCellSituationType( fileName, _initialSituationType );
      if( _highOrderPattern != 0 )
        fileName += "_HOP" + std::to_string( _highOrderPattern );
      fileName += ".csv";
      this->readDataProperties( fileName );
    }
    
    // if movies should be rendered, create in total 9 movies for NH, 2DC, and 2DCBase
    // for a visualization by cells, spheres, layercurves
    if( _renderMovies )
    {
      if( !_loadLastModel )
      {
        std::cerr << "Movies can only be rendered if the data was stored before!" << std::endl;
        return;
      }
        
      // start with NH
      _initialSituationType = 0;
      _drawCells = true;
      _drawSpheres = false;
      _drawLayerCurves = false;
    }
    
    _maxTimeSteps = 300;
    if( _realDataName == "130508_raw" )
      _maxTimeSteps = 350;
    else if( _bezierGrowthSurface && _realDataName == "none" )
      _maxTimeSteps = 500;
    
    // check if correct surface type was chosen
    if( _bezierGrowthSurface && SURFACETYPE == 1 )
    {
      std::cout << "Wrong surface type chosen!!!" << std::endl;
      return;
    }
        
    _initSurface.init( _lod, _lineageFileName,
                       _initialSituationType != 0,
                       _useAutomaticContourPoints,
                       _surfaceScale,
                       _realDataName );
    
    // set name strings
    _lineageFileName = "/tmp/model";
    _cellWallsFileName = "/tmp/modelCellWalls";
    _divisionFileName = "/tmp/divisionPropertiesModel";
    
    if( SURFACETYPE == 1 )
    {
      _lineageFileName += _realDataName;
      _cellWallsFileName += _realDataName;
      _divisionFileName += _realDataName;
    }
    
    ModelUtils::appendCellSituationType( _divisionFileName, _initialSituationType );
    ModelUtils::appendCellSituationType( _lineageFileName, _initialSituationType );
    ModelUtils::appendCellSituationType( _cellWallsFileName, _initialSituationType );
    
    if( _useAlternativeDT )
    {
      _divisionFileName += _divisionType;
      if( _divisionType != "Besson-Dumais" )
      {
        _divisionFileName += std::to_string( (unsigned int)_avoidTrianglesThreshold );
        _lineageFileName += std::to_string( (unsigned int)_avoidTrianglesThreshold );
        _cellWallsFileName += std::to_string( (unsigned int)_avoidTrianglesThreshold );
      }
    }
    else
    {
      _divisionFileName += "ShortestWall";
      _lineageFileName += "ShortestWall";
      _cellWallsFileName += "ShortestWall";
      
    }
    
    if( _highOrderPattern != 0 )
    {
      _divisionFileName += "_HOP" + std::to_string( _highOrderPattern );
      _lineageFileName += "_HOP" + std::to_string( _highOrderPattern );
      _cellWallsFileName += "_HOP" + std::to_string( _highOrderPattern );
    }
    
    _divisionFileName += ".csv";
    _lineageFileName += ".csv";
    _cellWallsFileName += ".csv";
    
    // single layer assignment
    //for( std::size_t l = 9; l < 14; l++ )
    // multiple layer assignment for each new daughter cell
    for( std::size_t l = 14; l < 45; l++ )
      _layerColorIndex.push_back( l );
    
    for( std::size_t l = 45; l < 53; l++ )
      _cellFileColorIndex.push_back( l );
    
    cell dummy;
    ModelExporter::exportLineageInformation( _lineageFileName, dummy, T, true );
    ModelExporter::exportCellWalls( _cellWallsFileName, dummy, T, true );
    
    std::pair<std::size_t, std::size_t> pair;
    ModelExporter::exportDivisionDaughterProperties( _divisionFileName,
                                                     dummy, dummy,
                                                     DivisionType::ANTICLINAL, 0.,
                                                     pair, true );
    
    /*
    _timeAgainstCellsFileName = "/tmp/timeAgainstCells";
    _timeAgainstCellsFileName += _realDataName;
    _timeAgainstCellsFileName += ".csv";
   
    if( dt > start-steps )
      ModelExporter::exportTimeAgainstCells( _timeAgainstCellsFileName, dt, 0, true );
    */
        
    // bezier
    if( SURFACETYPE == 0 )
    {
      // side
      if( _surfaceType == 0 )
      {
        _VLRBezierSurface.init( _parms, "Surface", _bezierGrowthSurface,
                          _interpolateBezierSurfaces,
                          _onlyGrowthInHeight, _realDataName,
                          _highOrderPattern );
        
        _VLRBezierSurface.growStep( 0 );
        
        // special cases of number of cells at the beginning
        if( _realDataName == "none" )
          _initSurface.initIdealizedCells( T, _initialCellNumber,
                                          _VLRBezierSurface );
        // else constellation of founder cells according to real data
        else
          _initSurface.initRealDataCells( T, _VLRBezierSurface );
      }
      // radial
      else
      {
        _VLRBezierSurface.initRadialSurface( _parms, "Surface" );
        _VLRBezierSurface.growStep( 0 );
        _initSurface.initRadialCells( T, _VLRBezierSurface );
      }
    }
    // real data points
    else
    {
      _VLRDataPointSurface.init( _surfaceScale,
                                 _realDataName,
                                 _useAutomaticContourPoints );
      
      std::vector<SurfacePoint> sps;
      _VLRDataPointSurface.growStep( 0, sps );
        
      _initSurface.initRealDataCells( T, _VLRDataPointSurface );
    }
    
    _divSetting = new DivisionSetting( T, _initialSituationType,
      _divisionType, _timeDelay, _firstDivisionsAreaRatio,
      _secondDivisionsAreaRatio, _useAlternativeDT, 
      _accurateCenterOfMass, _probabilityOfDecussationDivision,
      _useAreaRatio, _useCombinedAreaRatio, _useWallRatio,
      _divisionArea, _divisionAreaRatio, _divisionWallRatio,
      _LODThreshold, _avoidTrianglesThreshold, _loadLastModel,
      _cellPinch, _cellMaxPinch, _onlyGrowthInHeight,
      _VLRBezierSurface, _VLRDataPointSurface );
    
    // export initial cell constellation
    forall(const cell& c, T.C)
    {
      ModelExporter::exportLineageInformation( _lineageFileName, c, T, false );
      ModelExporter::exportCellWalls( _cellWallsFileName, c, T, false );
      
      this->checkLayerAppearance( c );
    }
    
    // determine the initial boundary between two founder cells
    _initialBoundary = Point3d( 0., 0., 0. );
    if( _initialSituationType != 0 )
    {
      Point3d max = Point3d( -5000., -5000., 0. );
      // consider the right most position
      forall(const cell& c, T.C)
      {
        if( c->id == 1 )
        {
          forall(const junction& j, T.S.neighbors(c))
          {
            Point3d pos = j->getPos();
            
            if( pos.i() > max.i() )
              max.i() = pos.i();
            
            if( pos.j() > max.j() )
              max.j() = pos.j();
          }
          
          _initialBoundary = max;
        }
      }
    }
    
    setStatus();
    this->renderImage();
  }

  //----------------------------------------------------------------
  
  ~MyModel()
  {
    delete _divSetting;
  }
  
  //----------------------------------------------------------------
  
  void renderImage()
  {
    if( _renderMovies )
    {
      unsigned int step = _initSurface.getTime()-1;
      QString file = _imagesFilename;
      
      // handle hard wired situations
      switch( _renderMoviesIndex )
      {
        case 0:
        case 1:
        case 2:
          file += "NH_";
        break;
        case 3:
        case 4:
        case 5:
          file += "1DC_";
        break;
        case 6:
        case 7:
        case 8:
          file += "2DC_";
        break;
        case 9:
        case 10:
        case 11:
          file += "2DCBase_";
        break;
      }
      
      if( _renderMoviesIndex%3 == 0 )
        file += "C";
      else if( _renderMoviesIndex%3 == 1 )
        file += "S";
      else
        file += "L";
      
      if( step < 10 )
        file += "000";
      else if( step < 100 )
        file += "00";
      else
        file += "0";
      
      file += QString( std::to_string(step).c_str() );
      file += ".jpg";
      this->screenshot( file, true );
    }  
  }
  
  //----------------------------------------------------------------
  
  void restartModel()
  {
    // reset the tissue to delete all previous cell information
    T = MyTissue( _palette, this );
    
    readParms();
    
    this->setMovieParameters();

    // Registering the configuration files
    registerFile("pal.map");
    registerFile("view.v");
    
    // reset division and layering results
    _divOccurrences = std::make_pair( 0, 0 );
    _firstLayerAppearances.clear();
    _totalLayerCount.clear();
    _randomChoices.clear();
    _choiceCounter = 0;
    _lastStep = false;
    
    // read the properties of the saved data model
    if( _loadLastModel )
    {
      std::string fileName = "/tmp/modelDataProp";
      ModelUtils::appendCellSituationType( fileName, _initialSituationType );
      if( _highOrderPattern != 0 )
        fileName += "_HOP" + std::to_string( _highOrderPattern );
      
      fileName += ".csv";
      this->readDataProperties( fileName );
    }
    
    _initSurface.init( _lod, _lineageFileName,
                       _initialSituationType != 0,
                       _useAutomaticContourPoints,
                       _surfaceScale,
                       _realDataName );
    
    // after each restart, reset the model and cell wall file because
    // it makes no sense to append these files
    cell dummy;
    ModelExporter::exportLineageInformation( _lineageFileName, dummy, T, true );
    ModelExporter::exportCellWalls( _cellWallsFileName, dummy, T, true );
    
    // bezier
    if( SURFACETYPE == 0 )
    {
      _VLRBezierSurface.init( _parms, "Surface", _bezierGrowthSurface,
                              _interpolateBezierSurfaces,
                              _onlyGrowthInHeight, _realDataName,
                              _highOrderPattern );
      
      _VLRBezierSurface.growStep( 0 );
      
      // special cases of number of cells at the beginning
      if( _realDataName == "none" )
        _initSurface.initIdealizedCells( T, _initialCellNumber,
                                         _VLRBezierSurface );
      // else constellation of founder cells according to real data
      else
        _initSurface.initRealDataCells( T, _VLRBezierSurface );
    }
    // real data points
    else
    {
      _VLRDataPointSurface = RealSurface( _parms, "Surface" );
      _VLRDataPointSurface.init( _surfaceScale,
                                 _realDataName,
                                 _useAutomaticContourPoints );
      
      std::vector<SurfacePoint> sps;
      _VLRDataPointSurface.growStep( 0, sps );
        
      _initSurface.initRealDataCells( T, _VLRDataPointSurface );
    }
    
    forall(const cell& c, T.C)
      this->checkLayerAppearance( c );
    
    // determine the initial boundary between two founder cells
    _initialBoundary = Point3d( 0., 0., 0. );
    if( _initialSituationType != 0 )
    {
      Point3d max = Point3d( -5000., -5000., 0. );
      // consider the right most position
      forall(const cell& c, T.C)
      {
        if( c->id == 1 )
        {
          forall(const junction& j, T.S.neighbors(c))
          {
            Point3d pos = j->getPos();
            
            if( pos.i() > max.i() )
              max.i() = pos.i();
            
            if( pos.j() > max.j() )
              max.j() = pos.j();
          }
          
          _initialBoundary = max;
        }
      }
    }
    
    setStatus();
    this->renderImage();
  }

  //----------------------------------------------------------------
  
  void setMovieParameters()
  {
    // if movies should be created then set
    // the corresponding parameters
    if( _renderMovies )
    {
      // continue with the other settings depending on
      // the parameter of _renderMoviesIndex
      // Cell vis
      if( _renderMoviesIndex%3 == 0 )
      {
        _drawCells = true;
        _drawSpheres = false;
        _drawLayerCurves = false;
        // if the same model results should be rendered with another
        // color model for example then just set _loadLastModel to true
        // here
        _loadLastModel = false;
      }
      // sphere vis
      else if( _renderMoviesIndex%3 == 1 )
      {
        _drawCells = false;
        _drawSpheres = true;
        _drawLayerCurves = false;
        _loadLastModel = true;
      }
      // layer line vis
      else
      {
        _drawCells = false;
        _drawSpheres = false;
        _drawLayerCurves = true;
        _loadLastModel = true;
      }
      
      // handle hard wired situations
      switch( _renderMoviesIndex )
      {
        case 0:
        case 1:
        case 2:
          _initialSituationType = 0;
        break;
        case 3:
        case 4:
        case 5:
          _initialSituationType = 1;
        break;
        case 6:
        case 7:
        case 8:
          _initialSituationType = 2;
        break;
        case 9:
        case 10:
        case 11:
          _initialSituationType = 3;
        break;
      }
      
      _divSetting->setAreaRatios( _loadLastModel,
                                  _initialSituationType );
    }
  }
 
  //----------------------------------------------------------------
  
  void setStatus()
  {
    std::size_t time = _initSurface.getTime();
    
    QString status = QString( "Vertices: %1\t "
                              "Cells: %2\t").
                              arg(T.W.size()).
                              arg(T.C.size());
    
    if( SURFACETYPE == 1 )
    {
      std::size_t timeStep = _VLRDataPointSurface.getCurTimeStep();
      status += QString( "TS: %1\t" ).arg(timeStep);
    }
    
    status += QString( "MS: %1\t" ).arg(time);
    if( _surfaceType == 0 )
    {
      status += QString( "AD: %1\t" ).arg(_divOccurrences.first);
      status += QString( "PD: %1\t" ).arg(_divOccurrences.second);
    }
    else
    {
      status += QString( "RD: %1\t" ).arg(_divOccurrences.first);
      status += QString( "PD: %1\t" ).arg(_divOccurrences.second);
    }
    
    if( _useCombinedAreaRatio )
    {
      QString type = "Type: ";
      ModelUtils::appendCellSituationType( type, _initialSituationType );
      type += "\t";
      status += QString( type );
      status += QString( "Area div ratio: %1\t" ).arg(_divisionAreaRatio);
    }
    else if( _useAreaRatio )
      status += QString( "Area div ratio: %1\t" ).arg(_divisionAreaRatio);
    else
      status += QString( "Area div: %1\t" ).arg(_divisionArea);
    
    if( _useWallRatio )
      status += QString( "Wall div ratio: %1\t" ).arg(_divisionWallRatio);
    
    if( _divisionType == "Decussation" )
      status += QString( "Dec prop: %1\%\t" ).arg(_probabilityOfDecussationDivision);
    
    if( _useAlternativeDT )
      status += QString( "Div Type: %1\t" ).arg( QString::fromStdString(_divisionType) );
    else
      status += QString( "Div Type: ShortestWall\t" );
    
    if( _useLoop )
      status += QString( "Loop count: %1" ).arg(_loopCounter);
    
    setStatusMessage( status );
  }
  
  //----------------------------------------------------------------
  
  void updateFromOld( const cell& cl, const cell& cr, const cell& c,
                      const MyTissue::division_data& ddata, MyTissue& )
  {
    double angle = ModelUtils::determineDivisionAngle( ddata );
    DivisionType::type divType = ModelUtils::determineDivisionType( ddata,
                                                                    _angleThreshold );
    
    // in the radial view an anticlinal division is a radial one
    if( _surfaceType == 1 && divType == DivisionType::ANTICLINAL )
      divType = DivisionType::RADIAL;
    
    // set properties of dividing cell
    c->angle = angle;
    c->divType = divType;
    
    // set cell properties for left cell
    this->setCellProperties( cl, c );
    
    // set cell properties for right cell
    this->setCellProperties( cr, c );
    
    // for the forced situation check which of the four cells
    // are the inner ones
    if( _initialSituationType != 0 && T.C.size() > 2 && T.C.size() < 5 )
    {
      double xIC = _initialBoundary.i();
      double yIC = _initialBoundary.j();
      double xLC = cl->center.i();
      double yLC = cl->center.j();
      double xRC = cr->center.i();
      double yRC = cr->center.j();
      double distL = std::sqrt( (xIC-xLC)*(xIC-xLC) + (yIC-yLC)*(yIC-yLC) );
      double distR = std::sqrt( (xIC-xRC)*(xIC-xRC) + (yIC-yRC)*(yIC-yRC) );
      if( distL < distR )
      {
        cl->innerCell = true;
        cr->innerCell = false;
      }
      else
      {
        cl->innerCell = false;
        cr->innerCell = true;
      }
    }
    
    // store division sequences
    if( divType == DivisionType::RADIAL )
    {
      auto iter = _divisionSequences.find( cl->divisionLetterSequence );
      if( iter != _divisionSequences.end() )
        iter->second.at( c->cellFile + 2 )++;
      else
      {
        std::vector<std::size_t> files;
        files.resize( 5, 0 );
        files.at( c->cellFile + 2 )++;
        _divisionSequences.insert( std::make_pair( cl->divisionLetterSequence, files ) );
      }
    }
    
    // check which cell is the upper one and only increase the layer value
    // of the upper one in the case of having a periclinal division if 
    // updateBothLayers is set to false
    // else both daughter cells are assigned a new layer value
    this->setLayerValues( cl, cr, c, divType, true );
    this->checkLayerAppearance( cl );
    this->checkLayerAppearance( cr );
    //this->applyLevelOfDetailToCell( cl );
    //this->applyLevelOfDetailToCell( cr );
    
    this->setCellFileSequence( cl, cr, c, divType, true );
    
    // update precursors
    cl->precursors.insert( c->id );
    cr->precursors.insert( c->id );
    
    for( std::set<std::size_t>::iterator setIter = c->precursors.begin();
         setIter != c->precursors.end(); setIter++ )
    {
      cl->precursors.insert( *setIter );
      cr->precursors.insert( *setIter );
    }
    
    // export division properties
    ModelExporter::exportDivisionDaughterProperties( _divisionFileName,
                                                     cl, cr, c->divType,
                                                     _angleThreshold,
                                                     _divOccurrences,
                                                     false );
  }

  //----------------------------------------------------------------
  
  void applyLevelOfDetailToCell( const cell &c )
  {
    std::vector<junction> juncs;
    forall(const junction& j, T.S.neighbors(c))
      juncs.push_back( j );
    
    if( juncs.size() > _lod * 4 )
      return;
    
    for( auto j = 0; j < juncs.size(); j++ )
    {
      junction newJ;
      
      if( j == juncs.size() - 1 )
      {
        Point3d pos = ( juncs.at(j)->getPos() + juncs.at(0)->getPos() )/2.;
        _VLRBezierSurface.setPos( newJ->sp, pos );
        T.splitWall( juncs.at(j), juncs.at(0), newJ );
      }
      else
      {
        Point3d pos = ( juncs.at(j)->getPos() + juncs.at(j+1)->getPos() )/2.;
        _VLRBezierSurface.setPos( newJ->sp, pos );
        T.splitWall( juncs.at(j), juncs.at(j+1), newJ );
      }
    } 
  }
  
  //----------------------------------------------------------------
  
  void checkLayerAppearance( const cell &c )
  {
    auto iter = _firstLayerAppearances.find( c->layerValue );
    
    // if the layer was not already inserted, it is definitely the
    // first one that will be created by a periclinal division
    if( iter == _firstLayerAppearances.end() )
      _firstLayerAppearances.insert( std::make_pair(
      c->layerValue, std::make_pair( c->timeStep, c->divisionSequence ) ) );
  }
  
  //----------------------------------------------------------------
  
  void setCellProperties( const cell &c, const cell &parentCell )
  {
    double area;
    Point3d center = ModelUtils::computeCellCenter( T, c, area, _accurateCenterOfMass );
    c->initialArea = area;
    c->center = center;
    c->divType = DivisionType::NONE;
    c->centerPos.push_back( center );
    c->area = c->initialArea;
    c->initialLongestWallLength = ModelUtils::determineLongestWallLength( c, T );
    c->longestWallLength = c->initialLongestWallLength;
    c->treeId = parentCell->treeId;
    c->id = _initSurface.getCellID();
    c->parentId = parentCell->id;
    c->timeStep = _initSurface.getTime();
    c->previousAngle = parentCell->angle;
    c->angle = parentCell->angle;
    c->previousDivDir = parentCell->divDir;
    c->divDir = c->center - parentCell->center;
    c->cellCycle = parentCell->cellCycle+1;
    c->cellFile = parentCell->cellFile;
    
    if( parentCell->divType == DivisionType::ANTICLINAL )
      c->divisionLetterSequence = parentCell->divisionLetterSequence + "A";
    else if( parentCell->divType == DivisionType::PERICLINAL )
      c->divisionLetterSequence = parentCell->divisionLetterSequence + "P";
    else
      c->divisionLetterSequence = parentCell->divisionLetterSequence + "R";
    
    if( parentCell->divType == DivisionType::PERICLINAL )
      c->periCycle = parentCell->periCycle+1;
    else
      c->periCycle = parentCell->periCycle;
    _initSurface.incrementCellID();
  }
  
  //----------------------------------------------------------------
  
  void setLayerValues( const cell& cl, const cell& cr, const cell& c,
                       const DivisionType::type divType, const bool updateBothLayers )
  {
    // else check the y values of the center and choose the upper one to be assigned
    // an increased layer value
    if( divType == DivisionType::PERICLINAL )
    {
      if( cl->center.j() > cr->center.j() )
      {
        if( !updateBothLayers )
        {
          cl->layerValue = c->layerValue+1;
          cr->layerValue = c->layerValue;
        }
        else
        {
          cl->layerValue = 2*c->layerValue;
          cr->layerValue = 2*c->layerValue+1;
        }
        
        cl->divisionSequence = c->divisionSequence + "1";
        cr->divisionSequence = c->divisionSequence + "2";
      }
      else
      {
        if( !updateBothLayers )
        {
          cr->layerValue = c->layerValue+1;
          cl->layerValue = c->layerValue;
        }
        else
        {
          cr->layerValue = 2*c->layerValue;
          cl->layerValue = 2*c->layerValue+1;
        }
        
        cr->divisionSequence = c->divisionSequence + "1";
        cl->divisionSequence = c->divisionSequence + "2";
      }
    }
    else
    {
      cl->layerValue = c->layerValue;
      cl->divisionSequence = c->divisionSequence;
      cr->layerValue = c->layerValue;
      cr->divisionSequence = c->divisionSequence;
    }
  }
  
  //----------------------------------------------------------------
  
  void setCellFileSequence( const cell& cl, const cell& cr, const cell& c,
                            const DivisionType::type divType,
                            const bool updateBothFiles )
  {
    if( divType == DivisionType::RADIAL )
    {
      if( cl->center.i() < cr->center.i() )
      {
        if( !updateBothFiles )
        {
          cl->cellFileColoringIndex = c->cellFileColoringIndex+1;
          cr->cellFileColoringIndex = c->cellFileColoringIndex;
        }
        else
        {
          cl->cellFileColoringIndex = 2*c->cellFileColoringIndex+2;
          cr->cellFileColoringIndex = 2*c->cellFileColoringIndex+3;
        }
        
        cl->cellFileSequence = c->cellFileSequence + "0";
        cr->cellFileSequence = c->cellFileSequence + "1";
      }
      else
      {
        if( !updateBothFiles )
        {
          cr->cellFileColoringIndex = c->cellFileColoringIndex+1;
          cl->cellFileColoringIndex = c->cellFileColoringIndex;
        }
        else
        {
          cr->cellFileColoringIndex = 2*c->cellFileColoringIndex+2;
          cl->cellFileColoringIndex = 2*c->cellFileColoringIndex+3;
        }
        
        cl->cellFileSequence = c->cellFileSequence + "1";
        cr->cellFileSequence = c->cellFileSequence + "0";
      }
    }
    else
    {
      cl->cellFileSequence = c->cellFileSequence;
      cl->cellFileColoringIndex = c->cellFileColoringIndex;
      cr->cellFileSequence = c->cellFileSequence;
      cr->cellFileColoringIndex = c->cellFileColoringIndex;
    }
  }
  
  //----------------------------------------------------------------
  
  void step()
  {
    std::size_t curTime, maxTime;
    if( SURFACETYPE == 0 )
    {
      curTime = _initSurface.getTime();
      maxTime = _maxTimeSteps;
    }
    else
    {
      curTime = _VLRDataPointSurface.getCurTimeStep();
      maxTime = _VLRDataPointSurface.getMaxTimeStep() - 1;
    }
   
    if( curTime == 200 )
    {
      // check different division sequences and store them
      this->storeDivisionSequences();
      this->stopModel();
    }
    
    if( curTime == maxTime && _avoidTrianglesThreshold <= 0.00001 )
    {
      ModelExporter::exportPropabilityDistribution( 
      "/tmp/propabilityDistribution.csv",
      _probValues, _lengths, _choices );
    }
    
    if( curTime == maxTime && !_loadLastModel )
    {
      std::string fileName = "/tmp/modelDataProp";
      ModelUtils::appendCellSituationType( fileName, _initialSituationType );
      if( _highOrderPattern != 0 )
        fileName += "_HOP" + std::to_string( _highOrderPattern );
      
      fileName += ".csv";
      this->exportDataProperties( fileName );
      
      unsigned int counter = 0;
      std::vector<unsigned int> cycleCounter;
      cycleCounter.resize( 10, 0 );
      forall(const cell& c, T.C)
      {
        if( c->periCycle < cycleCounter.size() )
          cycleCounter.at( c->periCycle )++;
      }
      
//       for( std::size_t c=0; c < cycleCounter.size(); c++ )
//         std::cout << cycleCounter.at(c) << " cells in cycle " << c << std::endl;
    }
    
    this->renderImage();
    
    if( _useLoop && curTime >= maxTime )
    {      
      // before restart the model, save the model information such
      // as the division occurrences as well as the layering
      this->updateLayerCount();
       
      // export information of periclinal divisions
      std::string filename = "/tmp/ModelDivProperties";
      
      if( SURFACETYPE == 1 )
        filename += _realDataName;
      
      ModelUtils::appendCellSituationType( filename, _initialSituationType );
      
      if( _useAlternativeDT )
      {
        filename += _divisionType;
        if( _divisionType != "Besson-Dumais" )
          filename += std::to_string( (unsigned int)_avoidTrianglesThreshold );
      }
      else
        filename += "ShortestWall";
      
      if( _highOrderPattern != 0 )
        filename += "_HOP" + std::to_string( _highOrderPattern );
      
      filename += ".csv";
      ModelExporter::exportModelProperties( filename,
                                            _loopCounter,
                                            _firstLayerAppearances,
                                            _totalLayerCount,
                                            _divOccurrences,
                                            T.C.size(),
                                            _amountLoops,
                                            _loopCounter == 1 );
                
      this->restartModel();
      if( _loopCounter < _amountLoops )
        _loopCounter++;
      else
        this->stopModel();
      return;
    }
    
    if( _renderMovies && curTime >= maxTime )
    {
      if( _renderMoviesIndex < 11 )
        _renderMoviesIndex++;
      else
        this->stopModel();
      
      if( !_loadLastModel )
      {
        std::string fileName = "/tmp/modelDataProp";
        ModelUtils::appendCellSituationType( fileName, _initialSituationType );
        if( _highOrderPattern != 0 )
          fileName += "_HOP" + std::to_string( _highOrderPattern );
        
        fileName += ".csv";
        this->exportDataProperties( fileName );
      }
      
      this->restartModel();
      return;
    }
     
    _initSurface.incrementTime();
    
    for(int i = 0 ; i < stepPerView ; ++i)
    {
      _divSetting->step_divisions( _initSurface.getTime() );
      
      if( !_lastStep )
      {
        _divSetting->step_tracking( _lineageFileName );
        _divSetting->step_cellWalls( _cellWallsFileName );
      }
      
      _divSetting->step_growth( dt );
    }
    this->setStatus();
    
    if( curTime == maxTime && !_lastStep )
    {
      //ModelExporter::exportTimeAgainstCells( _timeAgainstCellsFileName,
                                              //dt, T.C.size(), false );
      
      _lastStep = true;
      //std::cout << "dt: " << dt << std::endl;
      dt -= steps;
    }
  }

  //----------------------------------------------------------------
  
  void updateLayerCount()
  {
    for( auto iter = _firstLayerAppearances.begin();
         iter != _firstLayerAppearances.end(); ++iter )
    {
      std::string str = iter->second.second;
      // only consider sequence lengths longer than 2
      // because before no high order pattern can be observed
      if( str.length() > 2 )
      {
        // skip itself and avoid redundancies
        auto iter2 = iter;
        ++iter2;
        for( ;iter2 != _firstLayerAppearances.end(); ++iter2 )
        {
          std::string str2 = iter2->second.second;
          
          if( str.length() == str2.length() )
          {
            std::string substr = str;
            substr.erase( str.length()-2, 1 );
            std::string substr2 = str2;
            substr2.erase( str2.length()-2, 1 );
            // if the substrings are identical then compare the
            // time steps
            if( substr.compare( substr2 ) == 0 )
            {
              // if the first string sequence appeared earlier
              if( iter->second.first <= iter2->second.first )
              {
                auto iter3 = _totalLayerCount.find( str );
                if( iter3 != _totalLayerCount.end() )
                  iter3->second++;
                else
                  _totalLayerCount.insert( std::make_pair( str, 1 ) );
                
                auto iter4 = _totalLayerCount.find( str2 );
                if( iter4 == _totalLayerCount.end() )
                  _totalLayerCount.insert( std::make_pair( str2, 0 ) );
              }
              else
              {
                auto iter3 = _totalLayerCount.find( str2 );
                if( iter3 != _totalLayerCount.end() )
                  iter3->second++;
                else
                  _totalLayerCount.insert( std::make_pair( str2, 1 ) );
                
                auto iter4 = _totalLayerCount.find( str );
                if( iter4 == _totalLayerCount.end() )
                  _totalLayerCount.insert( std::make_pair( str, 0 ) );
              }
              
              // break the inner loop when we have made the comparison
              break;
            }
          }
        }
      }
    }
  }
  
  //----------------------------------------------------------------
  
  void storeDivisionSequences()
  {
    for( auto iter = _divisionSequences.begin(); iter != _divisionSequences.end(); iter++ )
    {
      std::cout << iter->first << " ";
      for( std::size_t f = 0; f < iter->second.size(); f++ )
        std::cout << iter->second.at(f) << " ";
      
      std::cout << std::endl;
    }
  }
  
 
  
  //----------------------------------------------------------------
  
  void initDraw( Viewer* viewer )
  {
    viewer->setSceneBoundingBox(Vec(-10.0, -10.0, -1.0), Vec(10.0, 10.0, 1.0));
  }

  //----------------------------------------------------------------
  
  void preDraw()
  {
    util::Palette::Color bg = _palette.getColor(bgColor);
    glClearColor(bg.r(), bg.g(), bg.b(), 1);
    T.preDraw();
  }

  //----------------------------------------------------------------
  
  void postDraw()
  {
    T.postDraw();
  }

  //----------------------------------------------------------------
  
  void draw( Viewer* viewer )
  {
    forall(const cell& c, T.C)
    {
      //T.drawBorders = false;
      if( SURFACETYPE == 0 )
      {
        if( !_bezierGrowthSurface )
          T.cellWallWidth = 0.001;
        else
          T.cellWallWidth = 0.3;
      }
      else
        T.cellWallWidth = 0.3;
          
      //T.cellWallMin = 0.0001;
      //T.strictCellWallMin = true;
      if( _drawSpheres )
        GraphicsClass::drawSphere( c->getPos(), 8., 20., 20.,
                                this->cellColor(c), quadratic );
        
        
      if( _drawCells )
      {
        T.drawCell(c, this->cellColor(c)*0.45, Colorf(this->cellColor(c)) );
        //T.drawCell(c, this->cellColor(c), Colorf(this->cellColor(c)*0.65) );
      }
      
      if( _drawCenter )
        GraphicsClass::drawControlPoint( c->center, _palette.getColor(3) );
      
      if( _drawPCLine &&
          _divisionType == "PerToGrowth" &&
          c->centerPos.size() > 2 )
      {
        for( std::size_t j = 0; j < c->centerPos.size(); j++ )
          GraphicsClass::drawControlPoint( c->centerPos.at(j), _palette.getColor(7) );
        
        for( std::size_t l = 0; l < _pcLines.size(); l++ )
          GraphicsClass::drawLine( _pcLines.at(l).first, _pcLines.at(l).second, _palette.getColor(3) );
      }
    }
    
    if( _drawJunctions )
    {
      forall(const junction& j, T.W)
      {
        GraphicsClass::drawControlPoint( j->getPos(),
                                      _palette.getColor(2) );
      }
    }
    
    if( _drawLayerCurves )
    {
      std::map<unsigned int, std::set<Point3d, lessXPos> > cellPos;
      std::map<unsigned int, std::set<cell, differentCell> > cellMap;
      std::map<unsigned int, util::Palette::Color> colors;
      
      bool oneLinePerLayer = false;
      
      forall(const cell& c, T.C)
      {
        unsigned int layer = c->layerValue-1;
        Point3d cp = c->center;
        cp.k() = layer;
        auto iter = cellPos.find( layer );
        if( iter != cellPos.end() )
        {
          iter->second.insert( cp );
          auto mIter = cellMap.find( layer );
          mIter->second.insert( c );
        }
        else
        {          
          // then initialize set of cell positions for new layer
          std::set<Point3d, lessXPos> cellp;
          std::set<cell, differentCell> ce;
          cellp.insert( cp );
          ce.insert( c );
          cellPos.insert( std::make_pair( layer, cellp ) );
          cellMap.insert( std::make_pair( layer, ce ) );
        }
        
        auto cIter = colors.find( layer );
        if( cIter == colors.end() )
        {
          // first insert color of cell layer
          util::Palette::Color color;
          color = _palette.getColor( _layerColorIndex.at(layer%_layerColorIndex.size()) );
          colors.insert( std::make_pair( layer, color ) ); 
        }
      }
        
      if( oneLinePerLayer )
      { 
        // at last draw the curves
        for( auto iter = cellPos.begin(); iter != cellPos.end(); iter++ )
          GraphicsClass::drawBezierCurve( iter->second, colors.at(iter->first), quadratic );
      }
      else
      {
        for( auto l : cellMap )
        {
          unsigned int layer = l.first;
          std::vector< std::set<cell, differentCell> > curves;
          //std::cout << "Layer " << layer << std::endl;
          ModelUtils::splitNonAdjacentCells( l.second, T, curves );
          for( auto i : curves )
          {
            std::set<Point3d, lessXPos> sortedXPos;
            for( auto s : i )
            {
              Point3d cp = s->center;
              cp.k() = layer;
              sortedXPos.insert( cp );
            }
            
            GraphicsClass::drawBezierCurve( sortedXPos, colors.at(layer), quadratic );
          }
        }
      }
    }
    
    // draw control points of bezier surface
    if( SURFACETYPE == 0 )
    {
      conpoi cps = _VLRBezierSurface.getCurrentControlPoints();
      if( _drawControlPoints )
      {
        for( auto i = 0; i < cps.size(); i++ )
          for( auto j = 0; j < cps.at(i).size(); j++ )
            GraphicsClass::drawControlPoint( cps.at(i).at(j),
                                          _palette.getColor(3) );
      }
      
      if( _drawBezierSurface )
        GraphicsClass::drawBezierSurface( cps, _palette.getColor(0) );
    }
  }
  
  //----------------------------------------------------------------
  
  // color for inner cells
  Colorf cellColor(const cell& c)
  {
    switch( _cellColoringType )
    {
      // coloring based on lineage trees/ founder cells
      case 0: return _palette.getColor(c->treeId);
      // coloring based on layer value
      case 1:
      default:
        return _palette.getColor(
          _layerColorIndex.at((c->layerValue-1)%_layerColorIndex.size()) );
      case 2:
        // boundary
        if( T.border( c ) )
          return _palette.getColor(1);
        // inner cell
        else
          return _palette.getColor(2);
      case 3:
        return _palette.getColor(
          _layerColorIndex.at((c->cellFileColoringIndex-1)%_layerColorIndex.size()) );
    }
  }

  //----------------------------------------------------------------
  
  // color for contour of cells
  Colorf contourColor(const cell& c)
  {
    return _palette.getColor(T.contourColor);
  }

  //----------------------------------------------------------------
  
  // Method needed by the tissue
  Point3d position(const cell& c) const
  { 
    return c->getPos();
  }
  
  //----------------------------------------------------------------
  
  // Method needed by the tissue
  Point3d position(const junction& c) const
  {
    return c->getPos();
  }

  //----------------------------------------------------------------
  
  // Method needed by the tissue
  void setPosition(const cell& c, const Point3d& p)
  {
    if( SURFACETYPE == 0 )
      _VLRBezierSurface.setPos( c->sp, p );
    else
      _VLRDataPointSurface.setPos( c->sp, p );
  }
  
  //----------------------------------------------------------------
  
  // Method needed by the tissue
  void setPosition(const junction& j, const Point3d& p)
  { 
    if( SURFACETYPE == 0 )
      _VLRBezierSurface.setPos( j->sp, p );
    else
      _VLRDataPointSurface.setPos( j->sp, p );
  }

  //----------------------------------------------------------------
  
  // Method needed by the tissue
  void setPositionHint(const junction&, const junction&, const junction&, double) {}

  //----------------------------------------------------------------
  
  // Method needed by the tissue
  Point3d normal(const junction& ) const { return Point3d(0,0,1); }
  
  //----------------------------------------------------------------
  
  // Method needed by the tissue
  Point3d normal(const cell& ) const { return Point3d(0,0,1); }
  
  //----------------------------------------------------------------

  void exportDataProperties( const std::string &filename )
  {
    std::ofstream out( filename.c_str(), std::ofstream::out );
    
    // export model properties
    out << _lod << "\n";
    out << _initialSituationType << "\n";
    out << _LODThreshold << "\n";
    out << _avoidTrianglesThreshold << "\n";
    out << _divisionArea << "\n";
    out << _divisionAreaRatio << "\n";
    out << _equalAreaRatio << "\n";
    out << _useAreaRatio << "\n";
    out << _useCombinedAreaRatio << "\n";
    out << _useWallRatio << "\n";
    out << _divisionWallRatio << "\n";
    out << _useAlternativeDT << "\n";
    out << _divisionType << "\n";
    out << _probabilityOfDecussationDivision << "\n";
    out << _angleThreshold << "\n";
    out << _firstDivisionsAreaRatio << "\n";
    out << _secondDivisionsAreaRatio << "\n";
    out << _timeDelay << "\n";
    out << _cellPinch << "\n";
    out << _cellMaxPinch << "\n";
    out << _interpolateBezierSurfaces << "\n";
    out << _randomChoices.size() << " ";
    
    for( std::size_t c = 0; c < _randomChoices.size(); c++ )
      out << _randomChoices.at(c) << " ";
    
    out << "\n";
    
    out.close();
  }
  
  //----------------------------------------------------------------

  void readDataProperties( const std::string &filename )
  {
    std::ifstream in( filename.c_str(), std::ofstream::in );
    
    // export model properties
    in >> _lod
       >> _initialSituationType
       >> _LODThreshold
       >> _avoidTrianglesThreshold
       >> _divisionArea
       >> _divisionAreaRatio
       >> _equalAreaRatio
       >> _useAreaRatio
       >> _useCombinedAreaRatio
       >> _useWallRatio
       >> _divisionWallRatio
       >> _useAlternativeDT
       >> _divisionType
       >> _probabilityOfDecussationDivision
       >> _angleThreshold
       >> _firstDivisionsAreaRatio
       >> _secondDivisionsAreaRatio
       >> _timeDelay
       >> _cellPinch
       >> _cellMaxPinch
       >> _interpolateBezierSurfaces;
       
    std::size_t size;
    in >> size;
    _randomChoices.resize( size );
    
    for( std::size_t c = 0; c < _randomChoices.size(); c++ )
      in >> _randomChoices.at(c);
    
    in.close();
  }
};

#include "model.moc"
DEFINE_MODEL(MyModel);
