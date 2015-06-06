#include "ModelHeader.h"
#include "ModelExporter.h"
#include "ModelUtils.h"
#include "SurfaceClass.h"
#include "PrincipalComponentAnalysis.h"

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
  util::Palette _palette;
  MyTissue T;
#ifndef TimeAgainstCellAnalysis
  double dt;
#endif
  bool useAreaRatio;
  bool useCombinedAreaRatio;
  bool useWallRatio;
  double divisionArea;
  double divisionAreaRatio;
  double divisionWallRatio;
  int bgColor;
  int stepPerView;
  std::size_t _initialCellNumber;
  std::string _initialCellsOfRealData;
  double probabilityOfDecussationDivision;
  double _angleThreshold;
  bool _bezierGrowthSurface;
  double _surfaceScale;
  bool _useAutomaticContourPoints;
  std::string _lineageFileName;
  std::string _cellWallsFileName;
  std::string _divisionFileName;
  //std::string _timeAgainstCellsFileName;
  std::vector<std::size_t> _layerColorIndex;
  double _cellPinch;
  double _cellMaxPinch;
  std::size_t _cellColoringType;
  std::pair<std::size_t, std::size_t> _divOccurrences;
  SurfaceClass _surfaceClass;
  std::size_t _lod;
  bool _lastStep;
  std::size_t _initialSituationType;
  Point3d _initialBoundary;
  double _firstDivisionsAreaRatio;
  double _secondDivisionsAreaRatio;
  std::size_t _timeDelay;
  std::size_t _timeSixCellStage;
  std::size_t _timeFourCellStage;
  bool _centerOfMassAfterLOD;
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
  bool _interpolateBezierSurfaces;
  bool _drawSpheres;
  GLUquadricObj *quadratic;
  
  std::vector<std::vector<double> > _probValues;
  std::vector<std::vector<double> > _lengths;
  std::vector<std::size_t> _choices;
  
  //----------------------------------------------------------------
  
  void readParms()
  {
    // read the parameters here
#ifndef TimeAgainstCellAnalysis
    _parms("Main", "Dt", dt);
#endif
    _parms("Main", "InitialCellNumber", _initialCellNumber);
    _parms("Main", "InitialCellsOfRealData", _initialCellsOfRealData);
    _parms("Main", "SubDivisionLevelOfCells", _lod);
    _parms("Main", "BezierGrowthSurface", _bezierGrowthSurface);
    _parms("Main", "SurfaceScale", _surfaceScale);
    _parms("Main", "UseAutomaticContourPoints", _useAutomaticContourPoints );
    _parms("Main", "InitialSituationType", _initialSituationType );
    _parms("Main", "CenterOfMassAfterLOD", _centerOfMassAfterLOD );
    _parms("Main", "LODThreshold", _LODThreshold );
    _parms("Main", "AvoidTrianglesThreshold", _avoidTrianglesThreshold );

    _parms("View", "StepPerView", stepPerView);
    _parms("View", "BackgroundColor", bgColor);

    _parms( "Division", "DivisionArea", divisionArea);
    _parms( "Division", "DivisionAreaRatio", divisionAreaRatio);
    _parms( "Division", "UseAreaRatio", useAreaRatio);
    _parms( "Division", "UseCombinedAreaRatio", useCombinedAreaRatio);
    _parms( "Division", "UseWallRatio", useWallRatio);
    _parms( "Division", "DivisionWallRatio", divisionWallRatio);
    _parms( "Division", "UseAlternativeDivisionType", _useAlternativeDT );
    _parms( "Division", "DivisionType", _divisionType );
    _parms( "Division", "ProbabilityOfDecussationDivision", probabilityOfDecussationDivision );
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
    _useLoop( false ),
    _amountLoops( 100 ),
    _drawControlPoints( false ),
    _drawSpheres( false ),
    _interpolateBezierSurfaces( true )
  {
    quadratic = gluNewQuadric();
    gluQuadricNormals( quadratic, GLU_SMOOTH );
    
    readParms();
    // Registering the configuration files
    registerFile("pal.map");
    registerFile("view.v");
    
    std::size_t t = 300;
    if( _initialCellsOfRealData == "130508_raw" )
      t = 350;
    else if( _initialCellsOfRealData.compare( 0, 7, "Average") == 0 )
      t = 150;
    else if( _bezierGrowthSurface && _initialCellsOfRealData == "none" )
      t = 500;
    
    // check if correct surface type was chosen
    if( _bezierGrowthSurface && SURFACETYPE == 1 )
    {
      std::cout << "Wrong surface type chosen!!!" << std::endl;
      return;
    }
    
    _surfaceClass.init( _lod, _lineageFileName, t, _initialSituationType != 0 );
    
    _VLRBezierSurface.init( _parms, "Surface", _bezierGrowthSurface,
                            _interpolateBezierSurfaces, _initialCellsOfRealData );
    
    // set name strings
    _lineageFileName = "/tmp/model";
    _cellWallsFileName = "/tmp/modelCellWalls";
    _divisionFileName = "/tmp/divisionPropertiesModel";
    
    switch( _initialSituationType )
    {
      case 0:
        _divisionFileName += "NH";
        _lineageFileName += "NH";
        _cellWallsFileName += "NH";
        break;
      case 1:
        _divisionFileName += "1DC";
        _lineageFileName += "1DC";
        _cellWallsFileName += "1DC";
        break;
      case 2:
        _divisionFileName += "2DC";
        _lineageFileName += "2DC";
        _cellWallsFileName += "2DC";
      break;
    }
    if( _useAlternativeDT )
    {
      _divisionFileName += _divisionType;
      if( _divisionType == "Energy" ||
          _divisionType == "Random" ||
          _divisionType == "Gibbs" )
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
    
    _divisionFileName += ".csv";
    _lineageFileName += ".csv";
    _cellWallsFileName += ".csv";
    
    // single layer assignment
    //for( std::size_t l = 9; l < 14; l++ )
    // multiple layer assignment for each new daughter cell
    for( std::size_t l = 14; l < 45; l++ )
      _layerColorIndex.push_back( l );
    
    cell dummy;
    MyTissue::division_data ddummy;
    
    ModelExporter::exportLineageInformation( _lineageFileName, dummy, T, true );
    ModelExporter::exportCellWalls( _cellWallsFileName, dummy, T, true );
    
    std::pair<std::size_t, std::size_t> pair;
    ModelExporter::exportDivisionDaughterProperties( _divisionFileName,
                                                     dummy, dummy,
                                                     ddummy, 0.,
                                                     pair, true );
    
    /*
    _timeAgainstCellsFileName = "/tmp/timeAgainstCells";
    _timeAgainstCellsFileName += _initialCellsOfRealData;
    _timeAgainstCellsFileName += ".csv";
   
    if( dt > start-steps )
      ModelExporter::exportTimeAgainstCells( _timeAgainstCellsFileName, dt, 0, true );
    */
        
    // bezier
    if( SURFACETYPE == 0 )
    {
      _VLRBezierSurface.GrowStep(0);
      
      // special cases of number of cells at the beginning
      if( _initialCellsOfRealData == "none" )
        _surfaceClass.initModelBasedOnBezier( T, _initialCellNumber,
                                              _VLRBezierSurface );
      // else constellation of founder cells according to real data
      else
        _surfaceClass.initLateralRootBasedOnBezier( T, _initialCellsOfRealData,
                                                    _VLRBezierSurface );
    }
    // real data points
    else
    {
      _VLRDataPointSurface.init( _surfaceScale,
                                 _initialCellsOfRealData,
                                 _useAutomaticContourPoints );
      
      std::vector<TrianglePoint> tps;
      _VLRDataPointSurface.growStep( 0, tps );
        
      _surfaceClass.initLateralRootBasedOnRealData( T, _VLRDataPointSurface,
                                                    _initialCellsOfRealData,
                                                    _surfaceScale,
                                                    _useAutomaticContourPoints );
    }
    
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
  }
  
  //----------------------------------------------------------------
  
  void restartModel()
  {
    // reset the tissue to delete all previous cell information
    T = MyTissue( _palette, this );
    
    readParms();
    // Registering the configuration files
    registerFile("pal.map");
    registerFile("view.v");
    
    // reset division and layering results
    _divOccurrences = std::make_pair( 0, 0 );
    _firstLayerAppearances.clear();
    _totalLayerCount.clear();
    
    std::size_t t = 300;
    if( _initialCellsOfRealData == "130508_raw" )
      t = 350;
    else if( _initialCellsOfRealData.compare( 0, 7, "Average") == 0 )
      t = 150;
    else if( _bezierGrowthSurface && _initialCellsOfRealData == "none" )
      t = 502;
    
    _surfaceClass.init( _lod, _lineageFileName, t, _initialSituationType != 0 );
    
    _VLRBezierSurface.init( _parms, "Surface", _bezierGrowthSurface,
                            _interpolateBezierSurfaces, _initialCellsOfRealData );
    
    // bezier
    if( SURFACETYPE == 0 )
    {
      _VLRBezierSurface.GrowStep(0);
      
      // special cases of number of cells at the beginning
      if( _initialCellsOfRealData == "none" )
        _surfaceClass.initModelBasedOnBezier( T, _initialCellNumber,
                                              _VLRBezierSurface );
      // else constellation of founder cells according to real data
      else
        _surfaceClass.initLateralRootBasedOnBezier( T, _initialCellsOfRealData,
                                                    _VLRBezierSurface );
    }
    // real data points
    else
    {
      _VLRDataPointSurface.init( _surfaceScale,
                                 _initialCellsOfRealData,
                                 _useAutomaticContourPoints );
      
      std::vector<TrianglePoint> tps;
      _VLRDataPointSurface.growStep( 0, tps );
        
      _surfaceClass.initLateralRootBasedOnRealData( T, _VLRDataPointSurface,
                                                    _initialCellsOfRealData,
                                                    _surfaceScale,
                                                    _useAutomaticContourPoints );
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
  }
  
  //----------------------------------------------------------------
  
  void setStatus()
  {
    //double time = _VLRBezierSurface.GetTime();
    std::size_t time = _surfaceClass.getTime();
//     if( time > _surfaceClass.getMaxTime() )
//       time = _surfaceClass.getMaxTime();
    
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
    status += QString( "AD: %1\t" ).arg(_divOccurrences.first);
    status += QString( "PD: %1\t" ).arg(_divOccurrences.second);
    
    if( useCombinedAreaRatio )
    {
      status += QString( "Area div: %1\t" ).arg(divisionArea);
      status += QString( "Area div ratio: %1\t" ).arg(divisionAreaRatio);
    }
    else if( useAreaRatio )
      status += QString( "Area div ratio: %1\t" ).arg(divisionAreaRatio);
    else
      status += QString( "Area div: %1\t" ).arg(divisionArea);
    
    if( useWallRatio )
      status += QString( "Wall div ratio: %1\t" ).arg(divisionWallRatio);
    
    if( _divisionType == "Decussation" )
      status += QString( "Dec prop: %1\%\t" ).arg(probabilityOfDecussationDivision);
    
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
    
    // check which cell is the upper one and only increase the layer value
    // of the upper one in the case of having a periclinal division if 
    // updateBothLayers is set to false
    // else both daughter cells are assigned a new layer value
    this->setLayerValues( cl, cr, c, divType, true );
    this->checkLayerAppearance( cl );
    this->checkLayerAppearance( cr );
    
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
                                                     cl, cr, ddata,
                                                     _angleThreshold,
                                                     _divOccurrences,
                                                     false );
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
    Point3d center;
    double area;
    this->setCellCenter( c, center, area );
    c->initialArea = area;
    c->center = center;
    c->divType = DivisionType::NONE;
    c->centerPos.push_back( center );
    c->area = c->initialArea;
    c->initialLongestWallLength = ModelUtils::determineLongestWallLength( c, T );
    c->longestWallLength = c->initialLongestWallLength;
    c->treeId = parentCell->treeId;
    c->id = _surfaceClass.getCellID();
    c->parentId = parentCell->id;
    c->timeStep = _surfaceClass.getTime();
    c->previousAngle = parentCell->angle;
    c->angle = parentCell->angle;
    c->previousDivDir = parentCell->divDir;
    c->divDir = c->center - parentCell->center;
    c->cellCycle = parentCell->cellCycle+1;
    _surfaceClass.incrementCellID();
  }

  //----------------------------------------------------------------
  
  void setCellCenter( const cell &c, Point3d &center, double &area )
  {
    center = Point3d( 0., 0., 0. );
    if( !_centerOfMassAfterLOD )
    {
      std::vector<Point2d> polygon;
      forall(const junction& j, T.S.neighbors(c))
      {
        Point3d pos = j->getPos();
        polygon.push_back( Point2d( pos.i(), pos.j() ) );
        center += pos;
      }
      
      center /= polygon.size();
      area = geometry::polygonArea(polygon);
    }
    else
    {
      center = ModelUtils::getCenterAfterApplyingLODToCell( c, T, area,
                                                            _LODThreshold );
    }
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
  
  void determineXMinMax( const cell &c )
  {
    double xMin = 5000000.;
    double xMax = -5000000.;
    forall(const junction& j, T.S.neighbors(c))
    {
      Point3d pos = j->getPos();
      
      if( pos.i() < xMin )
        xMin = pos.i();
      
      if( pos.i() >= xMax )
        xMax = pos.i();
    }
    
    // set min and max values for x
    c->xMin = xMin;
    c->xMax = xMax;
  }
  
  //----------------------------------------------------------------
  
  bool setNextDecussationDivision()
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    double p = probabilityOfDecussationDivision;
    std::discrete_distribution<int> d({p, 100.-p});
    
    if(d(gen) == 0)
      return true;
    else
      return false;
  }

  //----------------------------------------------------------------
  
  MyTissue::division_data getRandomDivisionData( const cell& c,
                                                 bool &empty )
  {
    std::vector<MyTissue::division_data> divData;
    divData = ModelUtils::determinePossibleDivisionData(
      c, _avoidTrianglesThreshold, _LODThreshold, T );
    
    //std::cout << "possible Divisions: " << divData.size() << std::endl;
    
    MyTissue::division_data ddata;
    if( divData.size() != 0 )
    {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<int> d( 0, divData.size()-1 );
      
      // get a random choice of all possible division data
      std::size_t choice = d(gen);
      //std::cout << "c: " << choice << std::endl;
      ddata = divData.at( choice );
      vvcomplex::testDivisionOnVertices(c, ddata, T, 0.01);
      
      // apply cell pinching
      tissue::CellPinchingParams params;
      params.cellPinch = _cellPinch;
      params.cellMaxPinch = _cellMaxPinch;
      tissue::cellPinching( c, T, ddata, params );
      empty = false;
    }
    else
      empty = true;
    
    return ddata;
  }
  
  //----------------------------------------------------------------
  
  MyTissue::division_data getEnergyDivisionData( const cell& c,
                                                 bool &empty )
  {
    std::vector<MyTissue::division_data> divData;
    divData = ModelUtils::determinePossibleDivisionData(
      c, _avoidTrianglesThreshold, _LODThreshold, T );
    
    // compute the lengths of all division lines and sort them
    // in ascending order
    std::vector<double> lengths;
    lengths.resize( divData.size() );
    double minLength = 5000.;
    double maxLength = 0.;
    double sumLength = 0.;
    std::size_t minIndex = 0;
    
    for( std::size_t l=0; l<lengths.size(); l++ )
    {
      Point3d pu = divData.at(l).pu;
      Point3d pv = divData.at(l).pv;
      //std::cout << "pu: " << pu << std::endl;
      //std::cout << "pv: " << pv << std::endl;
      double length = norm( pu - pv );
      lengths.at(l) = length;
      
      if( length < minLength )
      {
        minLength = length;
        minIndex = l;
      }
      
      if( length >= maxLength )
        maxLength = length;
      
      sumLength += length;
    }
    
    //std::cout << "possible Divisions: " << divData.size() << std::endl;
    //std::cout << "minLength: " << minLength << std::endl;
    //std::cout << "maxLength: " << maxLength << std::endl;
    //std::cout << "sumLength: " << sumLength << std::endl;
    
    std::size_t probDist = 1;
    
    MyTissue::division_data ddata;
    if( divData.size() != 0 )
    {
      std::vector<double> probValues;
      probValues.resize( divData.size() );
      
      if( probDist == 2 )
      {
        // get the probability distribution function for all lengths
        for( std::size_t l=0; l<probValues.size(); l++ )
        {
          double prob = (lengths.at(l)*100.)/sumLength;
          // the current probability value assigns a higher prob to longer lengths;
          // however, this should be inverse. Consequently, we substract the prob
          // from 100 percent and divide it by lengths.size()-1 to get the inverse
          // probability value
          if( lengths.size() > 1 )
            probValues.at(l) = (100. - prob);// /(lengths.size()-1);
          else
            probValues.at(l) = 100.;
          
          //std::cout << "prob value: " << probValues.at(l) << std::endl;
        }
      }
      else
      {
        double mu = 0.;
        double sd = std::sqrt(0.5);
        std::vector<double> normLengths;
        normLengths.resize( divData.size() );
        for( std::size_t l=0; l<normLengths.size(); l++ )
        {
          normLengths.at(l) = (lengths.at(l)-minLength)/(maxLength-minLength);
          double prob = ModelUtils::getNormalDistribution( normLengths.at(l), mu, sd );
      
          if( lengths.size() > 1 )
            probValues.at(l) = prob*100.;
          else
            probValues.at(l) = 100.;
          
          //std::cout << l << " prob value: " << probValues.at(l) << std::endl;
        }
      }
      
      std::size_t choice = ModelUtils::getRandomResultOfDistribution( probValues );
      //std::cout << "choice: " << choice << std::endl;
      ddata = divData.at( choice );
      vvcomplex::testDivisionOnVertices(c, ddata, T, 0.01);
      
      // apply cell pinching
      tissue::CellPinchingParams params;
      params.cellPinch = _cellPinch;
      params.cellMaxPinch = _cellMaxPinch;
      tissue::cellPinching( c, T, ddata, params );
      empty = false;
    }
    else
      empty = true;
    
    return ddata;
  }
  
  //----------------------------------------------------------------
  
  MyTissue::division_data getGibbsDivisionData( const cell& c,
                                                bool &empty )
  {
    std::vector<MyTissue::division_data> divData;
    divData = ModelUtils::determinePossibleDivisionData(
      c, _avoidTrianglesThreshold, _LODThreshold, T );
    
    std::size_t numLines = divData.size();
//    std::cout << "possible Divisions: " << numLines << std::endl;
    
    // fixed parameter from paper determined by experimental results
    const double beta = 20.6;
    
    // compute the lengths of all division lines and sort them
    // in ascending order
    std::vector<double> lengths;
    lengths.resize( numLines );
    
    // compute area size and mean cell diameter
    std::vector<Point2d> polygon;
    forall(const junction& j, T.S.neighbors(c))
    {
      Point3d pos = j->getPos();
      polygon.push_back( Point2d( pos.i(), pos.j() ) );
    }
    double area = geometry::polygonArea( polygon );
    // mean cell diameter
    double rho = std::sqrt( area );
    
    for( std::size_t l=0; l<numLines; l++ )
    {
      Point3d pu = divData.at(l).pu;
      Point3d pv = divData.at(l).pv;
      //std::cout << "pu: " << pu << std::endl;
      //std::cout << "pv: " << pv << std::endl;
      lengths.at(l) = norm( pu - pv );
      //std::cout << l << " length: " << lengths.at(l) << std::endl;
    }
    
    MyTissue::division_data ddata;
    if( numLines != 0 )
    {
      double sum = 0.;
      for( std::size_t l=0; l<numLines; l++ )
        sum += std::exp( (-1. * beta * lengths.at(l)) / rho );
      
      std::vector<double> probValues;
      probValues.resize( numLines );
      for( std::size_t l=0; l<numLines; l++ )
      {
        probValues.at(l) = 100. * std::exp( (-1. * beta * lengths.at(l)) / rho ) / sum;
        //std::cout << l << " prob value: " << probValues.at(l) << std::endl;
      }
      
      std::size_t choice = ModelUtils::getRandomResultOfDistribution( probValues );
      //std::cout << "choice: " << choice << std::endl;
      ddata = divData.at( choice );
      vvcomplex::testDivisionOnVertices(c, ddata, T, 0.01);
      
      _probValues.push_back( probValues );
      _lengths.push_back( lengths );
      _choices.push_back( choice );
      
      // apply cell pinching
      tissue::CellPinchingParams params;
      params.cellPinch = _cellPinch;
      params.cellMaxPinch = _cellMaxPinch;
      tissue::cellPinching( c, T, ddata, params );
      empty = false;
    }
    else
      empty = true;
    
    return ddata;
  }
  
  //----------------------------------------------------------------
  
  MyTissue::division_data setDivisionPoints( const cell& c )
  {
    MyTissue::division_data ddata;
    const Point3d& center = c->center;
    double a = M_PI/180. * c->angle;
    Point3d direction = Point3d(-sin(a), cos(a), 0);
    forall( const junction& j,T.S.neighbors(c) )
    {
      Point3d jpos, jnpos;
      const junction& jn = T.S.nextTo(c, j);
      jpos = j->getPos();
      jnpos = jn->getPos();
      
      Point3d u;
      double s;
      if(geometry::planeLineIntersection(u, s, center, direction, jpos, jnpos) and s >= 0 and s <= 1)
      {
        if((jpos - center)*direction > 0)
        {
          ddata.v1 = j;
          ddata.pv = u;
        }
        else if((jpos - center)*direction < 0)
        {
          ddata.u1 = j;
          ddata.pu = u;
        }
        if(ddata.u1 and ddata.v1)
          break;
      }
    }
    
    vvcomplex::testDivisionOnVertices(c, ddata, T, 0.01);
    
    // apply cell pinching
    tissue::CellPinchingParams params;
    params.cellPinch = _cellPinch;
    params.cellMaxPinch = _cellMaxPinch;
    tissue::cellPinching( c, T, ddata, params );
    
    return ddata;
  }
  
  //----------------------------------------------------------------
  
  bool step_divisions()
  {
    // for which number of cells should the division area ratio check apply
    // this is required since at the beginning only few divisions occur due
    // to the small increasing of area based on the initial area size
    std::size_t areaRatioStart;
    if( T.C.size() > 1 )
      areaRatioStart = 0;
    else
      areaRatioStart = 1;
    
    // Find cells to be divided
    std::list<cell> to_divide;
    
    // wait after the six cell stage until the predefined time has passed
    bool wait = false;
    if( T.C.size() == 4 && _initialSituationType == 1 )
    {
      std::size_t curTime = _surfaceClass.getTime();
      if( curTime - _timeFourCellStage > _timeDelay )
        wait = false;
      else
        wait = true;
    }
    
    if( T.C.size() == 6 && _initialSituationType == 2 )
    {
      std::size_t curTime = _surfaceClass.getTime();
      if( curTime - _timeSixCellStage > _timeDelay )
        wait = false;
      else
        wait = true;
    }
    
    // at first determine the min and max values for each cell
    forall(const cell& c, T.C)
      this->determineXMinMax(c);
    
    // wait for the four-cell stage some time steps such that the
    // the future divisions are not occurring too fast
    // and only update the center of cells
    if( T.C.size() == 4 &&
        _initialSituationType == 1 && wait )
    {
      forall(const cell& c, T.C)
      {
        // update center position
        Point3d center;
        double area;
        this->setCellCenter( c, center, area );
        this->setPosition( c, center );
        c->area = area;
        c->center = center;
        // add the current center position to the set
        c->centerPos.push_back( center );
        c->timeStep = _surfaceClass.getTime();
        c->longestWallLength = ModelUtils::determineLongestWallLength( c, T );
      }
    }
    // force the initial start of the VLR
    else if( T.C.size() < 6 && _initialSituationType != 0 )
    {
      forall(const cell& c, T.C)
      {
        // update center position
        Point3d center;
        double area;
        this->setCellCenter( c, center, area );
        if( T.C.size() < 4 )
        {
          double xLength = std::fabs( c->xMax - c->xMin );
          if( c->id == 1 )
            center.i() = c->xMin + 2.*xLength/3.; 
          else if( c->id == 2 )
            center.i() = c->xMin + 1.*xLength/3.;
        }
        this->setPosition( c, center );
        c->area = area;
        c->center = center;
        // add the current center position to the set
        c->centerPos.push_back( center );
        c->timeStep = _surfaceClass.getTime();
        c->longestWallLength = ModelUtils::determineLongestWallLength( c, T );
        
        double a = c->area;
        double l = c->longestWallLength;
        
        // divide cells if their area size has increased by a certain percentage amount
        double initialArea = c->initialArea;
        if( T.C.size() < 4 )
          initialArea += initialArea*_firstDivisionsAreaRatio;
        else
          initialArea += initialArea*_secondDivisionsAreaRatio;
        
        if( T.C.size() == 3 )
        {
          // when one of the two founder cells has already divided
          // then assure that the second founder cell will divide
          // before one of the two daughter cells of the first division
          // divide again
          if( c->id == 2 || c->id == 1 )
          {
            if( a > initialArea )
              to_divide.push_back(c);
          }
          
          _timeFourCellStage = _surfaceClass.getTime();
        }
        // else perform the "normal" division routine
        else
        {
          if( a > initialArea )
          {
            // if there are 5 cells and if we are in this loop
            // a next division will occur; so we store the time
            // of the six cell stage
            if( T.C.size() == 5 )
              _timeSixCellStage = _surfaceClass.getTime();
              
            to_divide.push_back(c);
          }
        }
      }
    }
    // wait for the six-cell stage some time steps such that the
    // the future divisions are not occurring too fast
    // and only update the center of cells
    else if( T.C.size() == 6 &&
             _initialSituationType != 0 && wait )
    {
      forall(const cell& c, T.C)
      {
        // update center position
        Point3d center;
        double area;
        this->setCellCenter( c, center, area );
        this->setPosition( c, center );
        c->area = area;
        c->center = center;
        // add the current center position to the set
        c->centerPos.push_back( center );
        c->timeStep = _surfaceClass.getTime();
        c->longestWallLength = ModelUtils::determineLongestWallLength( c, T );
      }
    }
    else
    {      
      forall(const cell& c, T.C)
      {
        // update center position
        Point3d center;
        double area;
        this->setCellCenter( c, center, area );
        this->setPosition( c, center );  
        c->area = area;
        c->center = center;
        // add the current center position to the set
        c->centerPos.push_back( center );
        c->timeStep = _surfaceClass.getTime();
        c->longestWallLength = ModelUtils::determineLongestWallLength( c, T );
        
        double a = c->area;
        double l = c->longestWallLength;
        if( useAreaRatio && useWallRatio && c->id > areaRatioStart )
        {
          // divide cells if their area size has increased by a certain percentage amount
          double initialArea = c->initialArea;
          initialArea += initialArea*divisionAreaRatio;
          // divide cell if its wall length has increased by a certain percentage amount
          double initialLongestLength = c->initialLongestWallLength;
          initialLongestLength += initialLongestLength*divisionWallRatio;
          
          if( a > initialArea || l > initialLongestLength )
            to_divide.push_back(c);
        } 
        // apply a division if the cells area exceeds a certain threshold or area ratio
        else if( useCombinedAreaRatio && c->id > areaRatioStart )
        {
          // divide cells if their area size has increased by a certain percentage amount
          double initialArea = c->initialArea;
          initialArea += initialArea*divisionAreaRatio;
          if( a > initialArea && a > divisionArea )
            to_divide.push_back(c);
        }
        // only apply the division based on ratio with at least areaRatioStart cells
        else if( useAreaRatio && c->id > areaRatioStart )
        {
          // divide cells if their area size has increased by a certain percentage amount
          double initialArea = c->initialArea;
          initialArea += initialArea*divisionAreaRatio;
          if( a > initialArea )
            to_divide.push_back(c);
        }
        // only apply the division based on ratio with at least areaRatioStart cells
        else if( useWallRatio && c->id > areaRatioStart )
        {
          // divide cell if its wall length has increased by a certain percentage amount
          double initialLongestLength = c->longestWallLength;
          initialLongestLength += initialLongestLength*divisionWallRatio;
          if( l > initialLongestLength )
            to_divide.push_back(c);
        }
        else
        {
          if( a > divisionArea )
            to_divide.push_back(c);
        }
      }
    }
    
    // Divide the cells
    forall(const cell& c, to_divide)
    {
      // perform the "normal" division routine for the initial two founder cells
      if( T.C.size() < 4 && _initialSituationType != 0 )
        T.divideCell(c);
      // set the division properties for the second division
      // periclinal division resulting in six cells
      // **********************************
      // *         *     *     *          *
      // *         *************          *
      // *         *     *     *          *
      // **********************************
      else if( T.C.size() > 3 && T.C.size() < 6 && _initialSituationType == 2 )
      {
        if( c->innerCell )
        {
          c->angle = 180.;
          MyTissue::division_data ddata = this->setDivisionPoints( c );
          T.divideCell( c, ddata );
        }
      }
      else if( _divisionType == "Decussation" &&
               c->id > areaRatioStart &&
               _useAlternativeDT )
      {
        // if true then the next division is perpendicular to the last one
        // else the division is collinear to the previous division direction
        if( this->setNextDecussationDivision() )
          c->angle = fmod( c->angle + 90., 360. );
        
        MyTissue::division_data ddata = this->setDivisionPoints( c );
        T.divideCell( c, ddata );
      }
      else if( _divisionType == "PerToGrowth" &&
               c->id > areaRatioStart &&
               _useAlternativeDT )
      {
        // if true then the next division is perpendicular to the principal
        // component growth of the center deformation since the cell was born
        // determine principal growth direction
        if( c->centerPos.size() > 1 )
          c->principalGrowthDir = this->determineLongestPCGrowth( c->centerPos );
        else
          std::cout << "Too few center positions!" << std::endl;
        
        Point3d xaxisDir = Point3d( 1., 0., 0. );
        c->angle = 180./M_PI * acos( c->principalGrowthDir*xaxisDir ) + 90.;
        // update the division angle
        MyTissue::division_data ddata = this->setDivisionPoints( c );
        T.divideCell( c, ddata );
      }
      else if( _divisionType == "Energy" &&
               c->id > areaRatioStart &&
               _useAlternativeDT )
      {
        // if true then all division planes are determined preserving
        // almost non-triangle cells with same area sizes and the shortest cell wall
        // has the lowest energy while the longest one has the highest energy;
        // probability values are then assigned like this:
        // lowest energy -> high probability
        // highest energy -> low probability
        bool empty = false;
        MyTissue::division_data ddata = this->getEnergyDivisionData( c, empty );
        // if the division data is empty then just use the default division
        // rule (e.g. ShortestWall)
        if( empty )
          T.divideCell( c );
        else
          T.divideCell( c, ddata );
      }
      else if( _divisionType == "Gibbs" &&
               c->id > areaRatioStart &&
               _useAlternativeDT )
      {
        // the division plane is chosen based on Gibbs measure as described
        // in the paper of Besson and Dumais, 2011 (Universal rule for the 
        // symmetric division of plant cells)
        bool empty = false;
        MyTissue::division_data ddata = this->getGibbsDivisionData( c, empty );
        // if the division data is empty then just use the default division
        // rule (e.g. ShortestWall)
        if( empty )
          T.divideCell( c );
        else
          T.divideCell( c, ddata );
      }
      else if( _divisionType == "Random" &&
               c->id > areaRatioStart &&
               _useAlternativeDT )
      {
        // if true then all division planes are determined randomly preserving
        // almost symmetric cells with same area sizes
        bool empty = false;
        MyTissue::division_data ddata = this->getRandomDivisionData( c, empty );
        // if the division data is empty then just use the default division
        // rule (e.g. ShortestWall)
        if( empty )
          T.divideCell( c );
        else
          T.divideCell( c, ddata );
      }
      else
        T.divideCell(c);
    }

    return !to_divide.empty();
  }
  
  //----------------------------------------------------------------
  
  void step_cellWalls()
  {    
    if( !_lastStep )
    {
      forall(const cell& c, T.C)
        ModelExporter::exportCellWalls( _cellWallsFileName, c, T, false );
    }
  }
  
  //----------------------------------------------------------------
  
  void step_tracking()
  {
    if( !_lastStep )
    {
      forall(const cell& c, T.C)
        ModelExporter::exportLineageInformation( _lineageFileName, c, T, false );
    }
  }
  
  //----------------------------------------------------------------
  
  void step_growth()
  {
    if( SURFACETYPE == 0 )
    {
      _VLRBezierSurface.GrowStep(dt);
      forall(const junction& v, T.W)
        _VLRBezierSurface.GetPos(v->sp);
    }
    else
    {
      // generate new triangulation surface based on real data points
      std::vector<TrianglePoint> tps;
      forall(const junction& j, T.W)
        tps.push_back(j->tp);

      _VLRDataPointSurface.growStep( dt, tps );
      
      int i = 0;
      forall(const junction& j, T.W)
      {
        j->tp = tps[i++];
        _VLRDataPointSurface.getPos( j->tp );
      }
    }
  }
  
  //----------------------------------------------------------------
  
  void step()
  {  
    if( _useLoop &&
        _surfaceClass.getTime() > _surfaceClass.getMaxTime() )
    {
      /*
      // before restart the model, save the model information such
      // as the division occurrences as well as the layering
      for( auto iter = _firstLayerAppearances.begin();
           iter != _firstLayerAppearances.end(); ++iter )
        std::cout << "Layer: " << iter->first << " t: "
                  << iter->second.first << " seq: "
                  << iter->second.second << std::endl;
      
      std::cout << "DivOcc: " << _divOccurrences.first << ", "
                << _divOccurrences.second << std::endl;
      */
      
      this->updateLayerCount();
        
      /*
      for( auto iter = _totalLayerCount.begin();
           iter != _totalLayerCount.end(); ++iter )
      {
        std::cout << "Seq: " << iter->first << " count: "
                  << iter->second << std::endl;
      }
      */
       
      // export information of periclinal divisions
      std::string filename = "/tmp/ModelDivProperties";
      switch( _initialSituationType )
      {
        case 0: filename += "NH"; break;
        case 1: filename += "1DC"; break;
        case 2: filename += "2DC"; break;
      }
      
      if( _useAlternativeDT )
      {
        filename += _divisionType;
        if( _divisionType == "Energy" ||
            _divisionType == "Random" ||
            _divisionType == "Gibbs" )
          filename += std::to_string( (unsigned int)_avoidTrianglesThreshold );
      }
      else
        filename += "ShortestWall";
      
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
    
    if( _surfaceClass.getTime() == _surfaceClass.getMaxTime() )
    {
      ModelExporter::exportPropabilityDistribution( 
      "/tmp/propabilityDistribution.csv",
      _probValues, _lengths, _choices );
    }
    
    _surfaceClass.incrementTime();
    
    for(int i = 0 ; i < stepPerView ; ++i)
    {
      this->step_divisions();
      this->step_tracking();
      this->step_cellWalls();
      this->step_growth();
    }
    this->setStatus();
              
    std::size_t curTime, maxTime;
    
    if( SURFACETYPE == 0 )
    {
      curTime = _surfaceClass.getTime();
      maxTime = _surfaceClass.getMaxTime();
    }
    else
    {
      curTime = _VLRDataPointSurface.getCurTimeStep();
      maxTime = _VLRDataPointSurface.getMaxTimeStep() - 1;
    }
    
    if( curTime == maxTime && !_lastStep )
    {
      //ModelExporter::exportTimeAgainstCells( _timeAgainstCellsFileName,
                                              //dt, T.C.size(), false );
      
      _lastStep = true;
      std::cout << "dt: " << dt << std::endl;
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
  
  Point3d determineLongestPCGrowth( const std::vector<Point3d> &positions )
  {
    Point3d principalGrowthDir;
    
    Eigen::MatrixXd points( positions.size(), 3 );
    
    for( std::size_t r=0;r<positions.size();r++ )
    {
      points( r, 0 ) = positions.at(r).i();
      points( r, 1 ) = positions.at(r).j();
      points( r, 2 ) = positions.at(r).k();
    }
      
    //std::cout << "Matrix: " << points << std::endl;
    std::vector<double> lpc = PCA::compute( points );
    principalGrowthDir = Point3d( lpc.at(0), lpc.at(1), lpc.at(2) );
    std::cout << "LPC: " << principalGrowthDir << std::endl;
    
    return principalGrowthDir;
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
          T.cellWallWidth = 0.2;
      }
      else
        T.cellWallWidth = 0.2;
          
      //T.cellWallMin = 0.0001;
      //T.strictCellWallMin = true;
      if( _drawSpheres )
        ModelUtils::drawSphere( c->getPos(), 10., 20., 20.,
                                this->cellColor(c), quadratic );
      else
        T.drawCell(c, this->cellColor(c)*0.5, Colorf(this->cellColor(c)) );
    }
    
    // draw control points of bezier surface
    if( _drawControlPoints && SURFACETYPE == 0 )
    {
      conpoi cps = _VLRBezierSurface.getCurrentControlPoints();
      for( auto i = 0; i < cps.size(); i++ )
        for( auto j = 0; j < cps.at(i).size(); j++ )
          ModelUtils::drawControlPoint( cps.at(i).at(j),
                                        _palette.getColor(3) );
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
      if( c->layerValue-1 < _layerColorIndex.size() )
        return _palette.getColor( _layerColorIndex.at(c->layerValue-1) );
      else
        return _palette.getColor( _layerColorIndex.at(_layerColorIndex.size()-1) );
      case 2:
        // boundary
        if( T.border( c ) )
          return _palette.getColor(1);
        // inner cell
        else
          return _palette.getColor(2);
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
      _VLRBezierSurface.SetPoint(c->sp, c->sp, p);
    else
      _VLRDataPointSurface.setPos( c->tp, p );
  }
  
  //----------------------------------------------------------------
  
  // Method needed by the tissue
  void setPosition(const junction& j, const Point3d& p)
  { 
    if( SURFACETYPE == 0 )
      _VLRBezierSurface.SetPoint(j->sp, j->sp, p);
    else
      _VLRDataPointSurface.setPos( j->tp, p );
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
};

#include "model.moc"
DEFINE_MODEL(MyModel);
