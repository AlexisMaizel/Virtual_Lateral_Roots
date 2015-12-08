#include "model.h"

#include "GraphicsClass.h"
#include "ModelExporter.h"

//#define TimeAgainstCellAnalysis

const double start = 0.05;
const double steps = 0.0005;
#ifdef TimeAgainstCellAnalysis
static double dt = start;
#endif

//----------------------------------------------------------------

void MyModel::readParms()
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
void MyModel::modifiedFiles( const std::set<std::string>& filenames )
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

MyModel::MyModel(QObject *parent) : Model(parent), _parms("view.v"),
  _VLRDataPointSurface( _parms, "Surface" ),
  _palette("pal.map"), T(_palette, this),
  _divOccurrences( std::make_pair( 0, 0 ) ),
  _lastStep( false ),
  _loopCounter( 1 ),
  _amountLoops( 2 ),
  _interpolateBezierSurfaces( true ),
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
  
  for( std::size_t l = 45; l < 48; l++ )
    _cellFileColorIndex.push_back( l );
  
  cell dummy;
  ModelExporter::exportLineageInformation( _lineageFileName, dummy, T, true );
  ModelExporter::exportCellWalls( _cellWallsFileName, dummy, T, true );
  
  std::pair<std::size_t, std::size_t> pair;
  std::vector<cell> dummies;
  dummies.push_back( dummy );
  dummies.push_back( dummy );
  ModelExporter::exportDivisionDaughterProperties( _divisionFileName,
                                                    dummies,
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
    _useAreaRatio, _equalAreaRatio,
    _useCombinedAreaRatio, _useWallRatio,
    _divisionArea, _divisionAreaRatio, _divisionWallRatio,
    _LODThreshold, _avoidTrianglesThreshold, _loadLastModel,
    _onlyGrowthInHeight, _VLRBezierSurface, _VLRDataPointSurface );
  
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

MyModel::~MyModel()
{
  delete _divSetting;
}

//----------------------------------------------------------------

void MyModel::renderImage()
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

void MyModel::restartModel()
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

void MyModel::setMovieParameters()
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

void MyModel::setStatus()
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

void MyModel::updateFromOld( const cell& cl, const cell& cr, const cell& c,
                    const MyTissue::division_data& ddata, MyTissue& )
{
  double angle = ModelUtils::determineDivisionAngle( ddata );
  DivisionType::type divType;
  
  // side model
  if( _surfaceType == 0 )
    divType = ModelUtils::determineSideModelDivisionType( ddata, _angleThreshold );
  // radial model
  else if( _surfaceType == 1 )
    divType = ModelUtils::determineRadialModelDivisionType( ddata, _angleThreshold, c->center );
  
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
  std::vector<cell> daughterCells;
  daughterCells.push_back( cl );
  daughterCells.push_back( cr );
  
  ModelExporter::exportDivisionDaughterProperties( _divisionFileName,
                                                   daughterCells,
                                                   c->divType,
                                                   _angleThreshold,
                                                   _divOccurrences,
                                                   false );
}

//----------------------------------------------------------------

void MyModel::checkLayerAppearance( const cell &c )
{
  auto iter = _firstLayerAppearances.find( c->layerValue );
  
  // if the layer was not already inserted, it is definitely the
  // first one that will be created by a periclinal division
  if( iter == _firstLayerAppearances.end() )
    _firstLayerAppearances.insert( std::make_pair(
    c->layerValue, std::make_pair( c->timeStep, c->divisionSequence ) ) );
}

//----------------------------------------------------------------

void MyModel::setCellProperties( const cell &c, const cell &parentCell )
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

void MyModel::setLayerValues( const cell& cl, const cell& cr, const cell& c,
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

void MyModel::setCellFileSequence( const cell& cl, const cell& cr, const cell& c,
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
      cl->radialDivisionSequence = c->radialDivisionSequence + "0";
      cr->radialDivisionSequence = c->radialDivisionSequence + "1";
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
      cl->radialDivisionSequence = c->radialDivisionSequence + "1";
      cr->radialDivisionSequence = c->radialDivisionSequence + "0";
    }
  }
  else
  {
    cl->cellFileSequence = c->cellFileSequence;
    cl->cellFileColoringIndex = c->cellFileColoringIndex;
    cr->cellFileSequence = c->cellFileSequence;
    cr->cellFileColoringIndex = c->cellFileColoringIndex;
    cl->radialDivisionSequence = c->radialDivisionSequence;
    cr->radialDivisionSequence = c->radialDivisionSequence;
  }
}

//----------------------------------------------------------------

void MyModel::step()
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
  
//   if( curTime == 200 )
//   {
//     // check different division sequences and store them
//     this->storeDivisionSequences();
//     this->stopModel();
//   }
  
//   if( curTime == maxTime && _avoidTrianglesThreshold <= 0.00001 )
//   {
//     ModelExporter::exportPropabilityDistribution( 
//     "/tmp/propabilityDistribution.csv",
//     _probValues, _lengths, _choices );
//   }
  
//   if( curTime == maxTime && !_loadLastModel )
//   {
//     std::string fileName = "/tmp/modelDataProp";
//     ModelUtils::appendCellSituationType( fileName, _initialSituationType );
//     if( _highOrderPattern != 0 )
//       fileName += "_HOP" + std::to_string( _highOrderPattern );
//     
//     fileName += ".csv";
//     this->exportDataProperties( fileName );
//     
//     unsigned int counter = 0;
//     std::vector<unsigned int> cycleCounter;
//     cycleCounter.resize( 10, 0 );
//     forall(const cell& c, T.C)
//     {
//       if( c->periCycle < cycleCounter.size() )
//         cycleCounter.at( c->periCycle )++;
//     }
//     
// //       for( std::size_t c=0; c < cycleCounter.size(); c++ )
// //         std::cout << cycleCounter.at(c) << " cells in cycle " << c << std::endl;
//   }
  
  this->renderImage();
  
  // side model using loop
  if( _useLoop && curTime >= maxTime && _surfaceType == 0 )
  {      
    // before restart the model, save the model information such
    // as the division occurrences as well as the layering
    this->updateLayerCount();
      
    // export information of periclinal divisions
    std::string filename = "/tmp/SideModelDivProperties";
    
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
  
  // radial model using loop
  if( _useLoop && curTime >= 200 && _surfaceType == 1 )
  {   
    this->storeDivisionSequences();
      
    // export information of periclinal divisions
    std::string filename = "/tmp/RadialModelDivProperties";
    
    ModelUtils::appendCellSituationType( filename, _initialSituationType );
    
    if( _useAlternativeDT )
    {
      filename += _divisionType;
      if( _divisionType != "Besson-Dumais" )
        filename += std::to_string( (unsigned int)_avoidTrianglesThreshold );
    }
    else
      filename += "ShortestWall";
    
    filename += ".csv";
              
    // export all division sequence information before a radial division
    // of the maximum amount of loops is reached
    if( _loopCounter == _amountLoops )
      ModelExporter::exportDivisionSequences( filename, _divisionSequences );
    
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
    this->step_divisions();
    
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
  
void MyModel::step_divisions()
{
  std::size_t curTime = _initSurface.getTime();
  
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
      if( curTime - _timeFourCellStage > _timeDelay )
        wait = false;
      else
        wait = true;
    }
    
    if( T.C.size() == 6 &&
        ( _initialSituationType == 2 || _initialSituationType == 3 ) )
    {
      if( curTime - _timeSixCellStage > _timeDelay )
        wait = false;
      else
        wait = true;
    }
    
    // at first determine the min and max values for each cell
    forall(const cell& c, T.C)
      ModelUtils::determineXMinMax( c, T );
    
    // wait for the four-cell stage some time steps such that the
    // the future divisions are not occurring too fast
    // and only update the center of cells
    if( T.C.size() == 4 && _initialSituationType == 1 && wait )
    {
      forall(const cell& c, T.C)
        _divSetting->setCellDivisionSettings( c, curTime );
    }
    // force the initial start of the VLR
    else if( T.C.size() < 6 && _initialSituationType != 0 )
    {
      forall(const cell& c, T.C)
      {
        _divSetting->setCellDivisionSettings( c, curTime, true );
        
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
          
          _timeFourCellStage = curTime;
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
              _timeSixCellStage = curTime;
              
            to_divide.push_back(c);
          }
        }
      }
    }
    // wait for the six-cell stage some time steps such that the
    // the future divisions are not occurring too fast
    // and only update the center of cells
    else if( T.C.size() == 6 && _initialSituationType != 0 && wait )
    {
      forall(const cell& c, T.C)
        _divSetting->setCellDivisionSettings( c, curTime );
    }
    else
    {      
      forall(const cell& c, T.C)
      {
        _divSetting->setCellDivisionSettings( c, curTime );
        
        double a = c->area;
        double l = c->longestWallLength;
        if( _useAreaRatio && _useWallRatio && c->id > areaRatioStart )
        {
          // divide cells if their area size has increased by a certain percentage amount
          double initialArea = c->initialArea;
          initialArea += initialArea*_divisionAreaRatio;
          // divide cell if its wall length has increased by a certain percentage amount
          double initialLongestLength = c->initialLongestWallLength;
          initialLongestLength += initialLongestLength*_divisionWallRatio;
          
          if( a > initialArea || l > initialLongestLength )
            to_divide.push_back(c);
        } 
        // apply a division if the cells area exceeds a certain threshold or area ratio
        else if( _useCombinedAreaRatio && c->id > areaRatioStart )
        {
          // divide cells if their area size has increased by a certain percentage amount
          double initialArea = c->initialArea;
          initialArea += initialArea*_divisionAreaRatio;
          if( a > initialArea && a > _divisionArea )
            to_divide.push_back(c);
        }
        // only apply the division based on ratio with at least areaRatioStart cells
        else if( _useAreaRatio && c->id > areaRatioStart )
        {
          // divide cells if their area size has increased by a certain percentage amount
          double initialArea = c->initialArea;
          initialArea += initialArea*_divisionAreaRatio;
          if( a > initialArea )
            to_divide.push_back(c);
        }
        // only apply the division based on ratio with at least areaRatioStart cells
        else if( _useWallRatio && c->id > areaRatioStart )
        {
          // divide cell if its wall length has increased by a certain percentage amount
          double initialLongestLength = c->longestWallLength;
          initialLongestLength += initialLongestLength*_divisionWallRatio;
          if( l > initialLongestLength )
            to_divide.push_back(c);
        }
        else
        {
          if( a > _divisionArea )
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
          MyTissue::division_data ddata = _divSetting->setDivisionPoints( c );
          this->applyDivision( c, ddata );
          T.divideCell( c, ddata );
        }
      }
      // set the division properties for the second division
      // anticlinal division resulting in six cells
      // **********************************
      // *         *  *  *  *  *          *
      // *         *  *  *  *  *          *
      // *         *  *  *  *  *          *
      // **********************************
      else if( T.C.size() > 3 && T.C.size() < 6 && _initialSituationType == 3 )
      {
        if( c->innerCell )
        {
          c->angle = 90.;
          MyTissue::division_data ddata = _divSetting->setDivisionPoints( c );
          this->applyDivision( c, ddata );
          T.divideCell( c, ddata );
        }
      }
      else if( _divisionType == "Decussation" && c->id > areaRatioStart &&
               _useAlternativeDT )
      {
        // if true then the next division is perpendicular to the last one
        // else the division is collinear to the previous division direction
        if( ModelUtils::setNextDecussationDivision( _probabilityOfDecussationDivision ) )
        {          
          double noise = _divSetting->generateNoiseInRange( -10., 10., 1000 );
          c->angle = fmod( c->angle + 90., 360. ) + noise;
        }
        
        MyTissue::division_data ddata = _divSetting->setDivisionPoints( c );
        this->applyDivision( c, ddata );
        T.divideCell( c, ddata );
      }
      else if( _divisionType == "PerToGrowth" && c->id > areaRatioStart &&
               _useAlternativeDT )
      {
        // if true then the next division is perpendicular to the principal
        // component growth of the center deformation since the cell was born
        // determine principal growth direction
        if( c->centerPos.size() > 2 )
          c->principalGrowthDir = ModelUtils::determineLongestPCGrowth( c->centerPos );
        else
          std::cout << "Too few center positions!" << std::endl;
        
        Point3d xaxisDir = Point3d( 1., 0., 0. );
        double noise = _divSetting->generateNoiseInRange( -10., 10., 1000 );
        c->angle = 180./M_PI * acos( c->principalGrowthDir*xaxisDir );
        if( ModelUtils::determineSlope( c->principalGrowthDir ) > 0 )
        {
          // only add the property that it is perpendicular to the PCA dir
          c->angle += 90. + noise;
        }
        else
        {
          // add the property AND handle the missing degrees for having a negative slope
          if( c->principalGrowthDir.i() > 0 )
            c->angle = 180. - c->angle;
          
          c->angle += 90. + noise;
        }
        // update the division angle
        MyTissue::division_data ddata = _divSetting->setDivisionPoints( c );
        this->applyDivision( c, ddata );
        
        // store the line positions for later drawing
        Point3d midDiv = (ddata.pu+ddata.pv)/2.;
        _pcLines.push_back( std::make_pair( midDiv - 5.*c->principalGrowthDir,
                                            midDiv + 5.*c->principalGrowthDir ) );
        
        T.divideCell( c, ddata );
      }
      else if( _divisionType == "Energy" && c->id > areaRatioStart &&
               _useAlternativeDT )
      {
        // if true then all division planes are determined preserving
        // almost non-triangle cells with same area sizes and the shortest cell wall
        // has the lowest energy while the longest one has the highest energy;
        // probability values are then assigned like this:
        // lowest energy -> high probability
        // highest energy -> low probability
        bool empty = false;
        MyTissue::division_data ddata = _divSetting->getEnergyDivisionData( c, empty );
        // if the division data is empty then just use the default division
        // rule (e.g. ShortestWall)
        if( empty )
          T.divideCell( c );
        else
        {
          this->applyDivision( c, ddata );
          T.divideCell( c, ddata );
        }
      }
      else if( _divisionType == "Besson-Dumais" && c->id > areaRatioStart &&
               _useAlternativeDT )
      {
        // the division plane is chosen based on Gibbs measure as described
        // in the paper of Besson and Dumais, 2011 (Universal rule for the 
        // symmetric division of plant cells)
        bool empty = false;
        MyTissue::division_data ddata = _divSetting->getBessonDumaisDivisionData( c, empty );
        // if the division data is empty then just use the default division
        // rule (e.g. ShortestWall)
        if( empty )
          T.divideCell( c );
        else
        {
          this->applyDivision( c, ddata );
          T.divideCell( c, ddata );
        }
      }
      else if( _divisionType == "RandomEqualAreas" && c->id > areaRatioStart &&
               _useAlternativeDT )
      {
        // if true then all division planes are determined randomly preserving
        // cells with almost same area sizes
        bool empty = false;
        MyTissue::division_data ddata = _divSetting->getRandomDivisionData( c, empty, true );
        // if the division data is empty then just use the default division
        // rule (e.g. ShortestWall)
        if( empty )
          T.divideCell( c );
        else
        {
          this->applyDivision( c, ddata );
          T.divideCell( c, ddata );
        }
      }
      else if( _divisionType == "Random" && c->id > areaRatioStart && _useAlternativeDT )
      {
        // if true then all division planes are determined randomly
        bool empty = false;
        MyTissue::division_data ddata = _divSetting->getRandomDivisionData( c, empty, false );
        // if the division data is empty then just use the default division
        // rule (e.g. ShortestWall)
        if( empty )
          T.divideCell( c );
        else
        {
          this->applyDivision( c, ddata );
          T.divideCell( c, ddata );
        }
      }
      else
        T.divideCell(c);
    }

    //return !to_divide.empty();  
}

//----------------------------------------------------------------

void MyModel::applyDivision( const cell &c,
                             MyTissue::division_data &ddata )
{
  vvcomplex::testDivisionOnVertices(c, ddata, T, 0.01);
  
  // apply cell pinching
  tissue::CellPinchingParams params;
  params.cellPinch = _cellPinch;
  params.cellMaxPinch = _cellMaxPinch;
  tissue::cellPinching( c, T, ddata, params );
}

//----------------------------------------------------------------

void MyModel::updateLayerCount()
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

void MyModel::storeDivisionSequences()
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

void MyModel::initDraw( Viewer* viewer )
{
  viewer->setSceneBoundingBox(Vec(-10.0, -10.0, -1.0), Vec(10.0, 10.0, 1.0));
}

//----------------------------------------------------------------

void MyModel::preDraw()
{
  util::Palette::Color bg = _palette.getColor(bgColor);
  glClearColor(bg.r(), bg.g(), bg.b(), 1);
  T.preDraw();
}

//----------------------------------------------------------------

void MyModel::draw( Viewer* viewer )
{
  // center sphere
  GraphicsClass::drawSphere( Point3d( 0., 0., 0. ), 8., 20., 20.,
                              Colorf( 0.5, 0.5, 0.5, 1. ), quadratic );
  
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
Colorf MyModel::cellColor(const cell& c)
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
      Colorf col = _palette.getColor(
        _cellFileColorIndex.at( ( std::abs(c->cellFile) )%_cellFileColorIndex.size() ) );
      
      if( c->cellCycle > 0 )
      {
        Point3d hsv = ModelUtils::transformRGBToHSV( Point3d( col.r(), col.g(), col.b() ), false );
        
        if( c->radialDivisionSequence.back() == '1' )
        {
          hsv[1] += 0.2 * c->radialDivisionSequence.size();
          hsv[2] += 0.2 * c->radialDivisionSequence.size();
        }
        else
        {
          hsv[1] -= 0.2 * c->radialDivisionSequence.size();
          hsv[2] -= 0.2 * c->radialDivisionSequence.size();
        }
        
        Point3d rgb = ModelUtils::transformHSVToRGB( hsv );
        
        return Colorf( rgb[0], rgb[1], rgb[2], 1. );
      }
      else
        return col;
  }
}

//----------------------------------------------------------------

// Method needed by the tissue
void MyModel::setPosition(const cell& c, const Point3d& p)
{
  if( SURFACETYPE == 0 )
    _VLRBezierSurface.setPos( c->sp, p );
  else
    _VLRDataPointSurface.setPos( c->sp, p );
}

//----------------------------------------------------------------

// Method needed by the tissue
void MyModel::setPosition(const junction& j, const Point3d& p)
{ 
  if( SURFACETYPE == 0 )
    _VLRBezierSurface.setPos( j->sp, p );
  else
    _VLRDataPointSurface.setPos( j->sp, p );
}

//----------------------------------------------------------------

void MyModel::exportDataProperties( const std::string &filename )
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
  
  std::vector<std::size_t> randomChoices = _divSetting->getRandomChoices();
  
  out << randomChoices.size() << " ";
  
  for( std::size_t c = 0; c < randomChoices.size(); c++ )
    out << randomChoices.at(c) << " ";
  
  out << "\n";
  
  out.close();
}

//----------------------------------------------------------------

void MyModel::readDataProperties( const std::string &filename )
{
  // TODO
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