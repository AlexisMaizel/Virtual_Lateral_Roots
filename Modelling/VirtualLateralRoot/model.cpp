#include "ModelHeader.h"
#include "ModelExporter.h"
#include "ModelUtils.h"
#include "SurfaceClass.h"

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
  util::Parms parms;
  Surface _VLRBezierSurface;
  RealSurface _VLRDataPointSurface;
  util::Palette palette;
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
  bool _exportLineage;
  bool _exportDivisionProperties;
  double probabilityOfDecussationDivision;
  bool useDecussationDivision;
  double _angleThreshold;
  std::size_t surfaceType;
  bool _bezierGrowthSurface;
  double _surfaceScale;
  bool _useAutomaticContourPoints;
  
  std::string _lineageFileName;
  std::string _cellWallsFileName;
  std::string _divisionFileName;
  std::string _timeAgainstCellsFileName;
  std::vector<std::size_t> _layerColorIndex;
  double _cellPinch;
  double _cellMaxPinch;
  std::size_t _cellColoringType;
  std::pair<std::size_t, std::size_t> _divOccurrences;
  
  SurfaceClass _surfaceClass;
  std::size_t _lod;
  
  bool _lastStep;
  
  SurfaceType::type _sType;
  
  bool _forceInitialSituation;
  
  Point3d _initialBoundary;
  
  double _firstDivisionsAreaRatio;
  double _secondDivisionsAreaRatio;
  std::size_t _timeDelay;
  std::size_t _timeSixCellStage;
  bool _smootherCells;
  
  //----------------------------------------------------------------
  
  void readParms()
  {
    // read the parameters here
#ifndef TimeAgainstCellAnalysis
    parms("Main", "Dt", dt);
#endif
    parms("Main", "InitialCellNumber", _initialCellNumber);
    parms("Main", "InitialCellsOfRealData", _initialCellsOfRealData);
    parms("Main", "SubDivisionLevelOfCells", _lod);
    parms("Main", "ExportLineage", _exportLineage);
    parms("Main", "ExportDivisionProperties", _exportDivisionProperties);
    parms("Main", "SurfaceType", surfaceType);
    parms("Main", "BezierGrowthSurface", _bezierGrowthSurface);
    parms("Main", "SurfaceScale", _surfaceScale);
    parms("Main", "UseAutomaticContourPoints", _useAutomaticContourPoints );
    parms("Main", "ForceInitialSituation", _forceInitialSituation );
    parms("Main", "SmootherCells", _smootherCells );

    parms("View", "StepPerView", stepPerView);
    parms("View", "BackgroundColor", bgColor);

    parms( "Division", "DivisionArea", divisionArea);
    parms( "Division", "DivisionAreaRatio", divisionAreaRatio);
    parms( "Division", "UseAreaRatio", useAreaRatio);
    parms( "Division", "UseCombinedAreaRatio", useCombinedAreaRatio);
    parms( "Division", "UseWallRatio", useWallRatio);
    parms( "Division", "DivisionWallRatio", divisionWallRatio);
    parms( "Division", "UseDecussationDivision", useDecussationDivision );
    parms( "Division", "ProbabilityOfDecussationDivision", probabilityOfDecussationDivision );
    parms( "Division", "DivisionAngleThreshold", _angleThreshold );
    parms( "Division", "CellColoringType", _cellColoringType );
    
    parms( "Division", "FirstDivisionsAreaRatio", _firstDivisionsAreaRatio);
    parms( "Division", "SecondDivisionsAreaRatio", _secondDivisionsAreaRatio);
    parms( "Division", "TimeDelay", _timeDelay);

    parms( "Tissue", "CellPinch", _cellPinch );
    parms( "Tissue", "CellMaxPinch", _cellMaxPinch );
    
    T.readParms(parms, "Tissue");
    T.readViewParms(parms, "TissueView");
  }

  //----------------------------------------------------------------
  
  // Here, reread the files when they are modified
  void modifiedFiles( const std::set<std::string>& filenames )
  {
    forall(const std::string& fn, filenames)
    {
      if(fn == "pal.map")
        palette.reread();
      else if(fn == "view.v")
        readParms();
    }
  }

  //----------------------------------------------------------------
  
  MyModel(QObject *parent) : Model(parent), parms("view.v"),
    _VLRDataPointSurface( parms, "Surface" ),
    palette("pal.map"), T(palette, this),
    _divOccurrences( std::make_pair( 0, 0 ) ),
    _lastStep( false )
  {
    readParms();
    // Registering the configuration files
    registerFile("pal.map");
    registerFile("view.v");

    if( surfaceType == 0 )
      _sType = SurfaceType::BEZIER;
    else
      _sType = SurfaceType::REALDATA;
    
    std::size_t t = 300;
    if( _initialCellsOfRealData == "130508_raw" )
      t = 350;
    else if( _initialCellsOfRealData.compare( 0, 7, "Average") == 0 )
      t = 150;
    
    //std::cout << "LOD: " << _lod << std::endl;
    
    _surfaceClass.init( _lod, _lineageFileName, _exportLineage, t, _sType, _forceInitialSituation );
    
    _VLRBezierSurface.init( parms, "Surface",
                            _bezierGrowthSurface, _initialCellsOfRealData );
    
    // set name strings
    _lineageFileName = "/tmp/model" + _initialCellsOfRealData + ".csv";
    _cellWallsFileName = "/tmp/modelCellWalls" + _initialCellsOfRealData + ".csv";
    _divisionFileName = "/tmp/divisionPropertiesModel" + _initialCellsOfRealData + ".csv";
    
    // single layer assignment
    //for( std::size_t l = 9; l < 14; l++ )
    // multiple layer assignment for each new daughter cell
    for( std::size_t l = 14; l < 45; l++ )
      _layerColorIndex.push_back( l );
    
    cell dummy;
    MyTissue::division_data ddummy;
    
    if( _exportLineage )
    {
      ModelExporter::exportLineageInformation( _lineageFileName, dummy, T,
                                               0, true, surfaceType );
      ModelExporter::exportCellWalls( _cellWallsFileName, dummy, T,
                                      true, surfaceType );
    }
    
    if( _exportDivisionProperties )
    {
      std::pair<std::size_t, std::size_t> pair;
      ModelExporter::exportDivisionDaughterProperties( _divisionFileName,
                                                       dummy, dummy,
                                                       ddummy, 0.,
                                                       pair, true );
    }
    
    _timeAgainstCellsFileName = "/tmp/timeAgainstCells";
    _timeAgainstCellsFileName += _initialCellsOfRealData;
    _timeAgainstCellsFileName += ".csv";
    
    if( dt > start-steps )
      ModelExporter::exportTimeAgainstCells( _timeAgainstCellsFileName, dt, 0, true );
        
    // bezier
    if( surfaceType == 0 )
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
    if( _exportLineage )
    {
      forall(const cell& c, T.C)
      {
        ModelExporter::exportLineageInformation( _lineageFileName, c, T,
                                                 _surfaceClass.getTime(),
                                                 false, surfaceType );
        
        ModelExporter::exportCellWalls( _cellWallsFileName, c,
                                        T, false, surfaceType );
      }
    }
    
    // determine the initial boundary between two founder cells
    _initialBoundary = Point3d( 0., 0., 0. );
    if( _forceInitialSituation )
    {
      Point3d max = Point3d( -5000., -5000., 0. );
      // consider the right most position
      forall(const cell& c, T.C)
      {
        if( c->id == 1 )
        {
          forall(const junction& j, T.S.neighbors(c))
          {
            Point3d pos;
            
            if( surfaceType == 1 )
              pos = Point3d( j->tp.Pos().i(), j->tp.Pos().j(), 0. );
            else
              pos = j->sp.Pos();
            
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
  
  bool setNextDecussationDivision()
  {
    // here we randomly decide which kind of division the current cell
    // will do based on a probability value given as a parameter
    srand( _surfaceClass.getTime() + _surfaceClass.getCellID() + time(NULL) );
    
    // generate a value between 1 and 10
    int val = rand() % 10 + 1;
    //std::cout << "val: " << val << std::endl;
    
    if( val <= (int)((1.-probabilityOfDecussationDivision)*10.) )
      return true;
    else
      return false;
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
      if( surfaceType == 0 )
      {
        jpos = j->sp.Pos();
        jnpos = jn->sp.Pos();
      }
      else
      {
        jpos = Point3d( j->tp.Pos().i(), j->tp.Pos().j(), 0. );
        jnpos = Point3d( jn->tp.Pos().i(), jn->tp.Pos().j(), 0. );
      }
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
  
  void setStatus()
  {
    //double time = _VLRBezierSurface.GetTime();
    std::size_t time = _surfaceClass.getTime();
//     if( time > _surfaceClass.getMaxTime() )
//       time = _surfaceClass.getMaxTime();
    
    QString status = QString( "Vertices: %1 \t "
                              "Cells: %2 \t").
                              arg(T.W.size()).
                              arg(T.C.size());
    
    if( surfaceType == 1 )
    {
      std::size_t timeStep = _VLRDataPointSurface.getCurTimeStep();
      status += QString( "TS: %1 \t" ).arg(timeStep);
    }
    
    status += QString( "MS: %1 \t" ).arg(time);
    status += QString( "AD: %1 \t" ).arg(_divOccurrences.first);
    status += QString( "PD: %1 \t" ).arg(_divOccurrences.second);
    
    if( useCombinedAreaRatio )
    {
      status += QString( "Area divison: %1 \t" ).arg(divisionArea);
      status += QString( "Area divison ratio: %1 \t" ).arg(divisionAreaRatio);
    }
    else if( useAreaRatio )
      status += QString( "Area divison ratio: %1 \t" ).arg(divisionAreaRatio);
    else
      status += QString( "Area divison: %1 \t" ).arg(divisionArea);
    
    if( useWallRatio )
      status += QString( "Wall divison ratio: %1 \t" ).arg(divisionWallRatio);
    
    if( useDecussationDivision )
      status += QString( "Decussation propability: %1" ).arg(probabilityOfDecussationDivision);
    
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
    
    // set cell properties for left cell
    this->setCellProperties( cl, c );
    
    // set cell properties for right cell
    this->setCellProperties( cr, c );
    
    // for the forced situation check which of the four cells
    // are the inner ones
    if( _forceInitialSituation && T.C.size() > 2 && T.C.size() < 5 )
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
    if( _exportDivisionProperties )
    {
      ModelExporter::exportDivisionDaughterProperties( _divisionFileName,
                                                       cl,
                                                       cr,
                                                       ddata,
                                                       _angleThreshold,
                                                       _divOccurrences,
                                                       false );
    }
  }

  //----------------------------------------------------------------
  
  void setCellProperties( const cell &c, const cell &parentCell )
  {
    Point3d center;
    double area;
    this->setCellCenter( c, center, area );
    c->initialArea = area;
    c->center = center;
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
    std::vector<Point2d> polygon;
    forall(const junction& j, T.S.neighbors(c))
    {
      Point3d pos;
      
      if( surfaceType == 0 )
        pos = j->sp.Pos();
      else
        pos = Point3d( j->tp.Pos().i(), j->tp.Pos().j(), 0. );
      
      polygon.push_back( Point2d( pos.i(), pos.j() ) );
      center += pos;
    }
    
    center /= polygon.size();
    area = geometry::polygonArea(polygon);
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
      }
    }
    else
    {
      cl->layerValue = c->layerValue;
      cr->layerValue = c->layerValue;
    }
  }
  
  //----------------------------------------------------------------
  
  void getxMinMax( const cell &c, double &xMin, double &xMax )
  {
    forall(const junction& j, T.S.neighbors(c))
    {
      if( surfaceType == 0 )
      {
        if( j->sp.Pos().i() < xMin )
          xMin = j->sp.Pos().i();
        
        if( j->sp.Pos().i() > xMax )
          xMax = j->sp.Pos().i();
      }
      else
      {
        if( j->tp.Pos().i() < xMin )
          xMin = j->tp.Pos().i();
        
        if( j->tp.Pos().i() > xMax )
          xMax = j->tp.Pos().i();
      }
    }  
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
    if( T.C.size() == 6 && _forceInitialSituation )
    {
      std::size_t curTime = _surfaceClass.getTime();
      if( curTime - _timeSixCellStage > _timeDelay )
        wait = false;
      else
        wait = true;
    }
    
    // force the initial start of the VLR
    if( T.C.size() < 6 && _forceInitialSituation )
    {
      forall(const cell& c, T.C)
      {
        // update center position
        Point3d center;
        double area;
        this->setCellCenter( c, center, area );
        if( T.C.size() < 4 )
        {
          double xMin = 5000.;
          double xMax = -5000.;
          this->getxMinMax( c, xMin, xMax );
          double xLength = std::fabs( xMax - xMin );
          Point3d divPos = center;
          if( c->id == 1 )
            divPos.i() = xMin + 2.*xLength/3.; 
          else if( c->id == 2 )
            divPos.i() = xMin + 1.*xLength/3.; 
          
          if( surfaceType == 0 )
            _VLRBezierSurface.SetPoint(c->sp, c->sp, divPos);
          else
            _VLRDataPointSurface.setPos(c->tp, divPos);
        }
        else
        {
          if( surfaceType == 0 )
            _VLRBezierSurface.SetPoint(c->sp, c->sp, center);
          else
            _VLRDataPointSurface.setPos(c->tp, center);
        }
      
        c->area = area;
        c->center = center;
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
    else if( T.C.size() == 6 && _forceInitialSituation && wait )
    {
      
    }
    else
    {      
      forall(const cell& c, T.C)
      {
        // update center position
        Point3d center;
        double area;
        this->setCellCenter( c, center, area );
      
        if( surfaceType == 0 )
          _VLRBezierSurface.SetPoint(c->sp, c->sp, center);
        else
          _VLRDataPointSurface.setPos(c->tp, center);
          
        c->area = area;
        c->center = center;
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
      if( useDecussationDivision && c->id > areaRatioStart )
      {
        // if true then the next division is perpendicular to the last one
        // else the division is collinear to the previous division direction
        if( this->setNextDecussationDivision() )
          c->angle = fmod( c->angle + 90., 360. );
        
        MyTissue::division_data ddata = this->setDivisionPoints( c );
        T.divideCell( c, ddata );
      }
      // set the division properties for the second division
      // periclinal division resulting in six cells
      // **********************************
      // *         *     *     *          *
      // *         *************          *
      // *         *     *     *          *
      // **********************************
      else if( T.C.size() > 3 && T.C.size() < 6 && _forceInitialSituation )
      {
        if( c->innerCell )
        {
          c->angle = 180.;
          MyTissue::division_data ddata = this->setDivisionPoints( c );
          T.divideCell( c, ddata );
        }
      }
      else
        T.divideCell(c);
    }

    return !to_divide.empty();
  }
  
  //----------------------------------------------------------------
  
  void step_cellWalls()
  {    
    if( _exportLineage && !_lastStep )
    {
      forall(const cell& c, T.C)
        ModelExporter::exportCellWalls( _cellWallsFileName, c,
                                        T, false, surfaceType );
    }
  }
  
  //----------------------------------------------------------------
  
  void step_tracking()
  {
    if( _exportLineage && !_lastStep )
    {
      forall(const cell& c, T.C)
        ModelExporter::exportLineageInformation( _lineageFileName, c, T,
                                                  _surfaceClass.getTime(),
                                                  false, surfaceType );
    }
  }
  
  //----------------------------------------------------------------
  
  void step_growth()
  {
    if( surfaceType == 0 )
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
    
    if( surfaceType == 0 )
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
      ModelExporter::exportTimeAgainstCells( _timeAgainstCellsFileName,
                                              dt, T.C.size(), false );
      
      _lastStep = true;
      std::cout << "dt: " << dt << std::endl;
      dt -= steps;
    }
  }

  //----------------------------------------------------------------
  
  void initDraw(Viewer* viewer)
  {
    viewer->setSceneBoundingBox(Vec(-10.0, -10.0, -1.0), Vec(10.0, 10.0, 1.0));
  }

  //----------------------------------------------------------------
  
  void preDraw()
  {
    util::Palette::Color bg = palette.getColor(bgColor);
    glClearColor(bg.r(), bg.g(), bg.b(), 1);
    T.preDraw();
  }

  //----------------------------------------------------------------
  
  void postDraw()
  {
    T.postDraw();
  }

  //----------------------------------------------------------------
  
  void draw(Viewer* viewer)
  {
    forall(const cell& c, T.C)
    {
      //T.drawBorders = false;
      if( surfaceType == 0 )
      {
        if( !_bezierGrowthSurface )
          T.cellWallWidth = 0.001;
        else
          T.cellWallWidth = 0.1;
      }
      else
        T.cellWallWidth = 0.2;
          
      //T.cellWallMin = 0.0001;
      //T.strictCellWallMin = true;
      T.drawCell(c, this->cellColor(c), Colorf(this->cellColor(c)*0.3) );
    }
  }

  //----------------------------------------------------------------
  
  // color for inner cells
  Colorf cellColor(const cell& c)
  {
    switch( _cellColoringType )
    {
      // coloring based on lineage trees/ founder cells
      case 0: return palette.getColor(c->treeId);
      // coloring based on layer value
      case 1:
      default:
      if( c->layerValue-1 < _layerColorIndex.size() )
        return palette.getColor( _layerColorIndex.at(c->layerValue-1) );
      else
        return palette.getColor( _layerColorIndex.at(_layerColorIndex.size()-1) );
      case 2:
        // boundary
        if( T.border( c ) )
          return palette.getColor(1);
        // inner cell
        else
          return palette.getColor(2);
    }
  }

  //----------------------------------------------------------------
  
  // color for contour of cells
  Colorf contourColor(const cell& c)
  {
    return palette.getColor(T.contourColor);
  }

  //----------------------------------------------------------------
  
  // Method needed by the tissue
  Point3d position(const cell& c) const
  { 
    if( surfaceType == 0 )
      return c->sp.Pos();
    else
      return Point3d( c->tp.Pos().i(), c->tp.Pos().j(), 0. );
  }
  
  //----------------------------------------------------------------
  
  // Method needed by the tissue
  Point3d position(const junction& c) const
  {
    if( surfaceType == 0 )
      return c->sp.Pos();
    else
      return Point3d( c->tp.Pos().i(), c->tp.Pos().j(), 0. );
  }

  //----------------------------------------------------------------
  
  // Method needed by the tissue
  void setPosition(const cell& c, const Point3d& p)
  {
    if( surfaceType == 0 )
      _VLRBezierSurface.SetPoint(c->sp, c->sp, p);
    else
      _VLRDataPointSurface.setPos( c->tp, p );
  }
  
  //----------------------------------------------------------------
  
  // Method needed by the tissue
  void setPosition(const junction& j, const Point3d& p)
  { 
    if( surfaceType == 0 )
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
