#include "ModelHeader.h"
#include "ModelExporter.h"
#include "ModelUtils.h"

class MyModel : public Model 
{
  Q_OBJECT
public: 
  util::Parms parms;
  util::Palette palette;
	Surface lateralRoot;
  MyTissue T;

  double dt;
  double divisionArea;
  bool useAreaRatio;
  double divisionAreaRatio;
  bool useWallRatio;
  double divisionWallRatio;
  int bgColor;
  int stepPerView;
  std::size_t _initialCellNumber;
  std::string _initialCellsOfRealData;
  bool exportLineage;
  bool exportDivisionProperties;
  double probabilityOfDecussationDivision;
  bool useDecussationDivision;
  double _angleThreshold;
  
  std::size_t _idCounter;
  std::size_t _time;
  std::size_t _maxTime;
  std::size_t _jId;
  std::string _lineageFileName;
  std::string _cellWallsFileName;
  std::string _divisionFileName;
  std::vector<std::size_t> _layerColorIndex;
  double _cellPinch;
  double _cellMaxPinch;
  std::size_t _cellColoringType;
  // defines the number of vertices for each wall of the init cells
  // #inner vertices of each wall = _cellSubdivisionLevel - 1
  std::size_t _cellSubdivisionLevel;
	
	std::pair<std::size_t, std::size_t> _divOccurrences;
  
  //----------------------------------------------------------------
  
  void readParms()
  {
    // read the parameters here
    parms("Main", "Dt", dt);
    parms("Main", "InitialCellNumber", _initialCellNumber);
    parms("Main", "InitialCellsOfRealData", _initialCellsOfRealData);
    parms("Main", "SubDivisionLevelOfCells", _cellSubdivisionLevel);
    parms("Main", "ExportLineage", exportLineage);
    parms("Main", "ExportDivisionProperties", exportDivisionProperties);

    parms("View", "StepPerView", stepPerView);
    parms("View", "BackgroundColor", bgColor);

    parms( "Division", "DivisionArea", divisionArea);
    parms( "Division", "UseAreaRatio", useAreaRatio);
    parms( "Division", "DivisionAreaRatio", divisionAreaRatio);
    parms( "Division", "UseWallRatio", useWallRatio);
    parms( "Division", "DivisionWallRatio", divisionWallRatio);
    parms( "Division", "UseDecussationDivision", useDecussationDivision );
    parms( "Division", "ProbabilityOfDecussationDivision", probabilityOfDecussationDivision );
    parms( "Division", "DivisionAngleThreshold", _angleThreshold );
    parms( "Division", "CellColoringType", _cellColoringType );

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
  
  MyModel(QObject *parent) : Model(parent), parms("view.v"), palette("pal.map"), 
    lateralRoot(parms, "Surface"), T(palette, this),
    _idCounter(1), _time(1), _maxTime(300), _jId( 0 ),
    _lineageFileName( "/tmp/model.csv" ),
    _cellWallsFileName( "/tmp/modelCellWalls.csv" ),
    _divisionFileName( "/tmp/divisionPropertiesModel.csv" ),
    _divOccurrences( std::make_pair( 0, 0 ) )
  {
    readParms();
    // Registering the configuration files
    registerFile("pal.map");
    registerFile("view.v");

    // single layer assignment
    //for( std::size_t l = 9; l < 14; l++ )
    // multiple layer assignment for each new daughter cell
    for( std::size_t l = 14; l < 45; l++ )
      _layerColorIndex.push_back( l );
    
    if( exportLineage )
    {
      ModelExporter::initExportFile( _lineageFileName );
      ModelExporter::initCellWallFile( _cellWallsFileName );
    }
    
    if( exportDivisionProperties )
      ModelExporter::initDivisionDaughterFile( _divisionFileName );
    
    if( _cellSubdivisionLevel < 1 )
      _cellSubdivisionLevel = 1;
    
    lateralRoot.GrowStep(0);
    idPairSet sharedJunctions;
    
    // special cases of number of cells at the beginning
    if( _initialCellsOfRealData == "none" )
    {
      if( _initialCellNumber == 1 )
      {
        this->generateCell( std::make_pair( 0., 0. ),
                            std::make_pair( 1., 1. ),
                            1, sharedJunctions );
      }
      else if( _initialCellNumber == 2 )
      {
        this->generateCell( std::make_pair( 0., 0. ),
                            std::make_pair( 0.5, 1. ),
                            1, sharedJunctions );
        
        // insert ids of shared junctions
        sharedJunctions.insert( std::make_pair( 5, 2 ) );
        sharedJunctions.insert( std::make_pair( 4, 3 ) );
        
        this->generateCell( std::make_pair( 0.5, 0. ),
                            std::make_pair( 0.5, 1. ),
                            2, sharedJunctions );
      }
      // TODO
      else if( _initialCellNumber == 8 )
        this->initLateralRoot();
      else
      {
        std::cerr << "Selected number of cells at the beginning is not implemented yet!" << std::endl;
      }
    }
    // else constellation of cells depending on the real data
    else
      this->initLateralRoot( _initialCellsOfRealData );
    
    setStatus();
  }
 
  //----------------------------------------------------------------

  void initLateralRoot( const std::string &dataset )
	{
    std::cout << "Init constellation of data: " << dataset << std::endl;
    idPairSet sharedJunctions;
    
		if( dataset == "120830" )
    {
      this->generateCell( std::make_pair( 0., 0. ),
                          std::make_pair( 0.5, 1. ),
                          1, sharedJunctions );
      
      // insert ids of shared junctions
      sharedJunctions.insert( std::make_pair( 5, 2 ) );
      sharedJunctions.insert( std::make_pair( 4, 3 ) );
      
      this->generateCell( std::make_pair( 0.5, 0. ),
                          std::make_pair( 0.5, 1. ),
                          2, sharedJunctions );
    }
		else if( dataset == "121204" )
    {
      std::size_t lineageCounter = 1;
   
      for( std::size_t c = 0; c < 2; c++ )
      {
        double u = 0. + c*2./5.;
        double v = 0.;
        double length;
        if( c == 0 )
          length = 2./5.;
        else
          length = 3./5.;
          
        if( c == 1 )
        {
          // insert ids of shared junctions
          sharedJunctions.insert( std::make_pair( 5, 2 ) );
          sharedJunctions.insert( std::make_pair( 4, 3 ) );
        }
        
        this->generateCell( std::make_pair( u, v ),
                            std::make_pair( length, 1. ),
                            lineageCounter, sharedJunctions );
        
        lineageCounter++;
      }
    }
    // TODO
    else if( dataset == "121211" )
    {
      std::size_t lineageCounter = 1;
   
      // outer bigger cells
      for( std::size_t c = 0; c < 2; c++ )
      {
        double u = 0. + c*2./3.;
        double v = 0.;
        this->generateCell( std::make_pair( u, v ),
                            std::make_pair( 1./3., 1. ),
                            lineageCounter, sharedJunctions );
        
        lineageCounter++;
      }
      
      // inner smaller cells
      for( std::size_t c = 0; c < 3; c++ )
      {
        double u = 1./3. + c*1./9.;
        double v = 0.;
        this->generateCell( std::make_pair( u, v ),
                            std::make_pair( 1./9., 1. ),
                            lineageCounter, sharedJunctions );
        
        lineageCounter++;
      }
    }
    // TODO
    else if( dataset == "130508" )
      this->generateCell( std::make_pair( 0., 0. ),
                          std::make_pair( 1., 1. ),
                          1, sharedJunctions );
    // TODO
    else if( dataset == "130607" )
    {
      std::size_t lineageCounter = 1;
   
      for( std::size_t c = 0; c < 2; c++ )
      {
        double u = 0. + c*2./5.;
        double v = 0.;
        double length;
        if( c == 0 )
          length = 2./5.;
        else
          length = 3./5.;
          
        this->generateCell( std::make_pair( u, v ),
                            std::make_pair( length, 1. ),
                            lineageCounter, sharedJunctions );
        
        lineageCounter++;
      }
    }
    else
      std::cerr << "Selected data set name is not supported!" << std::endl;
        
	}
  
	//----------------------------------------------------------------
	
  void initLateralRoot()
  {
    // render eight cells at the beginning for which each pair shares an
    // area ratio of 2:1. The bigger cell will always be located at the left
    // and right boundary
    
    idPairSet sharedJunctions;
    std::size_t lineageCounter = 1;
    
    // init the inner smaller cells
    for( std::size_t w = 0; w < 2; w++ )
      for( std::size_t h = 0; h < 2; h++ )
      {
        double u = 1./3. + w*1./6.;
        double v = 0. + h*1./2.;
        this->generateCell( std::make_pair( u, v ),
                            std::make_pair( 1./6., 1./2. ),
                            lineageCounter, sharedJunctions );
        
        lineageCounter++;
      }
    
    // then the outer bigger cells
    for( std::size_t w = 0; w < 2; w++ )
      for( std::size_t h = 0; h < 2; h++ )
      {
        double u = 0. + w*2./3.;
        double v = 0. + h*1./2.;
        this->generateCell( std::make_pair( u, v ),
                            std::make_pair( 1./3., 1./2. ),
                            lineageCounter, sharedJunctions );
        
        lineageCounter++;
      }
  }
  //----------------------------------------------------------------
  
  // generate a cell starting at the bottom left with coord u and v with different lengths
  void generateCell( const std::pair<double, double> &start,
                     const std::pair<double, double> &length,
                     const std::size_t treeId,
                     const idPairSet &sharedJunctions )
  {
    // set of junctions for the cell
    std::vector<junction> vs;
    std::vector< std::pair<double,double> > vertices;
    vertices.resize( 4 * _cellSubdivisionLevel );
    
    for( std::size_t w = 0; w < vertices.size(); w++ )
    {
      junction j;
      unsigned int uiw = (unsigned int)(w/_cellSubdivisionLevel);
      double cellLength;
      double curLength = (double)w/(double)_cellSubdivisionLevel - (double)uiw;
      double u,v;
      switch(uiw)
      {
        case 0:
          cellLength = length.second;
          u = start.first;
          v = start.second + curLength*cellLength;
          break;
        case 1:
          cellLength = length.first;
          u = start.first + curLength*cellLength;
          v = start.second + length.second;
          break;
        case 2:
          cellLength = length.second;
          u = start.first + length.first;
          v = start.second + length.second - curLength*cellLength;
          break;
        case 3:
          cellLength = length.first;
          u = start.first + length.first - curLength*cellLength;
          v = start.second;
          break;
        default: u = v = 0.; break;
      }
      
      //std::cout << "id: " << _jId << std::endl;
      //std::cout << u << " " << v << std::endl;
      vertices.at( w ) = std::make_pair( u, v );
      j->id = _jId;
      lateralRoot.InitPoint( j->sp, u, v );
      junctionAlreadyShared( j->id, j, sharedJunctions );
      
      vs.push_back(j);
      
      _jId++;
    }
    
    cell c;
    c->treeId = treeId;
    c->id = _idCounter;
    c->parentId = _idCounter;
    c->timeStep = _time;
    c->previousAngle = 0.;
    c->angle = 0.;
    c->previousDivDir = Point3d( 0., 0., 0. );
    c->divDir = Point3d( 0., 0., 0. );
    c->divType = DivisionType::NONE;
    c->layerValue = 1;
    c->cellCycle = 0;
    T.addCell( c, vs );
    
    std::vector<Point3d> polygon;
    Point3d center;
    forall(const junction& j, T.S.neighbors(c))
    {
      polygon.push_back(j->sp.Pos());
      center += j->sp.Pos();
    }
    center /= polygon.size();

    // store initial area for current cell
    c->center = center;
    c->initialArea = geometry::polygonArea(polygon);
    c->area = c->initialArea;
    c->initialLongestWallLength = ModelUtils::determineLongestWallLength( c, T );
    c->longestWallLength = c->initialLongestWallLength;
    lateralRoot.SetPoint(c->sp, c->sp, center);
    
    ModelExporter::exportLineageInformation( _lineageFileName, c, T, _time );
    
    // afterwards increment the id counter
    _idCounter++;
  }
  
  //----------------------------------------------------------------
  
  void junctionAlreadyShared( const std::size_t jId, junction &js,
                              const idPairSet &sharedJunctions )
  {
    forall( const cell& c, T.C )
    {
      forall(const junction& j, T.S.neighbors(c))
      {
        idPairSet::const_iterator iter =
        sharedJunctions.find( std::make_pair( jId, j->id ) );
        if( iter != sharedJunctions.end() )
        {
          js = j;
          std::cout << "true" << std::endl;
          return;
        }
      }
    }
  }
  
  //----------------------------------------------------------------
  
  bool setNextDecussationDivision()
  {
    // here we radomly decide which kind of division the current cell
    // will do based on a probability value given as a parameter
    srand( _time + _idCounter + time(NULL) );
    
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
      const junction& jn = T.S.nextTo(c, j);
      const Point3d& jpos = j->sp.Pos();
      const Point3d& jnpos = jn->sp.Pos();
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
		//double time = lateralRoot.GetTime();
    if( _time > _maxTime )
      _time = _maxTime;
    
    QString status = QString( "# Vertices: %1 \t "
                              "# Cells: %2 \t"
                              "Time step: %3 \t ").
                              arg(T.W.size()).
                              arg(T.C.size()).
                              arg(_time);
    
    if( useAreaRatio )
      status += QString( "Area divison ratio: %1 \t" ).arg(divisionAreaRatio);
    
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
    
		// set daughter cell properties
    // left cell
    std::vector<Point3d> polygon;
    Point3d center;
    forall(const junction& j, T.S.neighbors(cl))
    {
      polygon.push_back(j->sp.Pos());
      center += j->sp.Pos();
    }
    
    center /= polygon.size();
    cl->center = center;
    cl->initialArea = geometry::polygonArea(polygon);
    cl->area = cl->initialArea;
    cl->initialLongestWallLength = ModelUtils::determineLongestWallLength( cl, T );
    cl->longestWallLength = cl->initialLongestWallLength;
    
		cl->treeId = c->treeId;
    cl->id = _idCounter;
    cl->parentId = c->id;
    cl->timeStep = _time;
    cl->previousAngle = c->angle;
    cl->angle = angle;
    cl->previousDivDir = c->divDir;
    cl->divDir = cl->center - c->center;
    cl->cellCycle = c->cellCycle+1;
    _idCounter++;

    
    // right cell
    polygon.clear();
    center = Point3d( 0., 0., 0. );
    forall(const junction& j, T.S.neighbors(cr))
    {
      polygon.push_back(j->sp.Pos());
      center += j->sp.Pos();
    }
    
    center /= polygon.size();
    cr->center = center;
    cr->initialArea = geometry::polygonArea(polygon);
    cr->area = cr->initialArea;
    cr->initialLongestWallLength = ModelUtils::determineLongestWallLength( cr, T );
    cr->longestWallLength = cr->initialLongestWallLength;
    
    cr->treeId = c->treeId;
    cr->id = _idCounter;
    cr->parentId = c->id;
    cr->timeStep = _time;
    cr->previousAngle = c->angle;
    cr->previousDivDir = c->divDir;
    cr->divDir = cr->center - c->center;
    cr->cellCycle = c->cellCycle+1;
    _idCounter++;
        
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
    ModelExporter::exportDivisionDaughterProperties( _divisionFileName,
                                                     cl,
                                                     cr,
                                                     ddata,
                                                     _angleThreshold,
																										 _divOccurrences );
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
  
  bool step_divisions()
  {
    // for which number of cells should the division area ratio check apply
    // this is required since at the beginning only few divisions occur due
    // to the small increasing of area based on the inital area size
    std::size_t areaRatioStart;
    if( T.C.size() > 1 )
      areaRatioStart = 0;
    else
      areaRatioStart = 1;
    
    // Find cells to be divided
    std::list<cell> to_divide;
    forall(const cell& c, T.C)
    { 
      // update center position
      std::vector<Point3d> polygon;
      Point3d center;
      forall(const junction& j, T.S.neighbors(c))
      {
        polygon.push_back(j->sp.Pos());
        center += j->sp.Pos();
      }
      center /= polygon.size();
      lateralRoot.SetPoint(c->sp, c->sp, center);
      c->center = center;
      c->area = geometry::polygonArea(polygon);
      c->timeStep = _time;
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
      // only apply the division based on ratio with at least eight cells
      else if( useAreaRatio && c->id > areaRatioStart )
      {
        // divide cells if their area size has increased by a certain percentage amount
        double initialArea = c->initialArea;
        initialArea += initialArea*divisionAreaRatio;
        if( a > initialArea )
          to_divide.push_back(c);
      }
      // only apply the division based on ratio with at least eight cells
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
      else
        T.divideCell(c);
    }

    return !to_divide.empty();
  }
  
  //----------------------------------------------------------------
  
  void step_cellWalls()
  {
    if( _time <= _maxTime + 1 )
    {
      if( exportLineage )
      {
        forall(const cell& c, T.C)
        {
          if( c->timeStep > _maxTime )
            break;
            
          ModelExporter::exportCellWalls( _cellWallsFileName, c, T );
        }
      }
    }
  }
  
  //----------------------------------------------------------------
  
  void step_tracking()
  {
    if( _time <= _maxTime )
    {
      if( exportLineage )
      {
        forall(const cell& c, T.C)
          ModelExporter::exportLineageInformation( _lineageFileName, c, T, _time );
      }
    }
  }
  
  //----------------------------------------------------------------
  
  void step_growth()
  {
	  lateralRoot.GrowStep(dt);
    forall(const junction& v, T.W)
		  lateralRoot.GetPos(v->sp);
  }

  //----------------------------------------------------------------
  
  void step()
  {
    if( _time <= _maxTime )
      _time++;
    
    for(int i = 0 ; i < stepPerView ; ++i)
    {
      this->step_cellWalls();
      this->step_divisions();
      this->step_tracking();
      this->step_growth();
    }
    this->setStatus();
		
		if( _time == _maxTime )
		{
			std::cout << "Anticlinal divisions: " << _divOccurrences.first << std::endl;
			std::cout << "Periclinal divisions: " << _divOccurrences.second << std::endl;
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
      T.drawCell(c, this->cellColor(c), this->cellColor(c)*0.7 );
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
    }
  }

  //----------------------------------------------------------------
  
	// color for contour of cells
  Colorf contourColor(const cell& c)
  {
    return palette.getColor(T.contourColor);
  }

  //----------------------------------------------------------------
  
  // Methods needed by the tissue
  Point3d position(const cell& c) const { return c->sp.Pos(); }
  Point3d position(const junction& c) const { return c->sp.Pos(); }

  void setPosition(const cell& c, const Point3d& p) { lateralRoot.SetPoint(c->sp, c->sp, p); }
  void setPosition(const junction& j, const Point3d& p) { lateralRoot.SetPoint(j->sp, j->sp, p); }

  void setPositionHint(const junction&, const junction&, const junction&, double) {}

  Point3d normal(const junction& ) const { return Point3d(0,0,1); }
  Point3d normal(const cell& ) const { return Point3d(0,0,1); }
};

#include "model.moc"
DEFINE_MODEL(MyModel);
