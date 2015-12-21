#include "InitSurface.h"

#include "ModelExporter.h"

/**
  @file   InitSurface.cpp
  @brief  Contains class for setting initial surfaces of model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

// ---------------------------------------------------------------------

void InitSurface::init( const std::size_t lod,
                        const std::string &lineageFileName,
                        const bool forceInitialSituation,
                        const bool useAutomaticContourPoints,
                        const double surfaceScale,
                        const std::string &dataset )
{
  if( _lod == 0 )
    _lod = 1;
  else
    _lod = lod;
  
  _IDCounter = 1;
  _jIDCounter = 0;
  _time = 1;
  _lineageFileName = lineageFileName;
  _forceInitialSituation = forceInitialSituation;
  _useAutomaticContourPoints = useAutomaticContourPoints;
  _surfaceScale = surfaceScale;
  _dataset = dataset;
}

// ---------------------------------------------------------------------

void InitSurface::initIdealizedCells( MyTissue &T,
                                      const std::size_t founderCells,
                                      SurfaceBaseClass &surface )
{
  _radialModel = false;
  std::size_t lCounter = 1;
  switch( founderCells )
  {
    case 1:
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 1., 1. ),
                        lCounter, surface ); 
    break;
    case 2:
    for( std::size_t c=0; c < 2; c++, lCounter++ )
    {
      this->generateCell( T, std::make_pair( 0. + c*0.5, 0. ),
                          std::make_pair( 0.5, 1. ),
                          lCounter, surface );
    }
    
//     this->generateCell( T, std::make_pair( 0., 0. ),
//                           std::make_pair( 1./4., 1. ),
//                           lCounter, surface );
//     
//     lCounter++;
//     
//     this->generateCell( T, std::make_pair( 1./4., 0. ),
//                           std::make_pair( 3./4., 1. ),
//                           lCounter, surface );
    
    break;
    case 6:
    // 6 cells in total
    for( std::size_t t = 0; t < 2; t++ )
      for( std::size_t b = 0; b < 2; b++, lCounter++ )
      {
        double u = 1./3. + b*1./6.;
        double v = 0. + t*1./2.;
        double length = 1./6.;
          
        this->generateCell( T, std::make_pair( u, v ),
                            std::make_pair( length, 1./2. ),
                            lCounter, surface );
      }
      
    // outer cells
    for( std::size_t c = 0; c < 2; c++, lCounter++ )
    {
      double u = 0.;
      if( c == 1 )
        u = 2./3.;
      double v = 0.;
      double length = 1./3.;
        
      // set wall for which an additional junction has to be included
      std::size_t addJunctionToWall;
      if( c == 0 )
        addJunctionToWall = 2;
      else
        addJunctionToWall = 0;
      
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( length, 1. ),
                          lCounter, surface, addJunctionToWall );
    }
    
    /*
    // init the four inner smaller cells
    for( std::size_t w = 0; w < 2; w++ )
      for( std::size_t h = 0; h < 2; h++, lCounter++ )
      {
        double u = 1./4. + w*1./4.;
        double v = 0. + h*1./2.;
        this->generateCell( T, std::make_pair( u, v ),
                            std::make_pair( 1./4., 1./2. ),
                            lCounter, surface );
      }
    
    // outer left and right cell
    for( std::size_t h = 0; h < 2; h++, lCounter++ )
    {
      double u = 0. + h*3./4.;
      double v = 0.;
      // set wall for which an additional junction has to be included
      std::size_t addJunctionToWall;
      if( h == 0 )
        addJunctionToWall = 2;
      else
        addJunctionToWall = 0;
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( 1./4., 1. ),
                          lCounter, surface, addJunctionToWall );
    }
    */
    break;
    case 8:
    for( std::size_t h = 0; h < 2; h++, lCounter++ )
    {
      double u = 0.;
      double v = h*1./2.;
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( 1./3., 1./2. ),
                          lCounter, surface );
    }
    
    // init the inner smaller cells
    for( std::size_t w = 0; w < 2; w++ )
      for( std::size_t h = 0; h < 2; h++, lCounter++ )
      {
        double u = 1./3. + w*1./6.;
        double v = 0. + h*1./2.;
        this->generateCell( T, std::make_pair( u, v ),
                            std::make_pair( 1./6., 1./2. ),
                            lCounter, surface );
      }
    
    for( std::size_t h = 0; h < 2; h++, lCounter++ )
    {
      double u = 2./3.;
      double v = h*1./2.;
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( 1./3., 1./2. ),
                          lCounter, surface );
    }
    break;
    default:
    std::cerr << "Selected number of cells at the beginning is not implemented yet!" << std::endl;
    break;
  }
  
  forall(const cell& c, T.C)
    c->cellNumber = T.C.size();
}

// ---------------------------------------------------------------------

void InitSurface::initRadialCells( MyTissue &T,
                                   SurfaceBaseClass &surface )
{
  _radialModel = true;
  std::size_t lCounter = 1;
  int maxCells = 5;
  int cellFile = -2;
  double uPos = 0.;
  for( std::size_t c=1; c <= maxCells; c++, lCounter++, cellFile += 1 )
  {
    this->generateCell( T, std::make_pair( uPos, 0. ),
                        std::make_pair( 1./(double)maxCells, 1. ),
                        lCounter, surface, 4, cellFile );
    uPos += 1./(double)maxCells;
  }
  
  forall(const cell& c, T.C)
    c->cellNumber = T.C.size();
}

// ---------------------------------------------------------------------

void InitSurface::initRealDataCells( MyTissue &T, SurfaceBaseClass &surface )
{
  _radialModel = false;
  std::cout << "Lateral root constellation of data: " << _dataset << std::endl;
  std::size_t lCounter = 1;
  // 2 cells
  if( _dataset == "120830_raw" )
  {
    for( std::size_t c=0; c < 2; c++, lCounter++ )
    {
      this->generateCell( T, std::make_pair( 0. + c*0.5, 0. ),
                          std::make_pair( 0.5, 1. ),
                          lCounter, surface );
    }
  }
  // 2 cells
  else if( _dataset == "121204_raw_2014" )
  {
    for( std::size_t c = 0; c < 2; c++, lCounter++ )
    {
      double u = 0. + c*2./5.;
      double v = 0.;
      double length;
      if( c == 0 )
        length = 2./5.;
      else
        length = 3./5.;
        
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( length, 1. ),
                          lCounter, surface );
    }
  }
  // 5 cells
  else if( _dataset == "121211_raw" )
  {
    if( _forceInitialSituation )
    {
      for( std::size_t c=0; c < 2; c++, lCounter++ )
      {
        this->generateCell( T, std::make_pair( 0. + c*0.5, 0. ),
                            std::make_pair( 0.5, 1. ),
                            lCounter, surface );
      }
    }
    else
    {
      this->generateCell( T, std::make_pair( 0., 0. ),
                          std::make_pair( 1./3., 1. ),
                          lCounter, surface );
      
      lCounter++;
      
      // inner smaller cells
      for( std::size_t c = 0; c < 3; c++, lCounter++ )
      {
        double u = 1./3. + c*1./9.;
        double v = 0.;
        this->generateCell( T, std::make_pair( u, v ),
                            std::make_pair( 1./9., 1. ),
                            lCounter, surface );
      }
      
      this->generateCell( T, std::make_pair( 2./3., 0. ),
                          std::make_pair( 1./3., 1. ),
                          lCounter, surface );
    }
  }
  // 1 cell
  else if( _dataset == "130508_raw" )
  {
    if( _forceInitialSituation )
    {
      for( std::size_t c=0; c < 2; c++, lCounter++ )
      {
        this->generateCell( T, std::make_pair( 0. + c*0.5, 0. ),
                            std::make_pair( 0.5, 1. ),
                            lCounter, surface );
      }
    }
    else
    {
      this->generateCell( T, std::make_pair( 0., 0. ),
                          std::make_pair( 1., 1. ),
                          lCounter, surface );
    }
  }
  // 3 cells
  else if( _dataset == "130607_raw" )
  {
    if( _forceInitialSituation )
    {
      for( std::size_t c=0; c < 2; c++, lCounter++ )
      {
        this->generateCell( T, std::make_pair( 0. + c*0.5, 0. ),
                            std::make_pair( 0.5, 1. ),
                            lCounter, surface );
      }
    }
    else
    {
      for( std::size_t c = 0; c < 3; c++, lCounter++ )
      {
        double length,u,v;
        if( c == 0 )
        {
          length = 1./5.;
          u = 0.;
          v = 0.;
        }
        else
        {
          length = 2./5.;
          u = 1./5. + (c-1)*2./5.;
          v = 0.;
        }
          
        this->generateCell( T, std::make_pair( u, v ),
                            std::make_pair( length, 1. ),
                            lCounter, surface );
      }
    }
  }
  else
    std::cerr << "Selected data set is not supported!" << std::endl;
  
  forall(const cell& c, T.C)
    c->cellNumber = T.C.size();
}

// ---------------------------------------------------------------------

// generate a cell starting at the bottom left
// with different lengths
void InitSurface::generateCell( MyTissue &T,
                                const std::pair<double, double> &start,
                                const std::pair<double, double> &length,
                                const std::size_t treeId,
                                SurfaceBaseClass &surface,
                                const std::size_t addJunctionToWall,
                                const int cellFile )
{
  // first determine which kind of surface is used
  bool bezier = true;
  Surface *bez = dynamic_cast<Surface*>( &surface );
  if( !bez )
    bezier = false;
  
  // addJunctionToWall can be a number in [0,4] for which each number
  // represents a wall:
  // 0 -> left, 1 -> top, 2-> right, 3 -> down, 4 -> none
  // for one of these walls an additional junction is generated for which
  // an existing one is searched (if it was already generated before) such
  // that the old and new junctions are merged
  
  // set of junctions for the cell
  std::vector<junction> vs;
  
  // if _lod = 1 -> curLength = 0
  // order of generated points:
  // 1 --------- 2
  // |           |
  // |           |
  // 0 --------- 3
  
  for( std::size_t w = 0; w < 4 * _lod; w++ )
  {
    unsigned int uiw = (unsigned int)(w/_lod);
    double cellLength;
    double curLength = (double)w/(double)_lod - (double)uiw;
    double u,v;
    double mu,mv;
    switch(uiw)
    {
      case 0:
        u = start.first;
        v = start.second + curLength*length.second;
        if( addJunctionToWall == 0 )
        {
          mu = start.first;
          mv = start.second + length.second/2. + curLength*length.second;
        }
        break;
      case 1:
        u = start.first + curLength*length.first;
        v = start.second + length.second;
        if( addJunctionToWall == 1 )
        {
          mu = start.first + length.first/2. + curLength*length.first;
          mv = start.second + length.second;
        }
        break;
      case 2:
        u = start.first + length.first;
        v = start.second + length.second - curLength*length.second;
        if( addJunctionToWall == 2 )
        {
          mu = start.first + length.first;
          mv = start.second + length.second/2. - curLength*length.second;
        }
        break;
      case 3:
        u = start.first + length.first - curLength*length.first;
        v = start.second;
        if( addJunctionToWall == 3 )
        {
          mu = start.first + length.first/2. - curLength*length.first;
          mv = start.second;
        }
        break;
      default: u = v = 0.; break;
    }
    
    junction j;
    j->id = _jIDCounter;
    if( bezier )
    {
      j->sp.setUV( u, v );
      surface.initPos( j->sp );
    }
    else
    {
      Point3d curPos = this->determinePos( u, v );
      surface.setPos( j->sp, curPos );
    }
    ModelUtils::findNearestPointToMerge( T, j );
    vs.push_back( j );
    _jIDCounter++;
    
    if( uiw == addJunctionToWall )
    {
      junction j;
      j->id = _jIDCounter;
      if( bezier )
      {
        j->sp.setUV( mu, mv );
        surface.initPos( j->sp );
      }
      else
      {
        Point3d cPos = this->determinePos( mu, mv );
        surface.setPos( j->sp, cPos );
      }
      ModelUtils::findNearestPointToMerge( T, j );
      vs.push_back( j );
      _jIDCounter++;
    }
  }
  
  // TODO: only one cell for real data
//   for( std::size_t w = 0; w < conPoints.size(); w++ )
//   {
//     Point3d curPos = conPoints.at(w);
//     junction j;
//     j->id = _jIDCounter;
//     lateralRoot.setPos( j->sp, curPos );
//     ModelUtils::findNearestPointToMerge( T, j );
//     vs.push_back(j);
//     _jIDCounter++;
//   }
  
  cell c;
  c->treeId = treeId;
  c->id = _IDCounter;
  c->parentId = _IDCounter;
  c->timeStep = _time;
  c->previousAngle = 0.;
  c->angle = 0.;
  c->previousDivDir = Point3d( 0., 0., 0. );
  c->divDir = Point3d( 0., 0., 0. );
  c->divType = DivisionType::NONE;
  c->layerValue = 1;
  c->divisionSequence = "0";
  c->divisionLetterSequence = "";
  c->radialDivisionSequence = "";
  c->cellCycle = 0;
  c->cellFile = cellFile;
  c->cellFileColoringIndex = 1 + std::abs(c->cellFile);
  c->cellFileSequence = std::to_string( cellFile );
  c->periCycle = 0;
  T.addCell( c, vs );
  
  std::vector<Point3d> polygon;
  double xMin = 5000000.;
  double xMax = -5000000.;
  Point3d center;
  forall(const junction& j, T.S.neighbors(c))
  {
    polygon.push_back( j->getPos() );
    center += j->getPos();
    
    if( j->getPos().i() < xMin )
      xMin = j->getPos().i();
    
    if( j->getPos().i() > xMax )
      xMax = j->getPos().i();
  }
  center /= polygon.size();
  
  c->xMin = xMin;
  c->xMax = xMax;
  c->divType = DivisionType::NONE;
  c->centerPos.push_back( center );
  c->initialArea = geometry::polygonArea(polygon);
  c->area = c->initialArea;
  c->initialLongestWallLength = ModelUtils::determineLongestWallLength( c, T );
  c->longestWallLength = c->initialLongestWallLength;
  
  // set the division properties for the first division
  // anticlinal division in 1/3:2/3 ratio resulting in four cells
  // **********************************
  // *         *     *     *          *
  // *         *     *     *          *
  // *         *     *     *          *
  // **********************************
  if( _forceInitialSituation && c->id <= 2 && !_radialModel )
  {
    // determine xLength of cell
    double xLength = std::fabs( c->xMax - c->xMin );
    
    if( c->id == 1 )
      center.i() = c->xMin + 2.*xLength/3.; 
    else if( c->id == 2 )
      center.i() = c->xMin + 1.*xLength/3.; 
  }
  
  c->center = center;
  surface.setPos( c->sp, center );
  ModelExporter::exportLineageInformation( _lineageFileName, c, T, false );
  
  // afterwards increment the id counter
  _IDCounter++;
}

// ---------------------------------------------------------------------

Point3d InitSurface::determinePos( const double u, const double v )
{
  std::string name = "conPoints-";
  name += _dataset;
  if( _useAutomaticContourPoints )
    name += "_auto.txt";
  else
    name += ".txt";
  
  std::vector<Point3d> conPoints =
  ModelUtils::loadContourPoints( name, _surfaceScale );
  
  // u and v coordinates of start should correspond to the
  // following contour points
  // u = 0, v = 0 -> conPoints.at(14)
  // u = 0, v = 1 -> conPoints.at(0)
  // u = 1, v = 0 -> conPoints.at(8)
  // u = 1, v = 1 -> conPoints.at(6)
  
  Point3d pos;
  double h = v;
  double xPos = 6. * u;
  std::size_t xI = (std::size_t)xPos;
  double factor = xPos - (double)xI;

  Point3d pos1;
  Point3d pos2;
  if( xI < 6 )
  {
    // left to right: between conPoints at 0 and 6
    pos1 = (1.-factor) * conPoints.at(xI) + factor * conPoints.at(xI+1);
    // right to left: between conPoints at 14 and 8
    pos2 = (1.-factor) * conPoints.at(14-xI) + factor * conPoints.at(14-xI-1);
  }
  else
  {
    // left to right: between conPoints at 0 and 6
    pos1 = conPoints.at(xI);
    // right to left: between conPoints at 14 and 8
    pos2 = conPoints.at(14-xI);
  }
  
  pos = h * pos1 + (1.-h) * pos2;
  return pos;
}

// ---------------------------------------------------------------------
