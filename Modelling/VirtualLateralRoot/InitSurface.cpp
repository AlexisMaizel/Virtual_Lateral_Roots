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
                         const std::size_t timeSteps,
                         const bool forceInitialSituation )
{
  if( _lod == 0 )
    _lod = 1;
  else
    _lod = lod;
  
  _IDCounter = 1;
  _jIDCounter = 0;
  _time = 1;
  _maxTime = timeSteps;
  _lineageFileName = lineageFileName;
  _forceInitialSituation = forceInitialSituation;
}

// ---------------------------------------------------------------------

void InitSurface::initModelBasedOnBezier( MyTissue &T,
                                           const std::size_t cellNumber,
                                           Surface &lateralRoot )
{
  std::size_t lCounter = 1;
  switch( cellNumber )
  {
    case 1:
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 1., 1. ),
                        lCounter, lateralRoot ); 
    break;
    case 2:
    for( std::size_t c=0; c < 2; c++, lCounter++ )
    {
      this->generateCell( T, std::make_pair( 0. + c*0.5, 0. ),
                          std::make_pair( 0.5, 1. ),
                          lCounter, lateralRoot );
    }
    
//     this->generateCell( T, std::make_pair( 0., 0. ),
//                           std::make_pair( 1./4., 1. ),
//                           lCounter, lateralRoot );
//     
//     lCounter++;
//     
//     this->generateCell( T, std::make_pair( 1./4., 0. ),
//                           std::make_pair( 3./4., 1. ),
//                           lCounter, lateralRoot );
    
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
                            lCounter, lateralRoot );
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
                          lCounter, lateralRoot, addJunctionToWall );
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
                            lCounter, lateralRoot );
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
                          lCounter, lateralRoot, addJunctionToWall );
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
                          lCounter, lateralRoot );
    }
    
    // init the inner smaller cells
    for( std::size_t w = 0; w < 2; w++ )
      for( std::size_t h = 0; h < 2; h++, lCounter++ )
      {
        double u = 1./3. + w*1./6.;
        double v = 0. + h*1./2.;
        this->generateCell( T, std::make_pair( u, v ),
                            std::make_pair( 1./6., 1./2. ),
                            lCounter, lateralRoot );
      }
    
    for( std::size_t h = 0; h < 2; h++, lCounter++ )
    {
      double u = 2./3.;
      double v = h*1./2.;
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( 1./3., 1./2. ),
                          lCounter, lateralRoot );
    }
    break;
    default:
    std::cerr << "Selected number of cells at the beginning is not implemented yet!" << std::endl;
    break;
  }
}

// ---------------------------------------------------------------------

void InitSurface::initLateralRootBasedOnBezier( MyTissue &T,
                                                 const std::string &dataset,
                                                 Surface &lateralRoot )
{
  std::cout << "Lateral root constellation of data: " << dataset << std::endl;
  std::size_t lCounter = 1;
  // 2 cells
  if( dataset == "120830_raw" )
  {
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 0.55, 1. ),
                        lCounter, lateralRoot );
    
    lCounter++;
    
    this->generateCell( T, std::make_pair( 0.55, 0. ),
                        std::make_pair( 0.45, 1. ),
                        lCounter, lateralRoot );
  }
  // 2 cells
  else if( dataset == "121204_raw_2014" )
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
                          lCounter, lateralRoot );
    }
  }
  // 5 cells
  else if( dataset == "121211_raw" )
  {
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 1./3., 1. ),
                        lCounter, lateralRoot );
    
    lCounter++;
    
    // inner smaller cells
    for( std::size_t c = 0; c < 3; c++, lCounter++ )
    {
      double u = 1./3. + c*1./9.;
      double v = 0.;
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( 1./9., 1. ),
                          lCounter, lateralRoot );
    }
    
    this->generateCell( T, std::make_pair( 2./3., 0. ),
                        std::make_pair( 1./3., 1. ),
                        lCounter, lateralRoot );
  }
  // 1 cell
  else if( dataset == "130508_raw" )
  {
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 1., 1. ),
                        lCounter, lateralRoot );
  }
  // 3 cells
  else if( dataset == "130607_raw" )
  {
    for( std::size_t c = 0; c < 3; c++, lCounter++ )
    {
      double length,u,v;
      switch( c )
      {
        case 0:
        length = 1./12.;
        u = 0.;
        v = 0.;
        break;
        case 1:
        length = 5./12.;
        u = 1./12.;
        v = 0.;
        break;
        case 2:
        length = 6./12.;
        u = 6./12.;
        v = 0.;
        break;
      }
        
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( length, 1. ),
                          lCounter, lateralRoot );
    }
  }
  else
    std::cerr << "Selected data set is not supported!" << std::endl;  
}

// ---------------------------------------------------------------------

void InitSurface::initLateralRootBasedOnRealData( MyTissue &T,
                                                   RealSurface &lateralRoot,
                                                   const std::string &dataset,
                                                   const double surfaceScale,
                                                   const bool useAutomaticContourPoints )
{
  std::cout << "Lateral root constellation of data: " << dataset << std::endl;
  
  std::string fileName = dataset;
  if( dataset.compare( 0, 7, "Average") == 0 )
    fileName = "Average";
  
  std::string name = "conPoints-";
  name += fileName;
  if( useAutomaticContourPoints )
    name += "_auto.txt";
  else
    name += ".txt";
  
  std::vector<Point3d> conPoints = ModelUtils::loadContourPoints( name, surfaceScale );
  std::size_t lCounter = 1;
    
  // 2 cells
  if( dataset == "120830_raw" )
  {
    for( std::size_t c=0; c < 2; c++, lCounter++ )
    {
      this->generateCell( T, std::make_pair( 0. + c*0.5, 0. ),
                          std::make_pair( 0.5, 1. ),
                          lCounter, lateralRoot, conPoints );
    }
  }
  // 2 cells
  else if( dataset == "121204_raw_2014" )
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
                          lCounter, lateralRoot, conPoints );
    }
  }
  // 5 cells
  else if( dataset == "121211_raw" )
  {
    if( _forceInitialSituation )
    {
      for( std::size_t c=0; c < 2; c++, lCounter++ )
      {
        this->generateCell( T, std::make_pair( 0. + c*0.5, 0. ),
                            std::make_pair( 0.5, 1. ),
                            lCounter, lateralRoot, conPoints );
      }
    }
    else
    {
      this->generateCell( T, std::make_pair( 0., 0. ),
                          std::make_pair( 1./3., 1. ),
                          lCounter, lateralRoot, conPoints );
      
      lCounter++;
      
      // inner smaller cells
      for( std::size_t c = 0; c < 3; c++, lCounter++ )
      {
        double u = 1./3. + c*1./9.;
        double v = 0.;
        this->generateCell( T, std::make_pair( u, v ),
                            std::make_pair( 1./9., 1. ),
                            lCounter, lateralRoot, conPoints );
      }
      
      this->generateCell( T, std::make_pair( 2./3., 0. ),
                          std::make_pair( 1./3., 1. ),
                          lCounter, lateralRoot, conPoints );
    }
  }
  // 1 cell
  else if( dataset == "130508_raw" )
  {
    if( _forceInitialSituation )
    {
      for( std::size_t c=0; c < 2; c++, lCounter++ )
      {
        this->generateCell( T, std::make_pair( 0. + c*0.5, 0. ),
                            std::make_pair( 0.5, 1. ),
                            lCounter, lateralRoot, conPoints );
      }
    }
    else
    {
      this->generateCell( T, std::make_pair( 0., 0. ),
                          std::make_pair( 1., 1. ),
                          lCounter, lateralRoot, conPoints, 4, true );
    }
  }
  // 3 cells
  else if( dataset == "130607_raw" )
  {
    if( _forceInitialSituation )
    {
      for( std::size_t c=0; c < 2; c++, lCounter++ )
      {
        this->generateCell( T, std::make_pair( 0. + c*0.5, 0. ),
                            std::make_pair( 0.5, 1. ),
                            lCounter, lateralRoot, conPoints );
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
                            lCounter, lateralRoot, conPoints );
      }
    }
  }
  else if( dataset == "Average2" )
  {
    // 2 cells
    for( std::size_t c = 0; c < 2; c++, lCounter++ )
    {
      double length = 1./2.;
      this->generateCell( T, std::make_pair( c*1./2., 0. ),
                             std::make_pair( length, 1. ),
                             lCounter, lateralRoot, conPoints );
    }
  }
  else if( dataset == "Average4" )
  {
    // 4 cells
    for( std::size_t c = 0; c < 4; c++, lCounter++ )
    {
      double u = 0. + c*1./4.;
      double v = 0.;
      double length = 1./4.;
        
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( length, 1. ),
                          lCounter, lateralRoot, conPoints );
    }
  }
  else if( dataset == "Average6" )
  {
    // 6 cells in total
    for( std::size_t t = 0; t < 2; t++ )
      for( std::size_t b = 0; b < 2; b++, lCounter++ )
      {
        double u = 1./3. + b*1./6.;
        double v = 0. + t*1./2.;
        double length = 1./6.;
          
        this->generateCell( T, std::make_pair( u, v ),
                            std::make_pair( length, 1./2. ),
                            lCounter, lateralRoot, conPoints );
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
      // right wall
      if( c == 0 )
        addJunctionToWall = 2;
      // left wall
      else if( c == 1 )
        addJunctionToWall = 0;
      // none
      else
        addJunctionToWall = 4;
      
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( length, 1. ),
                          lCounter, lateralRoot,
                          conPoints, addJunctionToWall );
    }
  }
  else if( dataset == "Average8" )
  {
    // 8 cells in total
    for( std::size_t t = 0; t < 2; t++ )
      for( std::size_t b = 0; b < 2; b++, lCounter++ )
      {
        double u = 1./4. + b*1./4.;
        double v = 0. + t*1./2.;
        double length = 1./4.;
          
        this->generateCell( T, std::make_pair( u, v ),
                            std::make_pair( length, 1./2. ),
                            lCounter, lateralRoot, conPoints );
      }
      
    for( std::size_t c = 0; c < 4; c++, lCounter++ )
    {
      double u = 0. + c*1./8.;
      if( c > 1 )
        u = 1./2. + c*1./8.;
      double v = 0.;
      double length = 1./8.;
        
      // set wall for which an additional junction has to be included
      std::size_t addJunctionToWall = 4;
      // right wall
      if( c == 1 )
        addJunctionToWall = 2;
      // left wall
      else if( c == 2 )
        addJunctionToWall = 0;
      // none
      else
        addJunctionToWall = 4;
      
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( length, 1. ),
                          lCounter, lateralRoot,
                          conPoints, addJunctionToWall );
    }
  }
  else
    std::cerr << "Selected data set is not supported!" << std::endl;
}

// ---------------------------------------------------------------------

// generate a cell starting at the bottom left
// with different lengths
void InitSurface::generateCell( MyTissue &T,
                                 const std::pair<double, double> &start,
                                 const std::pair<double, double> &length,
                                 const std::size_t treeId,
                                 Surface &lateralRoot,
                                 const std::size_t addJunctionToWall )
{
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
    lateralRoot.InitPoint( j->sp, u, v );
    ModelUtils::findNearestPointToMerge( T, j );
    vs.push_back(j);
    _jIDCounter++;
    
    if( uiw == addJunctionToWall )
    {
      junction j;
      j->id = _jIDCounter;
      lateralRoot.InitPoint( j->sp, mu, mv );
      ModelUtils::findNearestPointToMerge( T, j );
      vs.push_back(j);
      _jIDCounter++;
    }
  }
  
  //std::cout << "cell" << std::endl;
  
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
  c->cellCycle = 0;
  c->periCycle = 0;
  T.addCell( c, vs );
  
  std::vector<Point3d> polygon;
  double xMin = 5000000.;
  double xMax = -5000000.;
  Point3d center;
  forall(const junction& j, T.S.neighbors(c))
  {
    polygon.push_back( j->sp.Pos() );
    center += j->sp.Pos();
    
    if( j->sp.Pos().i() < xMin )
      xMin = j->sp.Pos().i();
    
    if( j->sp.Pos().i() > xMax )
      xMax = j->sp.Pos().i();
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
  if( _forceInitialSituation && c->id <= 2 )
  {
    // determine xLength of cell
    double xLength = std::fabs( c->xMax - c->xMin );
    
    if( c->id == 1 )
      center.i() = c->xMin + 2.*xLength/3.; 
    else if( c->id == 2 )
      center.i() = c->xMin + 1.*xLength/3.; 
  }
  
  c->center = center;
  
  lateralRoot.SetPoint(c->sp, c->sp, center);

  ModelExporter::exportLineageInformation( _lineageFileName, c, T, false );
  
  // afterwards increment the id counter
  _IDCounter++;
}

// ---------------------------------------------------------------------

void InitSurface::generateCell( MyTissue &T,
                                 const std::pair<double, double> &start,
                                 const std::pair<double, double> &length,
                                 const std::size_t treeId,
                                 RealSurface &lateralRoot,
                                 const std::vector<Point3d> &conPoints,
                                 const std::size_t addJunctionToWall,
                                 const bool oneCell )
{
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
  
  if( !oneCell )
  {
    for( std::size_t w = 0; w < 4*_lod; w++ )
    {
      unsigned int uiw = (unsigned int)(w/_lod);
      double cellLength;
      double curLength = (double)w/(double)_lod - (double)uiw;
      std::pair<double, double> iPos;
      std::pair<double, double> mPos;
      
      switch(uiw)
      {
        case 0:
          iPos.first = start.first;
          iPos.second = start.second + curLength*length.second;
          if( addJunctionToWall == 0 )
          {
            mPos.first = start.first;
            mPos.second = start.second + length.second/2. + curLength*length.second;
          }
          break;
        case 1:
          iPos.first = start.first + curLength*length.first;
          iPos.second = start.second + length.second;
          if( addJunctionToWall == 1 )
          {
            mPos.first = start.first + length.first/2. + curLength*length.first;
            mPos.second = start.second + length.second;
          }
          break;
        case 2:
          iPos.first = start.first + length.first;
          iPos.second = start.second + length.second - curLength*length.second;
          if( addJunctionToWall == 2 )
          {
            mPos.first = start.first + length.first;
            mPos.second = start.second + length.second/2. - curLength*length.second;
          }
          break;
        case 3:
          iPos.first = start.first + length.first - curLength*length.first;
          iPos.second = start.second;
          if( addJunctionToWall == 3 )
          {
            mPos.first = start.first + length.first/2. - curLength*length.first;
            mPos.second = start.second;
          }
          break;
        default: iPos.first = iPos.second = 0.; break;
      }
      
      Point3d curPos = this->determinePos( iPos, conPoints );
      junction j;
      j->id = _jIDCounter;
      lateralRoot.setPos( j->tp, curPos );
      ModelUtils::findNearestPointToMerge( T, j );
      vs.push_back(j);
      _jIDCounter++;
      
      if( uiw == addJunctionToWall )
      {
        Point3d cPos = this->determinePos( mPos, conPoints );
        junction j;
        j->id = _jIDCounter;
        lateralRoot.setPos( j->tp, cPos );
        ModelUtils::findNearestPointToMerge( T, j );
        vs.push_back(j);
        _jIDCounter++;
      }
    }
  }
  else
  {
    for( std::size_t w = 0; w < conPoints.size(); w++ )
    {
      Point3d curPos = conPoints.at(w);
      junction j;
      j->id = _jIDCounter;
      lateralRoot.setPos( j->tp, curPos );
      ModelUtils::findNearestPointToMerge( T, j );
      vs.push_back(j);
      _jIDCounter++;
    }
  }
  
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
  c->cellCycle = 0;
  c->periCycle = 0;
  T.addCell( c, vs );
  
  std::vector<Point2d> polygon;
  double xMin = 5000000.;
  double xMax = -5000000.;
  Point3d center;
  forall(const junction& j, T.S.neighbors(c))
  {
    polygon.push_back( j->tp.Pos() );
    center += Point3d( j->tp.Pos().i(), j->tp.Pos().j(), 0. );
    
    if( j->tp.Pos().i() < xMin )
      xMin = j->tp.Pos().i();
    
    if( j->tp.Pos().i() > xMax )
      xMax = j->tp.Pos().i();
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
  if( _forceInitialSituation && c->id <= 2 )
  {
    // determine xLength of cell
    double xLength = std::fabs( xMax - xMin );
    
    if( c->id == 1 )
      center.i() = xMin + 2.*xLength/3.; 
    else if( c->id == 2 )
      center.i() = xMin + 1.*xLength/3.; 
  }
  
  c->center = center;
  
  lateralRoot.setPos(c->tp, center);
  
  ModelExporter::exportLineageInformation( _lineageFileName, c, T, false );
  
  // afterwards increment the id counter
  _IDCounter++;  
}

// ---------------------------------------------------------------------

Point3d InitSurface::determinePos( const std::pair<double, double> &coord,
                                    const std::vector<Point3d> &conPoints )
{
  // u and v coordinates of start should correspond to the
  // following contour points
  // u = 0, v = 0 -> conPoints.at(14)
  // u = 0, v = 1 -> conPoints.at(0)
  // u = 1, v = 0 -> conPoints.at(8)
  // u = 1, v = 1 -> conPoints.at(6)
  
  Point3d pos;
  double h = coord.second;
  double xPos = 6. * coord.first;
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
