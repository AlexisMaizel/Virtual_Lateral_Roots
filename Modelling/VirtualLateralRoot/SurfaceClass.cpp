#include "SurfaceClass.h"

#include "ModelExporter.h"

/**
  @file   SurfaceClass.cpp
  @brief  Contains class for setting initial surfaces of model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

// ---------------------------------------------------------------------

void SurfaceClass::init( const std::size_t lod,
                         const std::string &lineageFileName,
                         const bool exportLineage,
                         const std::size_t timeSteps )
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
  _exportLineage = exportLineage;
}

// ---------------------------------------------------------------------

void SurfaceClass::initModelBasedOnBezier( MyTissue &T,
                                           const std::size_t cellNumber,
                                           Surface &lateralRoot )
{
  // render eight cells at the beginning for which each pair shares an
  // area ratio of 2:1. The bigger cell will always be located at the left
  // and right boundary
  idPairSet sharedJunctions;
  std::size_t lineageCounter = 1;
  
  switch( cellNumber )
  {
    case 8:
    // insert ids of shared junctions
    for( std::size_t c = 1; c < 7; c++ )
    {
      sharedJunctions.insert( std::make_pair( 4*(c+1), 4*(c+1)-5 ) );
      sharedJunctions.insert( std::make_pair( 4*(c+1)+1, 4*(c+1)-6 ) );
      if( c < 5 )
      {
        sharedJunctions.insert( std::make_pair( 8*c-4, 8*c-7 ) );
        sharedJunctions.insert( std::make_pair( 8*c-1, 8*c-6 ) );
      }
      
      // vertices shared by 4 points
      if( c < 4 )
      {
        sharedJunctions.insert( std::make_pair( 8*c+1, 8*c-1 ) );
        sharedJunctions.insert( std::make_pair( 8*(c+1)-4, 8*(c-1)+2 ) );
      }
    }
    
    for( std::size_t h = 0; h < 2; h++ )
    {
      double u = 0.;
      double v = h*1./2.;
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( 1./3., 1./2. ),
                          lineageCounter, sharedJunctions,
                          lateralRoot );
      
      lineageCounter++;
    }
    
    // init the inner smaller cells
    for( std::size_t w = 0; w < 2; w++ )
    {
      for( std::size_t h = 0; h < 2; h++ )
      {
        double u = 1./3. + w*1./6.;
        double v = 0. + h*1./2.;
        this->generateCell( T, std::make_pair( u, v ),
                            std::make_pair( 1./6., 1./2. ),
                            lineageCounter, sharedJunctions,
                            lateralRoot );
        
        lineageCounter++;
      }
    }
    
    for( std::size_t h = 0; h < 2; h++ )
    {
      double u = 2./3.;
      double v = h*1./2.;
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( 1./3., 1./2. ),
                          lineageCounter, sharedJunctions,
                          lateralRoot );
      
      lineageCounter++;
    }
    break;
    case 1:
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 1., 1. ),
                        2, sharedJunctions, lateralRoot ); 
    break;
    case 2:
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 0.5, 1. ),
                        1, sharedJunctions, lateralRoot );
      
    // insert ids of shared junctions
    sharedJunctions.insert( std::make_pair( 5, 2 ) );
    sharedJunctions.insert( std::make_pair( 4, 3 ) );
    
    this->generateCell( T, std::make_pair( 0.5, 0. ),
                        std::make_pair( 0.5, 1. ),
                        2, sharedJunctions, lateralRoot );
    break;
    default:
    std::cerr << "Selected number of cells at the beginning is not implemented yet!" << std::endl;
    break;
  }
}

// ---------------------------------------------------------------------

void SurfaceClass::initLateralRootBasedOnBezier( MyTissue &T,
                                                 const std::string &dataset,
                                                 Surface &lateralRoot )
{
  std::cout << "Lateral root constellation of data: " << dataset << std::endl;
  idPairSet sharedJunctions;
  
  if( dataset == "120830_raw" )
  {
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 0.55, 1. ),
                        1, sharedJunctions, lateralRoot );
    
    // insert ids of shared junctions
    sharedJunctions.insert( std::make_pair( 5*_lod, 2*_lod ) );
    sharedJunctions.insert( std::make_pair( 4*_lod, 3*_lod ) );
    
    for( std::size_t l=1; l<_lod;l++ )
    {
      std::size_t u = 4*_lod + l;
      std::size_t v = 3*_lod - l;
      sharedJunctions.insert( std::make_pair( u, v ) );
    }
    
    this->generateCell( T, std::make_pair( 0.55, 0. ),
                        std::make_pair( 0.45, 1. ),
                        2, sharedJunctions, lateralRoot );
  }
  else if( dataset == "121204_raw_2014" )
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
        
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( length, 1. ),
                          lineageCounter, sharedJunctions, lateralRoot );
      
      lineageCounter++;
      
      // insert ids of shared junctions
      sharedJunctions.insert( std::make_pair( 5*_lod, 2*_lod ) );
      sharedJunctions.insert( std::make_pair( 4*_lod, 3*_lod ) );
    
      for( std::size_t l=1; l<_lod;l++ )
      {
        std::size_t u = 4*_lod + l;
        std::size_t v = 3*_lod - l;
        sharedJunctions.insert( std::make_pair( u, v ) );
      }
    }
  }
  else if( dataset == "121211_raw" )
  {
    std::size_t lineageCounter = 1;
  
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 1./3., 1. ),
                        lineageCounter, sharedJunctions, lateralRoot );
    
    // insert ids of shared junctions
    sharedJunctions.insert( std::make_pair( 5*_lod, 2*_lod ) );
    sharedJunctions.insert( std::make_pair( 4*_lod, 3*_lod ) );
  
    for( std::size_t l=1; l<_lod;l++ )
    {
      std::size_t u = 4*_lod + l;
      std::size_t v = 3*_lod - l;
      sharedJunctions.insert( std::make_pair( u, v ) );
    }
    
    lineageCounter++;
    
    // inner smaller cells
    for( std::size_t c = 0; c < 3; c++ )
    {
      double u = 1./3. + c*1./9.;
      double v = 0.;
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( 1./9., 1. ),
                          lineageCounter, sharedJunctions, lateralRoot );
      
      lineageCounter++;
      
      // insert ids of shared junctions
      sharedJunctions.insert( std::make_pair( (4*(c+2)+1)*_lod, (4*(c+2)-2)*_lod ) );
      sharedJunctions.insert( std::make_pair( 4*(c+2)*_lod, (4*(c+2)-1)*_lod ) );
    
      for( std::size_t l=1; l<_lod;l++ )
      {
        std::size_t u = 4*(c+2)*_lod + l;
        std::size_t v = (4*(c+2)-1)*_lod - l;
        sharedJunctions.insert( std::make_pair( u, v ) );
      }
    }
    
    this->generateCell( T, std::make_pair( 2./3., 0. ),
                        std::make_pair( 1./3., 1. ),
                        lineageCounter, sharedJunctions, lateralRoot );
  }
  else if( dataset == "130508_raw" )
  {
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 1., 1. ),
                        1, sharedJunctions, lateralRoot );
  }
  else if( dataset == "130607_raw" )
  {
    std::size_t lineageCounter = 1;
  
    for( std::size_t c = 0; c < 3; c++ )
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
                          lineageCounter, sharedJunctions, lateralRoot );
      
      lineageCounter++;
      
      // insert ids of shared junctions
      sharedJunctions.insert( std::make_pair( (4*(c+1)+1)*_lod, (4*(c+1)-2)*_lod ) );
      sharedJunctions.insert( std::make_pair( 4*(c+1)*_lod, (4*(c+1)-1)*_lod ) );
    
      for( std::size_t l=1; l<_lod;l++ )
      {
        std::size_t u = 4*(c+1)*_lod + l;
        std::size_t v = (4*(c+1)-1)*_lod - l;
        sharedJunctions.insert( std::make_pair( u, v ) );
      }
    }
  }
  else
    std::cerr << "Selected data set is not supported!" << std::endl;  
}

// ---------------------------------------------------------------------

void SurfaceClass::initLateralRootBasedOnRealData( MyTissue &T,
                                                   RealSurface &lateralRoot,
                                                   const std::string &dataset,
                                                   const double surfaceScale,
                                                   const bool useAutomaticContourPoints )
{
  std::cout << "Lateral root constellation of data: " << dataset << std::endl;
  idPairSet sharedJunctions;
  
  std::string name = "conPoints-";
  name += dataset;
  if( useAutomaticContourPoints )
    name += "_auto.txt";
  else
    name += ".txt";
  std::vector<Point3d> conPoints =
    ModelUtils::loadContourPoints( name, surfaceScale );
  
   if( dataset == "120830_raw" )
  {
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 0.5, 1. ),
                        1, sharedJunctions, lateralRoot, conPoints );
    
    // insert ids of shared junctions
    sharedJunctions.insert( std::make_pair( 5*_lod, 2*_lod ) );
    sharedJunctions.insert( std::make_pair( 4*_lod, 3*_lod ) );
    
    for( std::size_t l=1; l<_lod;l++ )
    {
      std::size_t u = 4*_lod + l;
      std::size_t v = 3*_lod - l;
      sharedJunctions.insert( std::make_pair( u, v ) );
    }
    
    this->generateCell( T, std::make_pair( 0.5, 0. ),
                        std::make_pair( 0.5, 1. ),
                        2, sharedJunctions, lateralRoot, conPoints );
  }
  else if( dataset == "121204_raw_2014" )
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
        
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( length, 1. ),
                          lineageCounter, sharedJunctions,
                          lateralRoot, conPoints );
      
      lineageCounter++;
      
      // insert ids of shared junctions
      sharedJunctions.insert( std::make_pair( 5*_lod, 2*_lod ) );
      sharedJunctions.insert( std::make_pair( 4*_lod, 3*_lod ) );
    
      for( std::size_t l=1; l<_lod;l++ )
      {
        std::size_t u = 4*_lod + l;
        std::size_t v = 3*_lod - l;
        sharedJunctions.insert( std::make_pair( u, v ) );
      }
    }
  }
  else if( dataset == "121211_raw" )
  {
    std::size_t lineageCounter = 1;
  
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 1./3., 1. ),
                        lineageCounter, sharedJunctions,
                        lateralRoot, conPoints );
    
    // insert ids of shared junctions
    sharedJunctions.insert( std::make_pair( 5*_lod, 2*_lod ) );
    sharedJunctions.insert( std::make_pair( 4*_lod, 3*_lod ) );
  
    for( std::size_t l=1; l<_lod;l++ )
    {
      std::size_t u = 4*_lod + l;
      std::size_t v = 3*_lod - l;
      sharedJunctions.insert( std::make_pair( u, v ) );
    }
    
    lineageCounter++;
    
    // inner smaller cells
    for( std::size_t c = 0; c < 3; c++ )
    {
      double u = 1./3. + c*1./9.;
      double v = 0.;
      this->generateCell( T, std::make_pair( u, v ),
                          std::make_pair( 1./9., 1. ),
                          lineageCounter, sharedJunctions,
                          lateralRoot, conPoints );
      
      lineageCounter++;
      
      // insert ids of shared junctions
      sharedJunctions.insert( std::make_pair( (4*(c+2)+1)*_lod, (4*(c+2)-2)*_lod ) );
      sharedJunctions.insert( std::make_pair( 4*(c+2)*_lod, (4*(c+2)-1)*_lod ) );
    
      for( std::size_t l=1; l<_lod;l++ )
      {
        std::size_t u = 4*(c+2)*_lod + l;
        std::size_t v = (4*(c+2)-1)*_lod - l;
        sharedJunctions.insert( std::make_pair( u, v ) );
      }
    }
    
    this->generateCell( T, std::make_pair( 2./3., 0. ),
                        std::make_pair( 1./3., 1. ),
                        lineageCounter, sharedJunctions,
                        lateralRoot, conPoints );
  }
  else if( dataset == "130508_raw" )
  {
    this->generateCell( T, std::make_pair( 0., 0. ),
                        std::make_pair( 1., 1. ),
                        1, sharedJunctions, lateralRoot,
                        conPoints, true );
  }
  else if( dataset == "130607_raw" )
  {
    std::size_t lineageCounter = 1;
  
    for( std::size_t c = 0; c < 3; c++ )
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
                          lineageCounter, sharedJunctions,
                          lateralRoot, conPoints );
      
      lineageCounter++;
      
      // insert ids of shared junctions
      sharedJunctions.insert( std::make_pair( (4*(c+1)+1)*_lod, (4*(c+1)-2)*_lod ) );
      sharedJunctions.insert( std::make_pair( 4*(c+1)*_lod, (4*(c+1)-1)*_lod ) );
    
      for( std::size_t l=1; l<_lod;l++ )
      {
        std::size_t u = 4*(c+1)*_lod + l;
        std::size_t v = (4*(c+1)-1)*_lod - l;
        sharedJunctions.insert( std::make_pair( u, v ) );
      }
    }
  }
  else
    std::cerr << "Selected data set is not supported!" << std::endl; 
  
  /*
  // set of junctions for the cell
  std::vector<junction> vs;
  
  //Point3d dataMean( 289.023540405678, -25.7548027981398, 0. );
      
  //std::vector<Point3d> conPoints;
  //conPoints.push_back( cPoints.at(0) );
  //conPoints.push_back( cPoints.at(6) );
  //conPoints.push_back( cPoints.at(8) );
  //conPoints.push_back( cPoints.at(14) );
  
  for( std::size_t w = 0; w < conPoints.size(); w++ )
  {
    //conPoints.at(w) -= dataMean;
    Point3d curPos = conPoints.at(w);
    Point3d nextPos;
    if( w != conPoints.size()-1 )
      nextPos = conPoints.at(w+1);
    else
      nextPos = conPoints.at(0);
    
    for( std::size_t d = 0; d < _lod; ++d )
    {
      double factor = (double)d/(double)_lod;
      Point3d pos = curPos + factor*(nextPos-curPos);
      junction j;
      j->id = _jIDCounter;
      lateralRoot.setPos( j->tp, pos );
      vs.push_back(j);
      _jIDCounter++;
    }
  }
  
  cell c;
  c->treeId = 2;
  c->id = _IDCounter;
  c->parentId = _IDCounter;
  c->timeStep = _time;
  c->previousAngle = 0.;
  c->angle = 0.;
  c->previousDivDir = Point3d( 0., 0., 0. );
  c->divDir = Point3d( 0., 0., 0. );
  c->divType = DivisionType::NONE;
  c->layerValue = 1;
  c->cellCycle = 0;
  T.addCell( c, vs );
  
  std::vector<Point2d> polygon;
  Point3d center;
  forall(const junction& j, T.S.neighbors(c))
  {
    polygon.push_back( j->tp.Pos() );
    center += Point3d( j->tp.Pos().i(), j->tp.Pos().j(), 0. );
  }
  center /= polygon.size();

  // store initial area for current cell
  c->center = center;
  c->initialArea = geometry::polygonArea(polygon);
  c->area = c->initialArea;
  c->initialLongestWallLength = ModelUtils::determineLongestWallLength( c, T );
  c->longestWallLength = c->initialLongestWallLength;
    
  lateralRoot.setPos(c->tp, center);
  
  if( _exportLineage )
    ModelExporter::exportLineageInformation( _lineageFileName, c, T, _time, false );
  
  // afterwards increment the id counter
  _IDCounter++;
  */
}

// ---------------------------------------------------------------------

void SurfaceClass::generateCell( MyTissue &T,
                                 const std::pair<double, double> &start,
                                 const std::pair<double, double> &length,
                                 const std::size_t treeId,
                                 const idPairSet &sharedJunctions,
                                 RealSurface &lateralRoot,
                                 const std::vector<Point3d> &conPoints,
                                 const bool oneCell )
{
  // set of junctions for the cell
  std::vector<junction> vs;
  
  if( !oneCell )
  {
    for( std::size_t w = 0; w < 4; w++ )
    {
      std::pair<double, double> iPos;
      
      switch(w)
      {
        case 0:
        iPos.first = start.first;
        iPos.second = start.second;
        break;
        case 1:
        iPos.first = start.first;
        iPos.second = start.second + length.second;
        break;
        case 2:
        iPos.first = start.first + length.first;
        iPos.second = start.second + length.second;
        break;
        case 3:
        iPos.first = start.first + length.first;
        iPos.second = start.second;
        break;
      }
      
      Point3d curPos = this->determinePos( iPos, conPoints );
      junction j;
      j->id = _jIDCounter;
      lateralRoot.setPos( j->tp, curPos );
      this->junctionAlreadyShared( T, j->id, j, sharedJunctions );
      vs.push_back(j);
      _jIDCounter++;
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
      this->junctionAlreadyShared( T, j->id, j, sharedJunctions );
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
  c->cellCycle = 0;
  T.addCell( c, vs );
  
  std::vector<Point2d> polygon;
  Point3d center;
  forall(const junction& j, T.S.neighbors(c))
  {
    polygon.push_back( j->tp.Pos() );
    center += Point3d( j->tp.Pos().i(), j->tp.Pos().j(), 0. );
  }
  center /= polygon.size();

  // store initial area for current cell
  c->center = center;
  c->initialArea = geometry::polygonArea(polygon);
  c->area = c->initialArea;
  c->initialLongestWallLength = ModelUtils::determineLongestWallLength( c, T );
  c->longestWallLength = c->initialLongestWallLength;
    
  lateralRoot.setPos(c->tp, center);
  
  if( _exportLineage )
    ModelExporter::exportLineageInformation( _lineageFileName, c, T, _time, false );
  
  // afterwards increment the id counter
  _IDCounter++;  
}

// ---------------------------------------------------------------------

Point3d SurfaceClass::determinePos( const std::pair<double, double> &coord,
                                    const std::vector<Point3d> &conPoints )
{
  // u and v coordinates of start should correspond to the
  // following contour points
  // u = 0, v = 0 -> conPoints.at(14)
  // u = 0, v = 1 -> conPoints.at(0)
  // u = 1, v = 0 -> conPoints.at(8)
  // u = 1, v = 1 -> conPoints.at(6)
  
  Point3d pos;
  
  // v will always be 0 or 1
  std::size_t v = (std::size_t)coord.second;
  
  double xPos = 6. * coord.first;
  std::size_t xI = (std::size_t)xPos;
  double factor = xPos - (double)xI;
  
  // left to right: between conPoints at 0 and 6
  if( v == 1 )
  {
    if( xI < 6 )
      pos = (1.-factor) * conPoints.at(xI) + factor * conPoints.at(xI+1);
    else
      pos = conPoints.at(xI);
  }
  // v == 0 -> right to left: between conPoints at 14 and 8
  else
  {
    if( xI < 6 )
      pos = (1.-factor) * conPoints.at(14-xI) + factor * conPoints.at(14-xI-1);
    else
      pos = conPoints.at(14-xI);
  }
  
  return pos;
}

// ---------------------------------------------------------------------

// generate a cell starting at the bottom left
// with different lengths
void SurfaceClass::generateCell( MyTissue &T,
                                 const std::pair<double, double> &start,
                                 const std::pair<double, double> &length,
                                 const std::size_t treeId,
                                 const idPairSet &sharedJunctions,
                                 Surface &lateralRoot )
{
  // TODO: use of mergeCells method instead of checking shared junctions
  // set of junctions for the cell
  std::vector<junction> vs;
  std::vector< std::pair<double,double> > vertices;
  vertices.resize( 4 * _lod );
  
  for( std::size_t w = 0; w < vertices.size(); w++ )
  {
    junction j;
    unsigned int uiw = (unsigned int)(w/_lod);
    double cellLength;
    double curLength = (double)w/(double)_lod - (double)uiw;
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
    
    vertices.at( w ) = std::make_pair( u, v );
    j->id = _jIDCounter;
    lateralRoot.InitPoint( j->sp, u, v );
      
    this->junctionAlreadyShared( T, j->id, j, sharedJunctions );
    
    vs.push_back(j);
    
    _jIDCounter++;
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
  c->cellCycle = 0;
  T.addCell( c, vs );
  
  std::vector<Point3d> polygon;
  Point3d center;
  forall(const junction& j, T.S.neighbors(c))
  {
    polygon.push_back( j->sp.Pos() );
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
  
  if( _exportLineage )
    ModelExporter::exportLineageInformation( _lineageFileName, c, T, _time, false );
  
  // afterwards increment the id counter
  _IDCounter++;
}

//----------------------------------------------------------------

void SurfaceClass::junctionAlreadyShared( MyTissue &T,
                                          const std::size_t jId,
                                          junction &js,
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
        return;
      }
    }
  }
}

// ---------------------------------------------------------------------
