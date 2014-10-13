#include "SurfacePoints.h"

/**
  @file   SurfacePoints.cpp
  @brief  Class for handling triangulation of real data points
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

//----------------------------------------------------------------

SurfacePoints::SurfacePoints()
{
}

//----------------------------------------------------------------

void SurfacePoints::readTriangulation( const std::string &fileName,
                                       const double surfaceScale )
{
  std::ifstream in( fileName.c_str(), std::ifstream::in );

  in >> _minTimeStep;
  in >> _maxTimeStep;
  std::size_t timeSteps = _maxTimeStep - _minTimeStep + 1;
  _points.resize( timeSteps );
  _subsequentPoints.resize( timeSteps - 1 );
  _triangles.resize( timeSteps );
  
  for( std::size_t t = 0; t < timeSteps; t++ )
  {
    std::size_t time, pSize, tSize;
    in >> time;
    in >> pSize;
    _points.at( t ).resize( pSize );
    
    // read points of time step t
    for( std::size_t p = 0; p < pSize; p++ )
    {
      double px, py;
      in >> px >> py;
      _points.at(t).at(p) = surfaceScale*Point2d( px, py );
    }
    
    // read triangle list info
    in >> tSize;
    _triangles.at( t ).resize( tSize );
    
    for( std::size_t tr = 0; tr < tSize; tr++ )
    {
      std::size_t t1, t2, t3;
      in >> t1 >> t2 >> t3;
      _triangles.at(t).at(tr) = Point3i( t1, t2, t3 );
    }
    
    // read subsequent mapped points of time step t+1
    if( t < timeSteps - 1 )
    {
      _subsequentPoints.at( t ).resize( pSize );
      for( std::size_t p = 0; p < pSize; p++ )
      {
        double px, py;
        in >> px >> py;
        _subsequentPoints.at(t).at(p) = surfaceScale*Point2d( px, py );
      }
    }
  }
  
  in.close();
  
  std::cout << "==============================================" << std::endl;
  std::cout << "Read data with " << timeSteps << " time steps." << std::endl;
  std::cout << "First Time step: " << _points.at(0).size() << " points, " << _triangles.at(0).size() << " triangles." << std::endl;
  std::cout << "Last Time step: " << _points.at( timeSteps-1 ).size() << " points, " << _triangles.at( timeSteps-1 ).size() << " triangles." << std::endl;
  std::cout << "==============================================" << std::endl;
}

//----------------------------------------------------------------

double SurfacePoints::checkBounds( double s, double min, double max ) 
{
  if(s < min)
    return(min);
  else if(s > max)
    return(max);
  else
    return(s);
}

//----------------------------------------------------------------

void SurfacePoints::interpolate( double timeStep,
                                 std::vector<TrianglePoint> &tps )
{
  timeStep = this->checkBounds(timeStep, 0., 1.);
  
  std::size_t timeStepRange = _maxTimeStep - _minTimeStep + 1;
  
  // get range of time steps that fit the current timeStep value
  double time = timeStep * (timeStepRange-1);
  std::size_t prevTimeStep = _curTimeStep;
  _curTimeStep = (std::size_t)time;
  bool newTriangulation = prevTimeStep != _curTimeStep and _curTimeStep != 0;
  
  if( newTriangulation )
  {
    std::cout << "New timepoint:" <<  prevTimeStep << " " << _curTimeStep << std::endl;
    for(size_t i = 0; i < tps.size(); ++i)
    {
      this->getBoundaryCoord( tps.at(i), prevTimeStep );
      //if( _curTimeStep == 238 )
        //tps.at(i).printProperties();
    }
  }
  
  // the triangle list is just copied by the lower time step
  _curTriangles = _triangles.at(_curTimeStep);
  
  // then interpolate linearly between these two time steps
  // depending on the factor difference
  double factor = time - (double)_curTimeStep;
  
  //std::cout << std::endl;
  std::cout << "curT: " << _curTimeStep << " factor: " << factor << std::endl;
  
  // if the last time step is reached then just set the last entry
  if( _curTimeStep == timeStepRange-1 )
    _curPoints = _points.at(_curTimeStep);
  else
  {
    _curPoints.clear();
    for( std::size_t p = 0; p < _points.at(_curTimeStep).size(); p++ )
    {
      Point2d pos = (1.-factor) * _points.at(_curTimeStep).at(p) + _subsequentPoints.at(_curTimeStep).at(p) * factor;
      _curPoints.push_back( pos );
    }
  }
  
  if( newTriangulation )
  {
    std::cout << "determine" << std::endl;
    for(size_t i = 0; i < tps.size(); ++i)
    {
      //this->determinePosProperties(tps[i], tps[i].Pos() );
      this->determineBoundaryPosProperties(tps[i], tps[i].Pos(), _curTimeStep);
      //if( _curTimeStep == 238 )
        //tps.at(i).printProperties();
    }
  }
}

//----------------------------------------------------------------
// compute cartesian coordinates depending on the index of triangle
void SurfacePoints::getBoundaryCoord( TrianglePoint &tp, size_t timeStep)
{ 
  Point2d p1,p2,p3;
  p1 = _subsequentPoints.at(timeStep).at( _curTriangles.at( tp.triIndex ).i() - 1 );
  p2 = _subsequentPoints.at(timeStep).at( _curTriangles.at( tp.triIndex ).j() - 1 );
  p3 = _subsequentPoints.at(timeStep).at( _curTriangles.at( tp.triIndex ).k() - 1 );
  
  // determine the smallest distance between old point and current point
  // such that the u,v,w parameters are satisfied
  Point2d pos1 = tp.u * p1 + tp.v * p2 + tp.w * p3;
  Point2d pos2 = tp.v * p1 + tp.w * p2 + tp.u * p3;
  Point2d pos3 = tp.w * p1 + tp.u * p2 + tp.v * p3;

  double dist1, dist2, dist3;
  dist1 = norm(pos1 - tp.pos);
  dist2 = norm(pos2 - tp.pos);
  dist3 = norm(pos3 - tp.pos);
  
  if( dist1 <= dist2 && dist1 <= dist3 )
    tp.pos = pos1;
  else if( dist2 <= dist1 && dist2 <= dist3 )
    tp.pos = pos2;
  else
    tp.pos = pos3;
}

//----------------------------------------------------------------

// compute cartesian coordinates depending on the index of triangle
void SurfacePoints::getCoord( TrianglePoint &tp )
{ 
  Point2d p1,p2,p3;
  p1 = _curPoints.at( _curTriangles.at( tp.triIndex ).i() - 1 );
  p2 = _curPoints.at( _curTriangles.at( tp.triIndex ).j() - 1 );
  p3 = _curPoints.at( _curTriangles.at( tp.triIndex ).k() - 1 );
  
  // determine the smallest distance between old point and current point
  // such that the u,v,w parameters are satisfied
  Point2d pos1 = tp.u * p1 + tp.v * p2 + tp.w * p3;
  Point2d pos2 = tp.v * p1 + tp.w * p2 + tp.u * p3;
  Point2d pos3 = tp.w * p1 + tp.u * p2 + tp.v * p3;

  double dist1, dist2, dist3;
  dist1 = norm(pos1 - tp.pos);
  dist2 = norm(pos2 - tp.pos);
  dist3 = norm(pos3 - tp.pos);
  
  if( dist1 <= dist2 && dist1 <= dist3 )
    tp.pos = pos1;
  else if( dist2 <= dist1 && dist2 <= dist3 )
    tp.pos = pos2;
  else
    tp.pos = pos3;
  
//   if( _curTimeStep == 238 )
//   {
//     std::cout << "newpos: ";
//     tp.printPos();
//   }
}

//----------------------------------------------------------------
// determine triangle index and barycentric coordinates of triangle
void SurfacePoints::determinePosProperties( TrianglePoint &tp, const Point2d &p )
{
  // set the position vector
  tp.pos = p;
  // iterate over all triangles in order to find the one in which cp is located
  for( std::size_t t = 0; t < _curTriangles.size(); t++ )
  {
    Point2d p1,p2,p3;
    p1 = _curPoints.at( _curTriangles.at( t ).i() - 1 );
    p2 = _curPoints.at( _curTriangles.at( t ).j() - 1 );
    p3 = _curPoints.at( _curTriangles.at( t ).k() - 1 );
    
    double u,v,w;
    // check if point is in the current triangle
    if( this->pointIsInTriangle( p, p1, p2, p3, u, v, w ) )
    {
      tp.triIndex = t;
      tp.u = u;
      tp.v = v;
      tp.w = w;
      
      if( _curTimeStep == 238 )
        tp.printProperties();
      
      return;
    }
  }
  
  std::cout << "not in triangle" << std::endl;
  tp.printProperties();
}
//----------------------------------------------------------------
// determine triangle index and barycentric coordinates of triangle
void SurfacePoints::determineBoundaryPosProperties( TrianglePoint &tp,
                                                    const Point2d &p, 
                                                    const std::size_t timeStep )
{
  // set the position vector
  tp.pos = p;
  // iterate over all triangles in order to find the one in which cp is located
  for( std::size_t t = 0; t < _curTriangles.size(); t++ )
  {
    Point2d p1,p2,p3;
    p1 = _points.at(timeStep).at( _curTriangles.at( t ).i() - 1 );
    p2 = _points.at(timeStep).at( _curTriangles.at( t ).j() - 1 );
    p3 = _points.at(timeStep).at( _curTriangles.at( t ).k() - 1 );
    
    double u,v,w;
    // check if point is in the current triangle
    if( this->pointIsInTriangle( p, p1, p2, p3, u, v, w ) )
    {
      tp.triIndex = t;
      tp.u = u;
      tp.v = v;
      tp.w = w;
      
      if( _curTimeStep == 238 )
        tp.printProperties();
      
      return;
    }
  }
  
  std::cout << "boundary not in triangle" << std::endl;
  tp.printProperties();
}

//----------------------------------------------------------------

void SurfacePoints::determineNormal( TrianglePoint &tp )
{
  Point2d p1,p2,p3;
  p1 = _curPoints.at( _curTriangles.at( tp.triIndex ).i() - 1 );
  p2 = _curPoints.at( _curTriangles.at( tp.triIndex ).j() - 1 );
  p3 = _curPoints.at( _curTriangles.at( tp.triIndex ).k() - 1 );
  
  Point3d t1 = Point3d( p2.i() - p1.i(), p2.j() - p1.j(), 0. );
  Point3d t2 = Point3d( p3.i() - p1.i(), p3.j() - p1.j(), 0. );
  
  tp.normal = t1 ^ t2;
  tp.normal.normalize();
}

//----------------------------------------------------------------

bool SurfacePoints::pointIsInTriangle( const Point2d &p, const Point2d &p1,
                                       const Point2d &p2, const Point2d &p3,
                                       double &u, double &v, double &w )
{
  // compute barycentric coordinates
  double detMat = (p2.j()-p3.j())*(p1.i()-p3.i()) + (p3.i()-p2.i())*(p1.j()-p3.j());
  u = ((p2.j()-p3.j())*(p.i()-p3.i()) + (p3.i()-p2.i())*(p.j()-p3.j()))/detMat;
  v = ((p3.j()-p1.j())*(p.i()-p3.i()) + (p1.i()-p3.i())*(p.j()-p3.j()))/detMat;
  w = 1. - u - v;
  
  return ( u >= -0.0001 && v >= -0.0001 && w >= -0.0001 );
}

//----------------------------------------------------------------

void SurfacePoints::printTriangleProperties( const std::size_t timeStep )
{
  if( timeStep < _points.size() )
  {
    std::cout << "Timestep: " << timeStep << std::endl;
    std::cout << "Points: " << _points.at(timeStep).size() << std::endl;
    for( std::size_t p = 0; p < _points.at(timeStep).size(); p++ )
    {
      std::cout << "x: " << _points.at(timeStep).at(p).i()
                << " y: " << _points.at(timeStep).at(p).j() << std::endl;
    }
    std::cout << "Triangle list: " << _triangles.at(timeStep).size() << std::endl;
    for( std::size_t t = 0; t < _triangles.at(timeStep).size(); t++ )
    {
      std::cout << "[ " << _triangles.at(timeStep).at(t).i()
                << " " << _triangles.at(timeStep).at(t).j()
                << " " << _triangles.at(timeStep).at(t).k()
                << " ]" << std::endl;
    }
  }
}

//----------------------------------------------------------------
