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

void SurfacePoints::readTriangulation( const std::string &fileName )
{
  std::ifstream in( fileName.c_str(), std::ifstream::in );

  in >> _maxTimeStep;
  _points.resize( _maxTimeStep );
  _subsequentPoints.resize( _maxTimeStep-1 );
  _triangles.resize( _maxTimeStep );
  
  for( std::size_t t = 0; t < _maxTimeStep; t++ )
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
      _points.at(t).at(p) = Point2d( px, py );
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
    if( t < _maxTimeStep - 1 )
    {
      _subsequentPoints.at( t ).resize( pSize );
      for( std::size_t p = 0; p < pSize; p++ )
      {
        double px, py;
        in >> px >> py;
        _subsequentPoints.at(t).at(p) = Point2d( px, py );
      }
    }
  }
  
  in.close();
}

//----------------------------------------------------------------

void SurfacePoints::interpolate( double timeStep )
{
  if( timeStep > 1. )
    timeStep = 1.;
  if( timeStep < 0. )
    timeStep = 0.;
  
  // get range of time steps that fit the current timeStep value
  double time = timeStep * (_maxTimeStep-1);
  std::size_t startT = (std::size_t)time;
  
  // the triangle list is just copied by the lower time step
  _curTriangles = _triangles.at(startT);
  
  // then interpolate linearly between these two time steps
  // depending on the factor difference
  double factor = time - (double)startT;
  
  // if the last time step is reached then just set the last entry
  if( startT == _maxTimeStep-1 )
    _curPoints = _points.at(startT);
  else
  {
    _curPoints.clear();
    for( std::size_t p = 0; p < _points.at(startT).size(); p++ )
    {
      Point2d pos = factor * _points.at(startT).at(p) + _subsequentPoints.at(startT).at(p) * (1.-factor);
      pos /= 2.;
      _curPoints.push_back( pos );
    }
  }
}

//----------------------------------------------------------------

// compute cartesian coordinates depending on the index of triangle
Point2d SurfacePoints::getCoord( const double l1, const double l2, const double l3, 
                                 const std::size_t triIndex, const std::size_t timeStep )
{
  Point3i indexes = _triangles.at(timeStep).at(triIndex);
  Point2d p1 = _points.at(timeStep).at( indexes.i() );
  Point2d p2 = _points.at(timeStep).at( indexes.j() );
  Point2d p3 = _points.at(timeStep).at( indexes.k() );
  
  return ( l1 * p1 + l2 * p2 + l3 * p3 );
}

//----------------------------------------------------------------

// compute area coordinates depending on the index of triangle
void SurfacePoints::getBarycentricCoord( double &l1, double &l2, double &l3,
                                         const Point2d &p,
                                         const std::size_t triIndex,
                                         const std::size_t timeStep )
{
  Point3i indexes = _triangles.at(timeStep).at(triIndex);
  Point2d p1 = _points.at(timeStep).at( indexes.i() );
  Point2d p2 = _points.at(timeStep).at( indexes.j() );
  Point2d p3 = _points.at(timeStep).at( indexes.k() );
  
  // compute barycentric coordinates
  double detMat = (p2.j()-p3.j())*(p1.i()-p3.i()) + (p3.i()-p2.i())*(p1.j()-p3.j());
  l1 = ((p2.j()-p3.j())*(p.i()-p3.i()) + (p3.i()-p2.i())*(p.j()-p3.j()))/detMat;
  l2 = ((p3.j()-p1.j())*(p.i()-p3.i()) + (p1.i()-p3.i())*(p.j()-p3.j()))/detMat;
  l3 = 1. - l1 - l2;
}

//----------------------------------------------------------------

bool SurfacePoints::pointIsInTriangle( const Point2d &p, const Point2d &p0,
                                       const Point2d &p1, const Point2d &p2 )
{
  double area = 0.5*(-p1.j()*p2.i() + p0.j()*(-p1.i() + p2.i()) + p0.i()*(p1.j() - p2.j()) + p1.i()*p2.j());
  double s = 1./(2.*area)*(p0.j()*p2.i() - p0.i()*p2.j() + (p2.j() - p0.j())*p.i() + (p0.i() - p2.i())*p.j());
  double t = 1./(2.*area)*(p0.i()*p1.j() - p0.j()*p1.i() + (p0.j() - p1.j())*p.i() + (p1.i() - p0.i())*p.j());
  
  return (s > 0. && t > 0. && 1.-s-t > 0);
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
