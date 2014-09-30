#include "RealSurface.h"

/**
  @file   RealSurface.cpp
  @brief  Class for handling triangulation of real data points
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

//----------------------------------------------------------------

RealSurface::RealSurface( util::Parms &parms, 
                          const std::string &section )
  : _time( 0. )
{
  // depending on data set -> TODO
  _curSurface.readTriangulation( "triangulation-120830_raw.txt" );
}

//----------------------------------------------------------------

void RealSurface::growStep( const double dt )
{ 
  _time += dt;
  _curSurface.interpolate( _time );
}

//----------------------------------------------------------------

void RealSurface::initPoint( TrianglePoint &tp )
{
  this->getPos( tp );
}

//----------------------------------------------------------------

void RealSurface::getPos( TrianglePoint &tp )
{
  calcPos(tp);
  calcNormal(tp);
}

//----------------------------------------------------------------

Point2d RealSurface::determinePos( const TrianglePoint &tp )
{
  return _curSurface.determineCoord( tp );
}

//----------------------------------------------------------------

void RealSurface::resetTriangleIndex( TrianglePoint &tp )
{
  _curSurface.determineTriangleIndex( tp );
}

//----------------------------------------------------------------

void RealSurface::setPos( TrianglePoint &tp, const Point3d &p )
{
  _curSurface.determinePosProperties( tp, Point2d( p.i(), p.j() ) );
  calcNormal(tp);
}

//----------------------------------------------------------------

void RealSurface::calcPos( TrianglePoint &tp )
{
  tp.u = _curSurface.checkBounds(tp.u, 0., 1.);
  tp.v = _curSurface.checkBounds(tp.v, 0., 1.);
  tp.w = _curSurface.checkBounds(tp.w, 0., 1.);
  _curSurface.getCoord( tp );
}

//----------------------------------------------------------------

std::size_t RealSurface::getCurTimeStep() const
{
  return _curSurface.getCurTimeStep();
}

//----------------------------------------------------------------

void RealSurface::calcNormal( TrianglePoint &tp )
{
  _curSurface.determineNormal( tp );
}

//----------------------------------------------------------------