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
}

//----------------------------------------------------------------

void RealSurface::init( const double surfaceScale,
                        const std::string &fileName,
                        const bool useAutomaticContourPoints )
{
  std::string name = "triangulation-";
  name += fileName;
  if( useAutomaticContourPoints )
    name += "_auto.txt";
  else
    name += ".txt";
  _curSurface.readTriangulation( name, surfaceScale );
}

//----------------------------------------------------------------

void RealSurface::growStep( const double dt,
                            std::vector<TrianglePoint> &tps )
{ 
  _time += dt;
  _curSurface.interpolate( _time , tps );
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

void RealSurface::calcNormal( TrianglePoint &tp )
{
  _curSurface.determineNormal( tp );
}

//----------------------------------------------------------------
