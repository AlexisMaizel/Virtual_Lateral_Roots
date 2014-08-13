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
  //parms(section.data(), "SurfTimeScale", _surfTimeScale);
  
  // depending on data set -> TODO
  _curSurface.readTriangulation( "/home/necrolyte/Data/VLR/Virtual_Lateral_Roots/FinalVLRForMatlab/triangulation-120830_raw.txt" );
  
  //this->growStep( 0.5 );
  //_curSurface.printTriangleProperties( 0 );
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
  //std::cout << "pos: " << tp.Pos() << std::endl;
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
  Point2d pos( p.i(), p.j() );
  _curSurface.determinePosProperties( tp, pos );
  calcNormal(tp);
}

//----------------------------------------------------------------

void RealSurface::calcPos( TrianglePoint &tp )
{
  tp.u = _curSurface.checkBounds(tp.u, 0., 1.);
  tp.v = _curSurface.checkBounds(tp.v, 0., 1.);
  tp.w = _curSurface.checkBounds(tp.w, 0., 1.);
  tp.pos = _curSurface.getCoord( tp.u, tp.v, tp.w, tp.triIndex );
}

//----------------------------------------------------------------

void RealSurface::calcNormal( TrianglePoint &tp )
{
  _curSurface.determineNormal( tp );
  //std::cout << "normal: " << tp.normal << std::endl;
}

//----------------------------------------------------------------