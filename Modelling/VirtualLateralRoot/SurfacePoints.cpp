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

bool SurfacePoints::pointIsInTriangle( const Point2d &p, const Point2d &p0,
                                       const Point2d &p1, const Point2d &p2 )
{
  double area = 0.5*(-p1.j()*p2.i() + p0.j()*(-p1.i() + p2.i()) + p0.i()*(p1.j() - p2.j()) + p1.i()*p2.j());
  double s = 1./(2.*area)*(p0.j()*p2.i() - p0.i()*p2.j() + (p2.j() - p0.j())*p.i() + (p0.i() - p2.i())*p.j());
  double t = 1./(2.*area)*(p0.i()*p1.j() - p0.j()*p1.i() + (p0.j() - p1.j())*p.i() + (p1.i() - p0.i())*p.j());
  
  return (s > 0. && t > 0. && 1.-s-t > 0);
}

//----------------------------------------------------------------
