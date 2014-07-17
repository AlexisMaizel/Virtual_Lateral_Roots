#ifndef SurfacePoints_HH
#define SurfacePoints_HH

/**
  @file   SurfacePoints.h
  @brief  Class for handling triangulation of real data points
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include <geometry/geometry.h>

class SurfacePoints
{
public:
  
  SurfacePoints();
  bool pointIsInTriangle( const Point2d &p, const Point2d &p0,
                          const Point2d &p1, const Point2d &p2 );
  
private:
  
};

#endif // SurfacePoints_HH