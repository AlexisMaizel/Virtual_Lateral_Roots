#ifndef RealSurface_HH
#define RealSurface_HH

#include <util/parms.h>

#include "SurfacePoints.h"

/**
  @file   RealSurface.h
  @brief  Class for handling triangulation of real data points
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

class RealSurface
{
public:
  RealSurface( util::Parms &parms, 
               const std::string &section );
  
  void growStep( const double dt );
private:
  double _time;
  double _surfTimeScale;
  double _maxTime;
  
  SurfacePoints _curSurface;
};

#endif // RealSurface_HH