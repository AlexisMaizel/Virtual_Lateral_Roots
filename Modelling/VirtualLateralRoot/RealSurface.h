#ifndef RealSurface_HH
#define RealSurface_HH

#include <util/parms.h>
#include <util/vector.h>

#include "SurfacePoints.h"
#include "ModelHeader.h"

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
  
  void init( const double surfaceScale,
             const std::string &fileName,
             const bool useAutomaticContourPoints );
  
  void initPoint( TrianglePoint &tp );
  
  void growStep( const double dt,
                 std::vector<TrianglePoint> &tps );
  
  void getPos( TrianglePoint &tp );
  
  void setPos( TrianglePoint &tp, const Point3d &p );
  
  void calcPos( TrianglePoint &tp );
  
  void calcNormal( TrianglePoint &tp );
  
  std::size_t getCurTimeStep() const
  { return _curSurface.getCurTimeStep(); }
  
  std::size_t getMaxTimeStep() const
  { return _curSurface.getMaxTimeStep(); }
  
private:
  double _time;
  double _surfTimeScale;
  double _maxTime;
  
  SurfacePoints _curSurface;
};

#endif // RealSurface_HH
