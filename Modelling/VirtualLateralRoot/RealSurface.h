#ifndef RealSurface_HH
#define RealSurface_HH

#include "SurfaceBaseClass.h"
#include "SurfacePoints.h"

/**
  @file   RealSurface.h
  @brief  Class for handling triangulation of real data points
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

class RealSurface : public SurfaceBaseClass
{
public:
  RealSurface( util::Parms &parms, 
               const std::string &section );
  
  void init( const double surfaceScale,
             const std::string &fileName,
             const bool useAutomaticContourPoints );
  
  void initPos( SurfacePoint &sp );
  
  void growStep( const double dt,
                 std::vector<SurfacePoint> &sps );
  
  void getPos( SurfacePoint &sp );
  
  void setPos( SurfacePoint &sp, const Point3d &p );
  
  void calcPos( SurfacePoint &sp );
  
  void calcNormal( SurfacePoint &sp );
  
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
