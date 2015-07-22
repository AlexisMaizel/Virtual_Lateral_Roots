#ifndef __SURFACE_H__
#define __SURFACE_H__

#include "SurfaceBaseClass.h"
#include "bezier.h"

const double DX = .0000001;//.0000001;
const double MAXSEARCHSTEPS = 1000;
const int MAXSURF = 1000;

class Surface : public SurfaceBaseClass
{
  public: 
    Surface() : _time( 0. ) {}
    
    void init( util::Parms &parms,
               string section,
               const bool bezierGrowthSurface,
               const bool interpolateBezierSurfaces,
               const bool onlyGrowthInHeight,
               const std::string &surfaceName,
               const std::size_t highOrderPattern );
    
    ~Surface() {};

    void determineSurfaceHeaderProperties( const std::string bezFile );
    
    void applyGrowthOnlyInHeight( Bezier &surfaceS, const Bezier &surfaceE );
    
    // Create inital point
    void initPos( SurfacePoint &sp );
   
    // Find closest point on surface to a given point with an initial guess
    void setPos( SurfacePoint &sp, const Point3d &cp );

    // Advance growth by dt
    void growStep( const double dt );

    // Compute pos of point based on current time
    void getPos( SurfacePoint &sp );

    conpoi getCurrentControlPoints() const
    { return surfCurr.getControlPoints(); }
    
    void increaseDomeTipHeight( Bezier &surface );
    
    void applyControlpointsVariation( Bezier &surface,
                                      const bool focusTop );
    
    friend class SurfacePoint;

  private:
    void calcPos( SurfacePoint &sp );
    void calcNormal( SurfacePoint &sp );

    double growthScale;			// Growth scaling constant

    int surfaces;                       // Surface for growth stages
    double surfMaxDist;                 // Max dist for closest point search
    double surfTimeScale;               // Surface time scale
    Bezier surface[MAXSURF];            // Bezier surfaces
    //BezierSurface surface[MAXSURF];            // Bezier surfaces
    double surfScale[MAXSURF];          // Surface scaling constants
    double surfTime[MAXSURF];           // Surdace time scale
    Bezier surfCurr;                    // Current surfaces
    //BezierSurface surfCurr;                    // Current surfaces
    
    double _time;
    bool _bezierGrowthSurface;
    std::size_t _numSurfaces;
    std::size_t _numPatches;
    std::size_t _numControlPoints;
};

#endif
