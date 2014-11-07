// Surface class
#ifndef __SURFACE_H__
#define __SURFACE_H__

#include <string>
#include <cmath>
#include <util/parms.h>

//#include "bezier.h"
#include "BezierSurface.h"

using std::string;

const int SAMPLES = 10000;	// Sample count for function/contours

const int ARCLENGTH = 0;	// Growth types
const int CONTOUR = 1;

const double PI = 3.14159265358979323846;
const double PIx2 = PI + PI;   // 2 PI

const double DX = .0001;//.00000001;
const double MAXSEARCHSTEPS = 1000;
const int MAXSURF = 10;

class SurfacePoint {
  public:
    SurfacePoint() : u(0.), v(0.){}
    ~SurfacePoint() {}
    
    Point3d Pos();
    Point3d Pos() const;
    Point3d Normal();
    double getU() const { return u; }
    double getV() const { return v; }
    void printPos() const
    { std::cout << "Point: [" << pos.i() << ", "
      << pos.j() << ", " << pos.k() << "]" << std::endl; }

    friend std::ostream& operator<<(std::ostream& os, const SurfacePoint &sp) {
      return os;
     };
    friend std::istream& operator>>(std::istream& is, SurfacePoint &sp) {
      return is;
    };
    friend class Surface;

  private:
    double u,v;                         // Coordinates in parametric space
    Point3d pos;		// x,y,z pos of vertex
    Point3d normal;		// normalized normal vector
};

class Surface {
  public: 
    Surface(){};
    Surface( util::Parms &parms,
             string section,
             const std::string &surfaceName );
    
    void init( util::Parms &parms,
               string section,
               const std::string &surfaceName );
    
    ~Surface() {};

    // Create zero point
    void Zero(SurfacePoint &p);

    // Create inital point
    void InitPoint(SurfacePoint &p, double u, double v);
   
    // Find closest point on surface to a given point with an initial guess
    bool SetPoint(SurfacePoint &p, SurfacePoint sp, Point3d cp);

    // Advance growth by dt
    void GrowStep(double dt);

    // Compute pos of point based on current time
    void GetPos(SurfacePoint &p);

    // Distance between two surface points
    double Distance(SurfacePoint &u, SurfacePoint &v);

		double GetTime() const;

    friend class SurfacePoint;

  private:
    void CalcPos(SurfacePoint &p);      // Calc x,y,x pos from u,v values
    void CalcNormal(SurfacePoint &p);   // Calc normal

    double growthScale;			// Growth scaling constant

    int surfaces;                       // Surface for growth stages
    double surfMaxDist;                 // Max dist for closest point search
    double surfTimeScale;               // Surface time scale
    //Bezier surface[MAXSURF];            // Bezier surfaces
    BezierSurface surface[MAXSURF];            // Bezier surfaces
    double surfScale[MAXSURF];          // Surface scaling constants
    double surfTime[MAXSURF];           // Surdace time scale
    //Bezier surfCurr;                    // Current surfaces
    BezierSurface surfCurr;                    // Current surfaces

    double time;                        // Time
};

#endif
