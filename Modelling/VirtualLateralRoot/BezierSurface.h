#ifndef BezierSurface_HH
#define BezierSurface_HH

#include <string>
#include <vector>

#include <util/parms.h>
#include <util/vector.h>

/**
  @file   BezierSurface.h
  @brief  Class for handling a bezier surface
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

typedef util::Vector<3, double> Point3d;

class BezierSurface
{
  public:
    BezierSurface() {};
    ~BezierSurface() {};

    void Load( const std::string &bezFile );
    
    // Copy and/or scale
    void Scale( BezierSurface &b, double Scale );

    // Create new interpolated (and scaled) surface
    void Interpolate( BezierSurface &src1, BezierSurface &src2, 
                      double scale1, double scale2, double s );

    // Return x,y,z point from u,v parameters
    Point3d EvalCoord( double u, double v );

    // Print points (debugging)
    void Print();

  private:
    int binom( unsigned int n, unsigned int k );

    // matrix storing unique control points
    std::vector< std::vector<Point3d> > _cpMatrix;
};

#endif // BezierSurface_HH