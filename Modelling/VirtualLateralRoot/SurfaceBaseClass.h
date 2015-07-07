#ifndef SurfaceBaseClass_HH
#define SurfaceBaseClass_HH

/**
  @file   SurfaceBaseClass.h
  @brief  Base class for handling surface properties
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include <vector>
#include <fstream>

#include <util/parms.h>
#include <util/vector.h>

typedef util::Vector<3, double> Point3d;
typedef util::Vector<3, std::size_t> Point3i;

class SurfacePoint
{
  public:
    SurfacePoint() : u(0.), v(0.), w(0.), triIndex(0){}
    ~SurfacePoint() {}
    
    Point3d getPos() { return pos; }
    Point3d getPos() const { return pos; }
    Point3d getNormal() { return normal; }
    
    double getU() const { return u; }
    double getV() const { return v; }
    double getW() const { return w; }
    
    void setUV( const double nu, const double nv )
    { u = nu; v = nv; }
    
    std::size_t getTriIndex() const { return triIndex; }
    
    void printPos( const std::string name = "Point" ) const
    { std::cout << name << ": [" << pos.i() << ", "
      << pos.j() << ", " << pos.k() << "]" << std::endl; }

    void printProperties() const
    {
      std::cout << "pos: " << pos << std::endl;
      std::cout << "u: " << u
                << " v: " << v
                << " w: " << w
                << " t: " << triIndex << std::endl;
    }
    
    void boundParameters()
    {
      u = this->inBounds( u );
      v = this->inBounds( v );
      w = this->inBounds( w );
    }
    
    double inBounds( const double s )
    {
      if( s < 0. )
        return 0.;
      else if( s > 1. )
        return 1.;
      else
        return s;
    }
      
    friend std::ostream& operator<<(std::ostream& os, const SurfacePoint &sp) { return os; }
    friend std::istream& operator>>(std::istream& is, SurfacePoint &sp)
    { return is; };
    friend class SurfaceBaseClass;
    // TODO
    friend class Surface;
    friend class RealSurface;
    friend class SurfacePoints;

  private:
    // coordinates in parametric space
    double u,v,w;
    std::size_t triIndex;
    // x,y,z pos of vertex
    Point3d pos;
    // normalized normal vector
    Point3d normal;
};

class SurfaceBaseClass
{
public:
  
  SurfaceBaseClass(){}
  virtual void initSurfaceProperties();
  virtual void initPos( SurfacePoint &sp );
  virtual void growStep( const double dt );
  virtual void growStep( const double dt, std::vector<SurfacePoint> &sps );
  virtual void setPos( SurfacePoint &sp, const Point3d &cp );
  virtual void getPos( SurfacePoint &sp );
  
private:
  virtual void calcPos( SurfacePoint &sp );
  virtual void calcNormal( SurfacePoint &sp );
};

#endif // SurfaceBaseClass_HH