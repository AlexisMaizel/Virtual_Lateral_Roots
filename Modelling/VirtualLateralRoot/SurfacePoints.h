#ifndef SurfacePoints_HH
#define SurfacePoints_HH

/**
  @file   SurfacePoints.h
  @brief  Class for handling triangulation of real data points
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include <geometry/geometry.h>

#include <vector>
#include <fstream>

typedef util::Vector<2, double> Point2d;
typedef util::Vector<3, std::size_t> Point3i;

class SurfacePoints
{
public:
  
  SurfacePoints();
  void readTriangulation( const std::string &fileName );
  
  Point2d getCoord( const double l1, const double l2, const double l3, 
                    const std::size_t triIndex, const std::size_t timeStep );
  
  void getBarycentricCoord( double &l1, double &l2, double &l3, const Point2d &p,
                            const std::size_t triIndex, const std::size_t timeStep );
  
  bool pointIsInTriangle( const Point2d &p, const Point2d &p0,
                          const Point2d &p1, const Point2d &p2 );
  void printTriangleProperties( const std::size_t timeStep );
  
private:
  
  std::vector< std::vector<Point2d> > _points;
  std::vector< std::vector<Point3i> > _triangles;
};

#endif // SurfacePoints_HH