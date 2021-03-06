#ifndef SurfacePoints_HH
#define SurfacePoints_HH

/**
  @file   SurfacePoints.h
  @brief  Class for handling triangulation of real data points
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SurfaceBaseClass.h"

const double eps = 0.1;

//----------------------------------------------------------------

class SurfacePoints
{
public:
  
  SurfacePoints();
  void readTriangulation( const std::string &fileName,
                          const double surfaceScale );
  
  void getCoord( SurfacePoint &sp );

  void getBoundaryCoord( SurfacePoint &sp, std::size_t timeStep );
  
  void determinePosProperties( SurfacePoint &sp, const Point3d &p );

  void determineBoundaryPosProperties( SurfacePoint &sp, const Point3d &p, 
                                       const std::size_t timeStep );
  
  void determineNormal( SurfacePoint &sp );
  
  bool pointIsInTriangle( const Point3d &p, const Point3d &p0,
                          const Point3d &p1, const Point3d &p2,
                          double &u, double &v, double &w );
  
  void printTriangleProperties( const std::size_t timeStep );
  
  void interpolate( double timeStep,
                    std::vector<SurfacePoint> &sps );

  std::size_t getNumTriangles() const
  { return _curTriangles.size(); }
  
  std::size_t getCurTimeStep() const
  { return _curTimeStep; }
  
  std::size_t getMaxTimeStep() const
  { return _maxTimeStep; }
  
private:
  
  // these are the 2D points for each cell position in time step t
  std::vector< std::vector<Point3d> > _points;
  // these are the mapped 2D points for each cell position in time step t+1
  // Note that if a cell bas been dividing in t+1 then we compute the average 
  //of the daughter's positions to get an unique mapping of cells from t to t+1
  std::vector< std::vector<Point3d> > _subsequentPoints;
  // these are the triangles of time step t generated by the list of indexes
  // to access _points
  std::vector< std::vector<Point3i> > _triangles;
  
  // these are the current 2D points for each cell position in the current time step
  std::vector<Point3d> _curPoints;
  // this is the current triangle list
  std::vector<Point3i> _curTriangles;
  
  std::size_t _minTimeStep;
  std::size_t _maxTimeStep;
  
  std::size_t _curTimeStep;
};

#endif // SurfacePoints_HH
