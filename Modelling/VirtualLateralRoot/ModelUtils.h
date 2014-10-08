#ifndef ModelUtils_HH
#define ModelUtils_HH

#include "ModelHeader.h"

/**
  @file   ModelUtils.h
  @brief  Contains namespace for util methods used in the model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

namespace ModelUtils
{
  double determineLongestWallLength( const cell& c,
                                     const MyTissue& T );
  
  DivisionType::type determineDivisionType( const MyTissue::division_data& ddata,
                                            const double angleThreshold );
  
  double getDivisionAngle( const MyTissue::division_data& ddata );
  
  double determineDivisionAngle( const MyTissue::division_data& ddata );
  
  std::vector<Point2d> determineConvexHull( const cell &c, const MyTissue& T );
  
  bool pointInHull( const Point2d &p, const std::vector<Point2d> &hull );
  
  bool doIntersect( const Point2d &p1, const Point2d &q1,
                    const Point2d &p2, const Point2d &q2 );
  
  bool onSegment( const Point2d &p1, const Point2d &p2, const Point2d &p3 );
  
  int orientation( const Point2d &p1, const Point2d &p2, const Point2d &p3 );
  
  double cross( const Point2d &p1, const Point2d &p2, const Point2d &p3 );
  
  std::vector<Point3d> loadContourPoints( const std::string &fileName,
                                          const double surfaceScale );
  
  bool xSort( const Point2d &p1, const Point2d &p2 );
}

#endif // ModelUtils_HH