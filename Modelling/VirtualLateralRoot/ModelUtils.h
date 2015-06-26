#ifndef ModelUtils_HH
#define ModelUtils_HH

#include "ModelHeader.h"

/**
  @file   ModelUtils.h
  @brief  Contains namespace for util methods used in the model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

struct differentCell
{
  bool operator()(const cell &c1, const cell &c2) const
  {
    return( c1->id < c2->id );
  }
};

namespace ModelUtils
{
  double determineLongestWallLength( const cell& c,
                                     const MyTissue& T );
  
  void findAppropriateJunctionPoint( Point3d &p,
                                   const std::set<junction> &juncs,
                                   const MyTissue& T,
                                   const cell& c,
                                   const bool clockwise );
  
  bool findPointInJunctionSet( const Point3d &p,
                             const std::set<junction> &juncs );
  
  void findNearestPointToMerge( MyTissue &T, junction &js );
  
  bool equalPoints( const Point3d &p1, const junction &j2 );
  
  bool equalPoints( const Point3d &p1, const Point3d &p2 );
  
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
  
  Point3d getCenterAfterApplyingLODToCell( const cell &c, const MyTissue& T,
                                           double &area, const double eps );
  
  bool shouldNodeBeErased( const Point3d &nodePos1,
                           const Point3d &nodePos2,
                           const double eps );
  
  std::vector<MyTissue::division_data> determinePossibleDivisionData(
                                      const cell& c,
                                      const double epsLength,
                                      const double epsLOD,
                                      const MyTissue& T );
  
  std::vector<MyTissue::division_data> determinePossibleDivisionDataAndPreserveEqualArea(
                                      const cell& c,
                                      const double epsLength,
                                      const double epsLOD,
                                      const MyTissue& T,
                                      const double equalAreaRatio );
  
  bool checkDivisionArea( const cell& c,
                          const MyTissue::division_data &ddata,
                          const MyTissue& T,
                          const double equalAreaRatio );
  
  std::vector<Point3d> loadContourPoints( const std::string &fileName,
                                          const double surfaceScale );
  
  bool xSort( const Point2d &p1, const Point2d &p2 );
  
  double getNormalDistribution( const double x, const double mu,
                                const double sd );
  
  double getSD( const std::vector<double> &vals,
                const double mu );
  
  std::set<junction> determineNeedlessJunctions( const cell &c,
                                                 const MyTissue& T,
                                                 Point3d &center,
                                                 const double eps );
  
  void drawBezierCurve( const std::set<Point3d, lessXPos> &cps,
                        const util::Palette::Color &color,
                        GLUquadricObj *quadratic );
  
  void drawLine( const Point3d &pos1, const Point3d &pos2,
                 const util::Palette::Color &color,
                 const double lineWidth = 2.5 );
  
  void drawCylinder( const Point3d &pos1, const Point3d &pos2,
                     const util::Palette::Color &color,
                     GLUquadricObj *quadratic,
                     const double r, std::size_t slices,
                     std::size_t stacks );
  
  void drawControlPoint( const Point3d &pos,
                         const util::Palette::Color &color );
  
  void drawBezierSurface( const conpoi &cps,
                          const util::Palette::Color &color );
  
  void splitNonAdjacentCells( std::set<cell, differentCell> curve,
                              const MyTissue& T,
                              std::vector< std::set<cell, differentCell> > &curves );
  
  void drawSphere( const Point3d &pos, double r, std::size_t lats,
                   std::size_t longs, const util::Palette::Color &color,
                   GLUquadricObj *quadratic );
  
  Point3d computeBezierPoint( const conpoi &cps,
                              const double u,
                              const double v );
  
  Point3d computeBezierPoint( const vector<Point3d> &cps,
                              const double t );
  
  double determineSlope( const Point3d &dir );
  
  bool emptyIntersection( const std::set<cell, differentCell> &s1,
                          const std::set<cell, differentCell> &s2 );
  
  int binom( unsigned int n, unsigned int k );
  
  int binomR( unsigned int n, unsigned int k );
  
  std::size_t getRandomResultOfDistribution( const std::vector<double> &probs );
}

#endif // ModelUtils_HH