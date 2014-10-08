#include "ModelUtils.h"

/**
  @file   ModelUtils.cpp
  @brief  Contains namespace for util methods used in the model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

namespace ModelUtils
{
  
// ---------------------------------------------------------------------

double determineLongestWallLength( const cell& c,
                                   const MyTissue& T )
{
  double maxLength = 0.;
  forall( const junction& j,T.S.neighbors(c) )
  {
    const junction& jn = T.S.nextTo(c, j);
    const Point3d& jpos = j->sp.Pos();
    const Point3d& jnpos = jn->sp.Pos();
    
    double length = 0.;
    for( std::size_t l = 0; l<3; l++ )
      length += (jpos[l] - jnpos[l])*(jpos[l] - jnpos[l]);
    
    if( length >= maxLength )
      maxLength = length;
  }
  return maxLength;
}

// ---------------------------------------------------------------------

DivisionType::type determineDivisionType( const MyTissue::division_data& ddata,
                                          const double angleThreshold )
{
  // get pair of points of division wall
  Point3d u = ddata.pu;
  Point3d v = ddata.pv;
  
  // y axis
  Point3d yaxisDir = Point3d( 0., 1., 0. );
  Point3d dir;
  
  if( u.j() <= v.j() )
    dir = v-u;
  else
    dir = u-v;
  
  dir.normalize();
  double angle = 180./M_PI * acos( dir*yaxisDir );
  // anticlinal
  if( angle <= angleThreshold )
    return DivisionType::ANTICLINAL;
  // periclinal
  else
    return DivisionType::PERICLINAL;
}

// ---------------------------------------------------------------------

double getDivisionAngle( const MyTissue::division_data& ddata )
{
  // get pair of points of division wall
  Point3d u = ddata.pu;
  Point3d v = ddata.pv;
  
  // y axis
  Point3d yaxisDir = Point3d( 0., 1., 0. );
  
  if( u.i() <= v.i() )
  {
    if( u.j() <= v.j() )
    {
      Point3d dir = v-u;
      dir.normalize();
      return ( 2.*M_PI - acos( dir*yaxisDir ) );
    }
    else
    {
      Point3d dir = u-v;
      dir.normalize();
      return ( acos( dir*yaxisDir ) );
    }
  }
  else
  {
    if( u.j() <= v.j() )
    {
      Point3d dir = v-u;
      dir.normalize();
      return ( acos( dir*yaxisDir ) );
    }
    else
    {
      Point3d dir = u-v;
      dir.normalize();
      return ( 2.*M_PI - acos( dir*yaxisDir ) );
    }
  }
}

// ---------------------------------------------------------------------

std::vector<Point2d> determineConvexHull( const cell &c, const MyTissue& T )
{
  std::vector<Point2d> points;
  // determine points
  forall(const junction& j, T.S.neighbors(c))
    points.push_back( j->tp.Pos() );
  
  std::size_t n = points.size();
  
  if( n < 3 )
    return points;
  
  std::size_t k = 0;
  std::vector<Point2d> H(2*n);
 
  // sort points lexicographically
  std::sort( points.begin(), points.end(), ModelUtils::xSort );
  
  // lower hull
  for(std::size_t i = 0; i < n; ++i)
  {
    while(k >= 2 && cross( H.at(k-2), H.at(k-1), points.at(i) ) <= 0)
      k--;
  
    H.at(k++) = points.at(i);
  }
 
  // upper hull
  std::size_t t = k+1;
  for(int i = n-2; i >= 0; i--)
  {
    while(k >= t && cross( H.at(k-2), H.at(k-1), points.at(i) ) <= 0)
      k--;
      
    H.at(k++) = points.at(i);
  }
  
  H.resize(k);
  return H;
}

// ---------------------------------------------------------------------

bool pointInHull( const Point2d &p, const std::vector<Point2d> &hull )
{
  // use ray casting check: if the number of intersections is odd -> inside hull
  // if even -> outside of hull
  
  if( hull.size() == 0 )
    return false;
  
  // create edges of point pairs of hull
  Point2d rightMostPoint( hull.at(0) );
  std::vector< std::pair<Point2d,Point2d> > edges;
  for( std::size_t p=0;p<hull.size()-1; ++p )
  {
    if( hull.at(p+1).i() > rightMostPoint.i() )
      rightMostPoint = hull.at(p+1);
    
    // ignore horizontal edges
    if( hull.at(p).j() != hull.at(p+1).j() )
      edges.push_back( std::make_pair( hull.at(p), hull.at(p+1) ) );
  }
  
  rightMostPoint.i() += 5.;
  
  std::size_t count = 0;
  for( std::size_t e=0; e< edges.size(); ++e )
  {
    if( doIntersect( p, rightMostPoint, edges.at(e).first, edges.at(e).second ) )
      count++;
  }
  
  if( count%2 )
    return true;
  else
    return false;
}

// ---------------------------------------------------------------------

bool doIntersect( const Point2d &p1, const Point2d &q1,
                  const Point2d &p2, const Point2d &q2 )
{
  int o1 = orientation( p1, q1, p2 );
  int o2 = orientation( p1, q1, q2 );
  int o3 = orientation( p2, q2, p1 );
  int o4 = orientation( p2, q2, q1 );

  // General case
  if(o1 != o2 && o3 != o4)
    return true;

  // Special Cases
  // p1, q1 and p2 are colinear and p2 lies on segment p1q1
  if(o1 == 0 && onSegment(p1, p2, q1)) return true;

  // p1, q1 and p2 are colinear and q2 lies on segment p1q1
  if(o2 == 0 && onSegment(p1, q2, q1)) return true;

  // p2, q2 and p1 are colinear and p1 lies on segment p2q2
  if(o3 == 0 && onSegment(p2, p1, q2)) return true;

    // p2, q2 and q1 are colinear and q1 lies on segment p2q2
  if(o4 == 0 && onSegment(p2, q1, q2)) return true;

  return false; // Doesn't fall in any of the above cases
}

// ---------------------------------------------------------------------

bool onSegment( const Point2d &p1, const Point2d &p2, const Point2d &p3 )
{
  if( p2.i() <= std::max(p1.i(), p3.i()) && p2.i() >= std::min(p1.i(), p3.i()) &&
      p2.j() <= std::max(p1.j(), p3.j()) && p2.j() >= std::min(p1.j(), p3.j()))
    return true;
 
  return false;  
}

// ---------------------------------------------------------------------

int orientation( const Point2d &p1, const Point2d &p2, const Point2d &p3 )
{
  int val = (p2.j() - p1.j()) * (p3.i() - p2.i()) - (p2.i() - p1.i()) * (p3.j() - p2.j());
 
  if(val == 0)
    return 0;  // colinear
 
  return (val > 0)? 1: 2; // clock or counterclock wise
}

// ---------------------------------------------------------------------

double cross( const Point2d &p1, const Point2d &p2, const Point2d &p3 )
{
  return (p2.i() - p1.i()) * (p3.j() - p1.j()) - (p2.j() - p1.j()) * (p3.i() - p1.i());
}

// ---------------------------------------------------------------------

bool xSort( const Point2d &p1, const Point2d &p2 )
{
  if( p1.i() != p2.i() )
    return (p1.i() < p2.i());
  else
    return ( p1.j() <= p2.j() );
}

// ---------------------------------------------------------------------

std::vector<Point3d> loadContourPoints( const std::string &fileName,
                                        const double surfaceScale )
{
  std::ifstream in( fileName.c_str(), std::ifstream::in );

  std::vector<Point3d> conPoints;  
  std::size_t numPoints;
  in >> numPoints;
  conPoints.resize(numPoints);
  for( std::size_t p = 0; p< numPoints; ++p )
  {
    Point3d pos(0., 0., 0.);
    in >> pos.i() >> pos.j();
    conPoints.at(p) = surfaceScale*pos;
  }
    
  in.close();
  return conPoints;
}

// ---------------------------------------------------------------------

double determineDivisionAngle( const MyTissue::division_data& ddata )
{
  // get pair of points of division wall
  Point3d u = ddata.pu;
  Point3d v = ddata.pv;
  
  // x axis
  Point3d xaxisDir = Point3d( 1., 0., 0. );
  Point3d dir;
  
  if( u.j() <= v.j() )
    dir = v-u;
  else
    dir = u-v;
  
  dir.normalize();
  return ( 180./M_PI * acos( dir*xaxisDir ) );
}

// ---------------------------------------------------------------------

}