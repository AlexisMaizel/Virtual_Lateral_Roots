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
  
  for( std::size_t p=0; p<n;++p )
   std::cout << points.at(p) << std::endl;
  
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