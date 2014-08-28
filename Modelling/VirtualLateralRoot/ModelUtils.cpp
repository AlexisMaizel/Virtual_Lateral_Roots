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

void determineConvexHull( const cell &c, const MyTissue& T )
{
  std::vector<Point2d> points;
  // determine points
  forall(const junction& j, T.S.neighbors(c))
    points.push_back( j->tp.Pos() );
  
  if( points.size() < 1 )
    return;
  
  // sort after y values
  std::sort( points.begin(), points.end(), ModelUtils::ySort );
  
  Point2d ymin = points.at(0);
  
  for( std::size_t i = 0; i < points.size(); ++i )
    std::cout << points.at(i) << std::endl;
  
  // TODO
  /*
  // compute orientation for all points with respect to
  // the first one
  for( std::size_t i = 1; i < points.size(); ++i )
  {
    Point2d pT = points.at(i);
    double aT = dotProduct( ymin, pT );
    std::size_t kT = i;
    for( std::size_t j = 0; j < points.size()-2; ++j )
    {
      if( dotProduct( ymin, points.at(i+j) ) < aT )
      {
        pT = points.at(i+j);
        aT = dotProduct( ymin, points.at(i+j) );
        kT = i+j;
      }
    }
    points.at(kT) = points.at(i);
    points.at(i) = pT;
  }
  */
  //points.at(0) = points.at(points.size()-1);
}

// ---------------------------------------------------------------------

bool ySort( const Point2d &p1, const Point2d &p2 )
{
  if( p1.j() != p2.j() )
    return (p1.j() < p2.j());
  else
    return ( p1.i() <= p2.i() );
}

// ---------------------------------------------------------------------

double dotProduct( const Point2d &p1, const Point2d &p2 )
{
  return (p2.i()-p1.i())/std::sqrt((p2.i()-p1.i())*(p2.i()-p1.i()) + (p2.j()-p1.j())*(p2.j()-p1.j()));
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