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