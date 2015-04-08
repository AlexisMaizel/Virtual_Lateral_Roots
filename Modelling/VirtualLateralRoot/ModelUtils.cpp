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

Point3d getCenterAfterApplyingLODToCell( const cell &c, const MyTissue& T,
                                         const std::size_t surfaceType,
                                         double &area,
                                         const double eps )
{
  Point3d center = Point3d( 0., 0., 0. );
  area = 0.;
  
  // first create linked information of junctions such that
  // later each triple of junctions can be compared to apply a lod
  std::vector<Point3d> juncs;
  forall( const junction& j, T.S.neighbors(c) )
  {
    if( surfaceType == 0 )
      juncs.push_back( j->sp.Pos() );
    else
      juncs.push_back( Point3d( j->tp.Pos().i(), j->tp.Pos().j(), 0. ) );
    
    center += juncs.back();
  }
  
  center /= juncs.size();
  area = geometry::polygonArea(juncs);
  
  // generate vector for storing node indexes which will be erased
  std::set<std::size_t> eraseNodes;
  // new vector of junctions for which later the center of mass is computed
  std::vector<Point3d> newJuncs;
  
  // loop over path making an edge criterion check for each triple of nodes
  // Note that we have a ring of a cell which means that for example for a 
  // cell with 4 junctions we also execute 4 checks of the edge criterion
  bool lastIndexWillBeErased = false;
  std::size_t erasedNodesInARow = 0;
  for( std::size_t p=0; p<juncs.size(); p++ )
  {
    std::size_t erasedIndex;
    
    // initialize ab and ac
    Point3d ac;
    Point3d ab;
    Point3d a,b,c;
    
    // generate direction vector AC and AB
    if( p == juncs.size()-2 )
    {
      a = juncs.at( p );
      b = juncs.at( p+1 );
      c = juncs.at( 0 );
      
      if( lastIndexWillBeErased )
      {
        lastIndexWillBeErased = false;
        a = juncs.at( p-erasedNodesInARow );
      }
      
      ac = c - a;
      ab = b - a;
      erasedIndex = p+1;
      //std::cout << "pl: " << a << " pm: " << b << " pr: " << c << std::endl;          
    }
    else if( p == juncs.size()-1 )
    {
      a = juncs.at( p );
      b = juncs.at( 0 );
      
      if( lastIndexWillBeErased )
      {
        lastIndexWillBeErased = false;
        a = juncs.at( p-erasedNodesInARow );
      }
      
      // for this special case check if the index=1 was not already
      // erased; if so then consider the next position which is not
      // erased
      std::size_t index = 1;
      while( eraseNodes.find(index) != eraseNodes.end() )
        index++;
      
      c = juncs.at( index );
      
      ac = c - a;
      ab = b - a;
      erasedIndex = 0;
      //std::cout << "pl: " << a << " pm: " << b << " pr: " << c << std::endl;
    }
    else
    {
      a = juncs.at( p );
      b = juncs.at( p+1 );
      c = juncs.at( p+2 );
      
      if( lastIndexWillBeErased )
      {
        lastIndexWillBeErased = false;
        a = juncs.at( p-erasedNodesInARow );
      }
      
      ac = c - a;
      ab = b - a;
      erasedIndex = p+1;
      //std::cout << "pl: " << a << " pm: " << b << " pr: " << c << std::endl;
    }
    
    // if node B does not pass test then add to vector eraseNodes
    if( shouldNodeBeErased( ac, ab, eps ) )
    {
      //std::cout << "erasedIndex: " << erasedIndex << std::endl;
      eraseNodes.insert( erasedIndex );
      lastIndexWillBeErased = true;
      
      if( erasedIndex - 1 == p )
        erasedNodesInARow++;
    }
    else
      erasedNodesInARow = 0;
  }
  
  if( eraseNodes.size() > 0 )
  {
    center = Point3d( 0., 0., 0. );
    area = 0.;
    
    for( std::size_t i=0; i<juncs.size(); i++ )
    {
      if( eraseNodes.find(i) == eraseNodes.end() )
      {
        newJuncs.push_back( juncs.at(i) );
        center += newJuncs.back();
      }
    }
    
    center /= newJuncs.size();
    area = geometry::polygonArea(newJuncs);
  }
  
  //std::cout << "center: " << center << std::endl;
  return center;
}

//----------------------------------------------------------------------------

bool shouldNodeBeErased( const Point3d &nodePos1,
                         const Point3d &nodePos2,
                         const double eps )
{
  // edge criterion:
  // shortest distance between AC=nodePos1 and node B which is
  // between node A and C in the path is calculated by:
  // shortest dist = norm( AC x AB ) / norm( AC )

  // temporal numerator vector
  Point3d numerator;

  // calculate cross product of AC=nodePos1 and AB=nodePos2
  // and write result in numerator
  numerator = nodePos1 ^ nodePos2;

  // compute norm of above result and of second
  double numeratorResult = norm( numerator );
  double denominatorResult = norm( nodePos1 );
  
  // result is shortest distance between line AC=nodePos1 and point B
  double result = numeratorResult/denominatorResult;
  //std::cout << "result: " << result << std::endl;

  // if the result is zero that means the direction vectors are colinear
  // then set the result to true because the inner node is then not required
  if( result == 0 )
    return true;

  if(result < eps)
    return true;
  else
    return false;
}

// ---------------------------------------------------------------------

std::vector<MyTissue::division_data> determinePossibleDivisionData(
                                      const cell& c,
                                      const std::size_t surfaceType,
                                      const double eps )
{
  std::vector<MyTissue::division_data> divisionData;
  double deltaAngle = 5.;
  for( double angle = 0.; angle < 180.; angle += deltaAngle )
  {
    MyTissue::division_data ddata;
    const Point3d& center = c->center;
    double a = M_PI/180. * angle;
    Point3d direction = Point3d(-sin(a), cos(a), 0);
    forall( const junction& j,T.S.neighbors(c) )
    {
      Point3d jpos, jnpos;
      const junction& jn = T.S.nextTo(c, j);
      if( surfaceType == 0 )
      {
        jpos = j->sp.Pos();
        jnpos = jn->sp.Pos();
      }
      else
      {
        jpos = Point3d( j->tp.Pos().i(), j->tp.Pos().j(), 0. );
        jnpos = Point3d( jn->tp.Pos().i(), jn->tp.Pos().j(), 0. );
      }
      Point3d u;
      double s;
      if(geometry::planeLineIntersection(u, s, center, direction, jpos, jnpos) and s >= 0 and s <= 1)
      {
        if((jpos - center)*direction > 0)
        {
          ddata.v1 = j;
          ddata.pv = u;
        }
        else if((jpos - center)*direction < 0)
        {
          ddata.u1 = j;
          ddata.pu = u;
        }
        if(ddata.u1 and ddata.v1)
          break;
      }
    }
    
    vvcomplex::testDivisionOnVertices(c, ddata, T, 0.01);
    
    // apply cell pinching
    tissue::CellPinchingParams params;
    params.cellPinch = _cellPinch;
    params.cellMaxPinch = _cellMaxPinch;
    tissue::cellPinching( c, T, ddata, params );
    
    // TODO: check if the distance between ddata.pv and v1 or v2 is smaller than
    // some eps OR if distance between ddata.pu and u1 or u2 is smaller than eps
    
    
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