#include "ModelUtils.h"

/**
  @file   ModelUtils.cpp
  @brief  Contains namespace for util methods used in the model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

const double EPS = 0.000000001;

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

void findAppropriateJunctionPoint( Point3d &p,
                                   const std::set<junction> &juncs,
                                   const MyTissue& T,
                                   const cell& c,
                                   const bool clockwise )
{
  Point3d temp = p;
  while( findPointInJunctionSet( temp, juncs ) )
  {
    for( auto iter = juncs.begin(); iter != juncs.end(); ++iter )
    {
      if( equalPoints( temp, *iter ) )
      {
        junction newJ;
        if( !clockwise )
          newJ = T.S.prevTo(c, *iter);
        else
          newJ = T.S.nextTo(c, *iter);
        
        temp = newJ->getPos();
        break;
      }
    }
  }
  
  p = temp;
}

// ---------------------------------------------------------------------

bool findPointInJunctionSet( const Point3d &p,
                             const std::set<junction> &juncs )
{
  for( auto iter = juncs.begin(); iter != juncs.end(); ++iter )
  {
    if( equalPoints( p, *iter ) )
      return true;
  }
  
  return false;
}

// ---------------------------------------------------------------------

bool equalPoints( const Point3d &p1, const junction &j2 )
{
  Point3d p2 = j2->getPos();
  return equalPoints( p1, p2 );
}

// ---------------------------------------------------------------------

bool equalPoints( const Point3d &p1, const Point3d &p2 )
{
  if( fabs( p1.i() - p2.i() ) <= EPS &&
      fabs( p1.j() - p2.j() ) <= EPS &&
      fabs( p1.k() - p2.k() ) <= EPS )
    return true;
  else
    return false;
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
    juncs.push_back( j->getPos() );
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

// ---------------------------------------------------------------------

std::set<junction> determineNeedlessJunctions( const cell &c,
                                               const MyTissue& T,
                                               Point3d &center,
                                               const double eps )
{
  // first create linked information of junctions such that
  // later each triple of junctions can be compared to apply a lod
  std::set<junction> junctions;
  std::vector<Point3d> juncs;
  //std::cout << "Orig Pos:" << std::endl;
  forall( const junction& j, T.S.neighbors(c) )
  {
    juncs.push_back( j->getPos() );
    center += juncs.back();
    //std::cout << juncs.back() << std::endl;
  }
  
  center /= juncs.size();
  
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
    
    //std::cout << "Removed Pos:" << std::endl;
    for( std::size_t i=0; i<juncs.size(); i++ )
    {
      if( eraseNodes.find(i) == eraseNodes.end() )
      {
        newJuncs.push_back( juncs.at(i) );
        center += newJuncs.back();
      }
      // else find the corresponding junction pointer and insert
      // it into the set of needless junctions
      else
      {
        forall( const junction& j,T.S.neighbors(c) )
        {
          // if we have found the corresponding junction position
          if( equalPoints( j->getPos(), juncs.at(i) ) )
          {
            //std::cout << jp << std::endl;
            junctions.insert( j );
            break;
          }
        }
      }
    }
    
    center /= newJuncs.size();
  }
  
  return junctions;
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
                                      const double epsLength,
                                      const double epsLOD,
                                      const MyTissue& T )
{
  std::vector<MyTissue::division_data> divisionData;
  
  //std::cout << std::endl;
  
  // we perform a LOD to the current cell such that for three colinear junctions
  // the inner one is not considered because then the method to avoid triangles
  // will not work any more
  Point3d center = Point3d( 0., 0., 0. );
  std::set<junction> juncs = ModelUtils::determineNeedlessJunctions(
    c, T, center, epsLOD );
  
  //std::cout << "orig size: " << T.S.neighbors(c).size() << std::endl;
  //std::cout << "size: " << juncs.size() << std::endl;
  
  double deltaAngle = 3.;
  for( double angle = 0.; angle < 180.; angle += deltaAngle )
  {
    // points for storing cell wall information
    Point3d u1,u2,v1,v2;
    MyTissue::division_data ddata;
    //const Point3d& center = c->center;
    double a = M_PI/180. * angle;
    Point3d direction = Point3d(-sin(a), cos(a), 0);
    
    // for each cell wall
    forall( const junction& j,T.S.neighbors(c) )
    {
      Point3d jpos, jnpos;
      const junction& jn = T.S.nextTo(c, j);
      jpos = j->getPos();
      jnpos = jn->getPos();
      
      Point3d u;
      double s;
      // compute intersection between line [jpos, jnpos] and the plane with
      // position center and normal direction; the result is stored in u and
      // s: Position of the intersection on the line in [0,1].
      if(geometry::planeLineIntersection(u, s, center, direction, jpos, jnpos) and s >= 0 and s <= 1)
      {
        if((jpos - center)*direction > 0)
        {
          ddata.v1 = j;
          ddata.pv = u;
          v1 = jpos;
          v2 = jnpos;
        }
        else if((jpos - center)*direction < 0)
        {
          ddata.u1 = j;
          ddata.pu = u;
          u1 = jpos;
          u2 = jnpos;
        }
        
        if(ddata.u1 and ddata.v1)
          break;
      }
    }
    
    // check if the distance between ddata.pv and v1 or v2 is smaller than
    // some eps OR if distance between ddata.pu and u1 or u2 is smaller than eps
    bool pass = true;
    
    // compute junction lengths depending on the simplified cell structure
    // after applying a level of detail such that our approach in avoiding
    // triangle-shaped cells are guaranteed
    // u1 -> right of pu -> ccw
    // u2 -> left of pu -> cw
    // v1 -> right of pv -> ccw
    // v2 -> left of pv -> cw
    findAppropriateJunctionPoint( u1, juncs, T, c, false );
    findAppropriateJunctionPoint( u2, juncs, T, c, true );
    findAppropriateJunctionPoint( v1, juncs, T, c, false );
    findAppropriateJunctionPoint( v2, juncs, T, c, true );
    
    double uJunctionLength = norm( u2 - u1 );
    double vJunctionLength = norm( v2 - v1 );
    double dist1 = norm( ddata.pu - u1 );
    double dist2 = norm( ddata.pu - u2 );
    double dist3 = norm( ddata.pv - v1 );
    double dist4 = norm( ddata.pv - v2 );
    
    // compute the percentage of length that is set by the user
    double uPercLength = (uJunctionLength*epsLength)/100.;
    double vPercLength = (vJunctionLength*epsLength)/100.;
    
    /*
    std::cout << "uJunctionLength: " << uJunctionLength << std::endl;
    std::cout << "vJunctionLength: " << vJunctionLength << std::endl;
    std::cout << "uPercLength: " << uPercLength << std::endl;
    std::cout << "vPercLength: " << vPercLength << std::endl;
    std::cout << "dist1: " << dist1 << std::endl;
    std::cout << "dist2: " << dist2 << std::endl;
    std::cout << "dist3: " << dist3 << std::endl;
    std::cout << "dist4: " << dist4 << std::endl;
    std::cout << "pu: " << ddata.pu << std::endl;
    std::cout << "pv: " << ddata.pv << std::endl;
    std::cout << "u1: " << u1 << " u2: " << u2 << std::endl;
    std::cout << "v1: " << v1 << " v2: " << v2 << std::endl;
    */
    
    if( dist1 < uPercLength || dist2 < uPercLength ||
        dist3 < vPercLength || dist4 < vPercLength )
      pass = false;
    
    if( pass )
      divisionData.push_back( ddata );
    
    //std::cout << "pass: " << pass << std::endl;
  }
  
  return divisionData;
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

double getNormalDistribution( const double x,
                              const double mu,
                              const double sd )
{
  double a = 1./std::sqrt(2.*M_PI*sd*sd);
  double e = std::exp( (-(x-mu)*(x-mu))/(2.*sd*sd) );
  return a * e;
}

// ---------------------------------------------------------------------

double getSD( const std::vector<double> &vals,
              const double mu )
{
  double sum = 0.;
  for( std::size_t v=0; v < vals.size(); v++ )
    sum += (vals.at(v)-mu)*(vals.at(v)-mu);
  
  sum /= vals.size();
  
  return std::sqrt(sum);
}

// ---------------------------------------------------------------------

void drawLine( const Point3d &pos1, const Point3d &pos2,
               const util::Palette::Color &color )
{
  glLineWidth(2.5); 
  glColor4fv( color.c_data() );
  glBegin( GL_LINES );
  glVertex3f( pos1.i(), pos1.j(), 2. );
  glVertex3f( pos2.i(), pos2.j(), 2. );
  glEnd();
}

// ---------------------------------------------------------------------

void drawControlPoint( const Point3d &pos,
                       const util::Palette::Color &color )
{
  double halfLength = 1.;
  glNormal3f( 0., 0., 1. );
  glPolygonMode( GL_FRONT, GL_FILL );
  glColor4fv( color.c_data() );
  glPushMatrix();
  glBegin( GL_QUADS );
  glVertex3f( pos.i()-halfLength, pos.j()+halfLength, 1. );
  glVertex3f( pos.i()+halfLength, pos.j()+halfLength, 1. );
  glVertex3f( pos.i()+halfLength, pos.j()-halfLength, 1. );
  glVertex3f( pos.i()-halfLength, pos.j()-halfLength, 1. );
  glEnd();
  glPopMatrix();
}

// ---------------------------------------------------------------------

void drawSphere( const Point3d &pos, double r, std::size_t lats,
                 std::size_t longs, const util::Palette::Color &color,
                 GLUquadricObj *quadratic )
{
  glColor4fv( color.c_data() );
  glShadeModel( GL_SMOOTH );
  GLfloat lmodel_ambient[] = { 0.01, 0.01, 0.01, 1.0 };
  glPushMatrix();
  glLightModelfv( GL_LIGHT_MODEL_AMBIENT, lmodel_ambient );
  glTranslatef( pos.i(), pos.j(), 0. );
  gluSphere( quadratic, r, lats, longs );
  glPopMatrix();
} 

 // ---------------------------------------------------------------------

std::size_t getRandomResultOfDistribution( const std::vector<double> &probs )
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<int> d(probs.begin(), probs.end());
  return d(gen);
}

// ---------------------------------------------------------------------

}