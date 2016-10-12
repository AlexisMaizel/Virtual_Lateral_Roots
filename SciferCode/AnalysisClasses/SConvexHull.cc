#include "SConvexHull.hh"

/**
  @file   SConvexHull.cc
  @brief  Contains class for generating a 3D convex hull
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/bounds.hpp>

#include <CGAL/convex_hull_3_to_polyhedron_3.h>

// ---------------------------------------------------------------------

SConvexHull::SConvexHull()
  : SSurfaceGenerator()
{
}

// ---------------------------------------------------------------------

SConvexHull::SConvexHull( const osg::ref_ptr<osg::Vec3Array> vertices,
                          const std::size_t timeStep,
                          const osg::Vec4 &color,
                          const bool render,
                          const bool wireframe )
  : SSurfaceGenerator( vertices, timeStep, wireframe )
{
  this->triangulate();

  // initialize colormap properties
  _colorMap = new SColorMap;
  //_colorMap->setCoolToWarm();
  SColorMap::ColorList &cl = _colorMap->getColorList();
  _colorMap->setColorListNormalized( true );
  cl.clear();
  cl[ 0.0/1.0 ] = SColor( color );
  cl[ 1.0/1.0 ] = SColor( color );
  _colorMap->scaleColors( 0., 1. );

  if( render )
    this->renderSurface();
}

// ---------------------------------------------------------------------

SConvexHull::SConvexHull( const osg::ref_ptr<osg::Vec3Array> vertices,
                          const std::size_t timeStep,
                          const curvatureType cType,
                          std::vector<std::pair<double,double> > &minMaxCurvatures,
                          const std::map<osg::Vec3,std::pair<std::size_t,std::size_t> > &cells,
                          const bool wireframe )
  : SSurfaceGenerator( vertices, timeStep, cType, cells, wireframe )
{
  this->triangulate();
  this->computeDiscreteMeanAndGaussianCurvature();

  // apply logarithmic scaling to all values since the
  // resulting values are strongly varying
  this->applyLogarithmicScaleToCurvature();

  // compute min and max of curvature values
  // first pair is for mean and second for gaussian min/max values
  _minMaxCurvatures.resize( 2, std::make_pair( boost::numeric::bounds<double>::highest(),
                                               boost::numeric::bounds<double>::lowest() ) );

  // determine min and max values of curvatures in order to scale the color map
  this->determineMinMaxCurvature();

  std::size_t type;
  switch( cType )
  {
  case MEANCURVATURE: type = 0; break;
  case GAUSSIANCURVATURE: type = 1; break;
  default: type = 0; break;
  }

  for( std::size_t c = 0; c < 2; c++ )
  {
    if( _minMaxCurvatures.at(c).first <= minMaxCurvatures.at(c).first )
      minMaxCurvatures.at(c).first = _minMaxCurvatures.at(c).first;

    if( _minMaxCurvatures.at(c).second > minMaxCurvatures.at(c).second )
      minMaxCurvatures.at(c).second = _minMaxCurvatures.at(c).second;
  }

  // initialize colormap properties
  _colorMap = new SColorMap;
  //_colorMap->setCoolToWarm();
  SColorMap::ColorList &cl = _colorMap->getColorList();
  _colorMap->setColorListNormalized( true );
  cl.clear();
  cl[ 0.0/4.0 ] = SColor(0,0,1);
  cl[ 1.0/4.0 ] = SColor(0,1,1);
  cl[ 2.0/4.0 ] = SColor(0,1,0);
  cl[ 3.0/4.0 ] = SColor(1,1,0);
  cl[ 4.0/4.0 ] = SColor(1,0,0);
  //_colorMap->setHSV();

  _colorMap->scaleColors( _minMaxCurvatures.at(type).first, _minMaxCurvatures.at(type).second );

  this->renderSurface();
}

// ---------------------------------------------------------------------

void SConvexHull::triangulate()
{
  _triangles->setName( "Convex Hull for t= " +
                       boost::lexical_cast<std::string>( _timeStep ) );

  // if the membrane should be approximated then compute
  // the centre of mass and extend all positions along the line
  // (pos - center) by a fixed chosen length
  if( _surfaceOffset > 0. )
  {
    osg::Vec3 cen( 0., 0., 0. );
    for( std::size_t p = 0; p < _vertices->size(); p++ )
      cen += _vertices->at(p);

    cen /= _vertices->size();

    for( std::size_t p = 0; p < _vertices->size(); p++ )
    {
      osg::Vec3 dir = _vertices->at(p) - cen;
      dir.normalize();
      _vertices->at(p) = _vertices->at(p) + dir * _surfaceOffset;
    }
  }

  std::vector< std::pair<Point, unsigned int> > p;

  unsigned int numPoints = 0;

  for( osg::Vec3Array::iterator iter = _vertices->begin();
       iter != _vertices->end(); ++iter )
  {
    p.push_back( std::make_pair( Point( (*iter)[0], (*iter)[1], (*iter)[2] ),
                                 numPoints ) );
    numPoints++;
  }

  // compute delaunay
  _d.insert( p.begin(), p.end() );
  assert( _d.is_valid() );

  // get the convex hull surface
  CGAL::convex_hull_3_to_polyhedron_3( _d, _poly );
}

// ---------------------------------------------------------------------

const double SConvexHull::computeVolume() const
{
  double volume = 0.;
  for( Delaunay::Finite_cells_iterator it = _d.finite_cells_begin();
       it != _d.finite_cells_end(); it++ )
  {
    Tetrahedron tetr = _d.tetrahedron( it );
    volume += tetr.volume();
  }

  //std::cout << "volume: " << volume << std::endl;
  return volume;
}

// ---------------------------------------------------------------------

void SConvexHull::renderSurface()
{
  for( Polyhedron::Facet_iterator iter = _poly.facets_begin();
       iter != _poly.facets_end(); ++iter )
  {
    // circulator for half edges around the facets
    FHIterator ci = iter->facet_begin();

    osg::ref_ptr<osg::Vec3Array> coords = new osg::Vec3Array;

    // get the three points of the facet
    for( std::size_t j=0;j<3;j++ )
    {
      osg::Vec3 p;

      for( std::size_t i=0;i<3;i++ )
        p[i] = ci->vertex()->point()[i];

      coords->push_back( p );
      ci++;
    }

//    osg::Vec3 normal = ( coords->at(1) - coords->at(0) ) ^ ( coords->at(2) - coords->at(1) );
//    normal.normalize();
//    this->renderSurfaceNormal( (coords->at(0)+coords->at(1)+coords->at(2))/3., normal, _triangles, osg::Vec4( 1., 0., 0., 1. ) );

    this->renderTriangles( coords );
  }
}

// ---------------------------------------------------------------------

void SConvexHull::computeDiscreteMeanAndGaussianCurvature()
{
  // The computation of the discrete mean and gaussian curvatures
  // is based on the following paper from Meyer et al.:
  // Discrete Differential-Geometry Operators for Triangulated 2-Manifolds

  // circulator for half edges around the facets
  FHIterator ci;

  double area[3];
  double angle[3];

  std::map<osg::Vec3,double> areaMap;
  std::map<osg::Vec3,osg::Vec3> contrMap;

  // first traverse over all vertices and initialize start values
  for( Polyhedron::Vertex_iterator iter = _poly.vertices_begin();
       iter != _poly.vertices_end(); ++iter )
  {
    osg::Vec3 p;

    for( std::size_t i=0;i<3;i++ )
      p[i] = iter->point()[i];

    // initial value for areas of a vertex
    std::map<osg::Vec3,double>::iterator mapIter = areaMap.find( p );
    if( mapIter == areaMap.end() )
      areaMap.insert( std::make_pair( p, 0. ) );

    // initial value for contribution of a vertex
    std::map<osg::Vec3,osg::Vec3>::iterator mapIter2 = contrMap.find( p );
    if( mapIter2 == contrMap.end() )
      contrMap.insert( std::make_pair( p, osg::Vec3( 0., 0., 0. ) ) );

    // initial value for gaussian curvature of a vertex
    std::map<osg::Vec3,double>::iterator mapIter3 = _gaussianCurvature.find( p );
    if( mapIter3 == _gaussianCurvature.end() )
      _gaussianCurvature.insert( std::make_pair( p, 2. * M_PI ) );

    // initial value for mean curvature of a vertex
    std::map<osg::Vec3,double>::iterator mapIter4 = _meanCurvature.find( p );
    if( mapIter4 == _meanCurvature.end() )
      _meanCurvature.insert( std::make_pair( p, 0. ) );
  }

  // traverse over all facets for computation
  // of the vornoi areas of obtuse and non-obtuse triangles
  for( Polyhedron::Facet_iterator iter = _poly.facets_begin();
       iter != _poly.facets_end(); ++iter )
  {
    ci = iter->facet_begin();

    osg::Vec3 p[3];

    // get the three points of the facet
    for( std::size_t j=0;j<3;j++ )
    {
      for( std::size_t i=0;i<3;i++ )
        p[j][i] = ci->vertex()->point()[i];

      ci++;
    }

    // compute angles
    angle[0] = this->computeAngle( p[1]-p[0], p[2]-p[0] );
    angle[1] = this->computeAngle( p[0]-p[1], p[2]-p[1] );
    angle[2] = M_PI - (angle[0] + angle[1]);

    // for non-obtuse triangles
    if( (angle[0] < M_PI/2.) && (angle[1] < M_PI/2.) && (angle[2] < M_PI/2.) )
    {
      double e01 = (p[1] - p[0]).length2();
      double e12 = (p[2] - p[1]).length2();
      double e20 = (p[0] - p[2]).length2();

      area[0] = ( e20*( 1./tan(angle[1]) ) + e01*( 1./tan(angle[2]) ) ) / 8.;
      area[1] = ( e01*( 1./tan(angle[2]) ) + e12*( 1./tan(angle[0]) ) ) / 8.;
      area[2] = ( e12*( 1./tan(angle[0]) ) + e20*( 1./tan(angle[1]) ) ) / 8.;

      // store area values for the three vertices
      for( std::size_t i = 0; i < 3; i++ )
        areaMap.at(p[i]) += area[i];
    }
    // for obtuse triangles
    else
    {
      for( std::size_t i = 0; i < 3; i++ )
      {
        if( angle[i] >= M_PI/2. )
        {
          areaMap.at(p[i]) += this->computeArea(*iter) / 2.;
          areaMap.at(p[(i+1)%3]) += this->computeArea(*iter) / 4.;
          areaMap.at(p[(i+2)%3]) += this->computeArea(*iter) / 4.;

          break;
        }
      }
    }
  }

  // second traversal over all facets to compute the gaussian curvature
  for( Polyhedron::Facet_iterator iter = _poly.facets_begin();
       iter != _poly.facets_end(); ++iter )
  {
    ci = iter->facet_begin();

    osg::Vec3 p[3];

    // get the three points of the facet
    for( std::size_t j=0;j<3;j++ )
    {
      for( std::size_t i=0;i<3;i++ )
        p[j][i] = ci->vertex()->point()[i];

      ci++;
    }

    // compute angles
    angle[0] = this->computeAngle( p[1]-p[0], p[2]-p[0] );
    angle[1] = this->computeAngle( p[0]-p[1], p[2]-p[1] );
    angle[2] = M_PI - (angle[0] + angle[1]);

    if( angle[0] == 0. || angle[1] == 0. || angle[2] == 0. )
    //if( angle[0] <= 0.035 || angle[1] <= 0.035 || angle[2] <= 0.035 )
    {
      std::cout << "ignore" << std::endl;
      continue;
    }

    osg::Vec3 e01Dir = p[1] - p[0];
    osg::Vec3 e12Dir = p[2] - p[1];
    osg::Vec3 e20Dir = p[0] - p[2];

    contrMap.at(p[0]) += ( e20Dir * ( 1./tan( angle[1] ) ) - e01Dir * ( 1./tan( angle[2] ) ) ) / 4.;
    contrMap.at(p[1]) += ( e01Dir * ( 1./tan( angle[2] ) ) - e12Dir * ( 1./tan( angle[0] ) ) ) / 4.;
    contrMap.at(p[2]) += ( e12Dir * ( 1./tan( angle[0] ) ) - e20Dir * ( 1./tan( angle[1] ) ) ) / 4.;

    for( std::size_t i = 0; i < 3; i++ )
      _gaussianCurvature.at(p[i]) -= angle[i];

    // TODO: not sure if this is correct
//    {
//      osg::Vec3 e1 = p[1] - p[0];
//      osg::Vec3 e2 = p[2] - p[0];
//      e1.normalize();
//      e2.normalize();
//      double angle = acos( e1 * e2 );
//      _gaussianCurvature.at(p[0]) -= angle;

//      e1 = p[2] - p[1];
//      e2 = p[0] - p[1];
//      e1.normalize();
//      e2.normalize();
//      angle = acos( e1 * e2 );
//      _gaussianCurvature.at(p[1]) -= angle;

//      e1 = p[0] - p[2];
//      e2 = p[1] - p[2];
//      e1.normalize();
//      e2.normalize();
//      angle = acos( e1 * e2 );
//      _gaussianCurvature.at(p[2]) -= angle;
//    }
  }

  // last traversal over all vertices
  for( Polyhedron::Vertex_iterator iter = _poly.vertices_begin();
       iter != _poly.vertices_end(); ++iter )
  {
    osg::Vec3 p;

    for( std::size_t i=0;i<3;i++ )
      p[i] = iter->point()[i];

    if( areaMap.at(p) <= std::numeric_limits<double>::epsilon() )
    {
      _gaussianCurvature.at(p) = 0.;
      _meanCurvature.at(p) = 0.;
    }
    else
    {
      _gaussianCurvature.at(p) /= areaMap.at(p);

      double query = 0.;
      osg::Vec3 n = this->getVertexNormal(p);
      for( std::size_t i = 0; i < 3; i++ )
        query += contrMap.at(p)[i] * n[i];

      if( query > 0. )
        _meanCurvature.at(p) = 1. * ( contrMap.at(p) / areaMap.at(p) ).length();
      else
        _meanCurvature.at(p) = -1. * ( contrMap.at(p) / areaMap.at(p) ).length();
    }
  }
}

// ---------------------------------------------------------------------

const osg::Vec3 SConvexHull::getVertexNormal( const osg::Vec3 &vertex )
{
  // circulator for half edges around the facets
  FHIterator ci;

  //std::cout << "vertex: " << vertex << std::endl;

  // initialize vertex normal
  osg::Vec3 nv( 0., 0., 0. );
  for( Polyhedron::Facet_iterator iter = _poly.facets_begin();
       iter != _poly.facets_end(); ++iter )
  {
    ci = iter->facet_begin();

    osg::Vec3 p[3];

    // get the three points of the facet
    for( std::size_t j=0;j<3;j++ )
    {
      for( std::size_t i=0;i<3;i++ )
        p[j][i] = ci->vertex()->point()[i];

      ci++;
    }

    if( p[0] == vertex ||
        p[1] == vertex ||
        p[2] == vertex )
    {
      osg::Vec3 dir1,dir2;

      // depending on the current vertex, determine
      // the two direction vectors defining the plane
      // and compute the normal pointing outwards
      // Note that if vertex == p[1], the directions
      // are swapped
      if( p[0] == vertex )
      {
        dir1 = p[1] - vertex; dir2 = p[2] - vertex;
      }
      else if( p[1] == vertex )
      {
        dir2 = p[0] - vertex; dir1 = p[2] - vertex;
      }
      else if( p[2] == vertex )
      {
        dir1 = p[0] - vertex; dir2 = p[1] - vertex;
      }

      // compute normal and add it
      osg::Vec3 locNormal = dir1 ^ dir2;
      locNormal.normalize();

      //this->renderSurfaceNormal( (p[0]+p[1]+p[2])/3., locNormal, _tetrahedrons );

      nv += locNormal;
      //std::cout << "locNormal: " << locNormal << std::endl;
    }
  }

  // at last normalize normal
  nv.normalize();

  return nv;
}

// ---------------------------------------------------------------------

const osg::Vec3 SConvexHull::getSurfaceNormal( const osg::Vec3 &vertex,
                                               osg::Vec3 &closestPoint )
{
  // constructs AABB tree and computes internal KD-tree
  // data structure to accelerate distance queries
  // TODO: use this for later release of CGAL (>4.5) because
  // AABB_polyhedron_triangle_primitive is deprectated
  Tree tree( _poly.facets_begin(), _poly.facets_end(), _poly );
  //Tree tree( _poly.facets_begin(), _poly.facets_end() );
  tree.accelerate_distance_queries();

  // query point
  Point query( vertex[0], vertex[1], vertex[2] );

  // compute closest point and primitive id
  Point_and_primitive_id pp = tree.closest_point_and_primitive( query );
  // pp is a pair storing first the closest point and
  // second the primitive id which we use to get the closest triangle
  Polyhedron::Face_handle f = pp.second;

  // store closest point
  closestPoint = osg::Vec3( pp.first.x(), pp.first.y(), pp.first.z() );

  osg::Vec3 p1;
  for( std::size_t i=0;i<3;i++ )
    p1[i] = f->halfedge()->vertex()->point()[i];

  osg::Vec3 p2;
  for( std::size_t i=0;i<3;i++ )
    p2[i] = f->halfedge()->next()->vertex()->point()[i];

  osg::Vec3 p3;
  for( std::size_t i=0;i<3;i++ )
    p3[i] = f->halfedge()->next()->next()->vertex()->point()[i];

  osg::Vec3 dir1,dir2;

  dir1 = p2 - p1;
  dir2 = p3 - p1;

  // compute normal and add it
  osg::Vec3 normal = dir1 ^ dir2;
  normal.normalize();
  return normal;
}

// ---------------------------------------------------------------------

const std::vector<osg::Vec3> SConvexHull::getOuterBoundary( const int viewingType,
                                                            const osg::Matrix &rotMat )
{
  // TODO
  std::vector<osg::Vec3> boundary;
  for( Polyhedron::Facet_iterator iter = _poly.facets_begin();
       iter != _poly.facets_end(); ++iter )
  {
    // circulator for half edges around the facets
    FHIterator ci = iter->facet_begin();

    osg::ref_ptr<osg::Vec3Array> coords = new osg::Vec3Array;

    // get the three points of the facet
    for( std::size_t j=0;j<3;j++ )
    {
      osg::Vec3 p;

      for( std::size_t i=0;i<3;i++ )
        p[i] = ci->vertex()->point()[i];

      boundary.push_back( p );
      ci++;
    }
  }
  return boundary;
}

// ---------------------------------------------------------------------
