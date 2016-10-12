#include "SAlphaShape.hh"

/**
  @file   SAlphaShape.cc
  @brief  Contains class for generating a 3D alpha shape
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/bounds.hpp>

// ---------------------------------------------------------------------

SAlphaShape::SAlphaShape()
  : SSurfaceGenerator()
{
}

// ---------------------------------------------------------------------

SAlphaShape::SAlphaShape( const osg::ref_ptr<osg::Vec3Array> vertices,
                          const std::size_t timeStep,
                          const osg::Vec4 &color,
                          const bool render,
                          const bool wireframe,
                          const double surfaceOffset )
  : SSurfaceGenerator( vertices, timeStep, wireframe, surfaceOffset ),
    _renderNormals( false )
{
  this->triangulate();

  // initialize single color in colormap
  _colorMap = new SColorMap;
  SColorMap::ColorList &cl = _colorMap->getColorList();
  _colorMap->setColorListNormalized( true );
  cl.clear();
  cl[ 0.0/1.0 ] = SColor( color );
  cl[ 1.0/1.0 ] = SColor( color );
  _colorMap->scaleColors( 0., 1. );

  // store the triangles of the alpha shape for a later processing
  this->storeTriangles();

  if( render )
    this->renderSurface();

  if( _renderNormals )
  {
    for( std::size_t v=0;v<vertices->size();v++ )
      this->getVertexNormal( vertices->at(v) );
  }
}

// ---------------------------------------------------------------------

SAlphaShape::SAlphaShape( const osg::ref_ptr<osg::Vec3Array> vertices,
                          const std::size_t timeStep,
                          const curvatureType cType,
                          std::vector<std::pair<double,double> > &minMaxCurvatures,
                          const std::map<osg::Vec3,std::pair<std::size_t,std::size_t> > &cells,
                          const bool wireframe,
                          const double surfaceOffset )
  : SSurfaceGenerator( vertices, timeStep, cType, cells, wireframe, surfaceOffset )
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
  case MEANCURVATURE:
  case NOCURVATURE:
  default:
    type = 0; break;
  case GAUSSIANCURVATURE:
    type = 1; break;
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

  // store the triangles of the alpha shape for a later processing
  this->storeTriangles();

  this->renderSurface();
}

// ---------------------------------------------------------------------

void SAlphaShape::storeTriangles()
{
  for( Alpha_Shape_3::Alpha_shape_facets_iterator iter = _as.Alpha_shape_facets_begin();
       iter != _as.Alpha_shape_facets_end(); ++iter )
  {
    Alpha_Shape_3::Cell_handle cell = iter->first;
    int vtx_idx = iter->second;

    if( _as.classify( *iter ) == Alpha_Shape_3::REGULAR &&
        _as.classify( cell ) == Alpha_Shape_3::INTERIOR )
    {
      cell = cell->neighbor( vtx_idx );
      vtx_idx = cell->index( iter->first );
    }

    // compute normal vector of surface
    int idx0 = Alpha_Shape_3::Triangulation_utils_3::vertex_triple_index( vtx_idx, 0 );
    int idx1 = Alpha_Shape_3::Triangulation_utils_3::vertex_triple_index( vtx_idx, 1 );
    int idx2 = Alpha_Shape_3::Triangulation_utils_3::vertex_triple_index( vtx_idx, 2 );

    PointK &vtx0 = cell->vertex(idx0)->point();
    PointK &vtx1 = cell->vertex(idx1)->point();
    PointK &vtx2 = cell->vertex(idx2)->point();

    if( _renderNormals )
    {
      osg::Vec3 p0( vtx0[0], vtx0[1], vtx0[2] );
      osg::Vec3 p1( vtx1[0], vtx1[1], vtx1[2] );
      osg::Vec3 p2( vtx2[0], vtx2[1], vtx2[2] );
      osg::Vec3 normal = ( p1 - p0 ) ^ ( p2 - p1 );
      normal.normalize();
      this->renderSurfaceNormal( (p0+p1+p2)/3., normal, _triangles, osg::Vec4( 0., 1., 0., 1. ) );
    }

    // store each triangle in the list
    _triangleList.push_back( Triangle( vtx0, vtx1, vtx2 ) );
  }
}

// ---------------------------------------------------------------------

void SAlphaShape::triangulate()
{
  _triangles->setName( "Alpha Shape for t= " +
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

  std::list<PointK> lp;

  unsigned int numPoints = 0;

  for( osg::Vec3Array::iterator iter = _vertices->begin();
       iter != _vertices->end(); ++iter )
  {
    lp.push_back( PointK( (*iter)[0], (*iter)[1], (*iter)[2] ) );
    numPoints++;
  }

  // compute alpha shape
  _as.make_alpha_shape( lp.begin(), lp.end() );

  // find optimal alpha value
  Alpha_iterator opt = _as.find_optimal_alpha( 1 );

  _as.set_alpha( *opt );
  // this mode decides if singular simplexes should be considered or not
  // GENERAL -> do consider singular simplexes
  // REGULARIZED -> only single simplex
  // In fact in the case of the primordium there are no singular simplexes
  // but CGAL determines the alpha value too small such that single faces
  // aka triangles are not drawn and the surface will get holes. Thus,
  // the GENERAL mode is activated to avoid these holes.
  _as.set_mode( Alpha_Shape_3::GENERAL );
}

// ---------------------------------------------------------------------

const double SAlphaShape::computeVolume() const
{
  // alpha shape volume
  double alpha_volume = 0.;
  for( Alpha_Shape_3::Finite_cells_iterator it = _as.finite_cells_begin();
       it != _as.finite_cells_end(); it++ )
  {
    if( _as.classify(it) == Alpha_Shape_3::INTERIOR )
    {
      Tetrahedron tetr = _as.tetrahedron( it );
      alpha_volume += tetr.volume();
    }
  }

  //std::cout << "alpha volume: " << alpha_volume << std::endl;

  return alpha_volume;
}

// ---------------------------------------------------------------------

void SAlphaShape::renderSurface()
{
  for( Alpha_Shape_3::Alpha_shape_facets_iterator iter = _as.Alpha_shape_facets_begin();
       iter != _as.Alpha_shape_facets_end(); ++iter )
  {
    Triangle tr = _as.triangle( *iter );

    osg::ref_ptr<osg::Vec3Array> coords = new osg::Vec3Array;

    // get the three points of the triangle
    for( int i=0;i<3;i++ )
    {
      PointK p = tr.vertex(i);
      coords->push_back( osg::Vec3( p.x(), p.y(), p.z() ) );
    }
    this->renderTriangles( coords );
  }
}

// ---------------------------------------------------------------------

void SAlphaShape::computeDiscreteMeanAndGaussianCurvature()
{
  // The computation of the discrete mean and gaussian curvatures
  // is based on the following paper from Meyer et al.:
  // Discrete Differential-Geometry Operators for Triangulated 2-Manifolds

  double area[3];
  double angle[3];

  std::map<osg::Vec3,double> areaMap;
  std::map<osg::Vec3,osg::Vec3> contrMap;

  // first traverse over all facets and initialize start values of unique vertices
  for( Alpha_Shape_3::Alpha_shape_facets_iterator iter = _as.Alpha_shape_facets_begin();
       iter != _as.Alpha_shape_facets_end(); ++iter )
  {
    Triangle tr = _as.triangle( *iter );

    osg::Vec3 p[3];

    // get the three points of the triangle
    for( int i=0;i<3;i++ )
    {
      PointK point = tr.vertex(i);
      p[i] = osg::Vec3( point.x(), point.y(), point.z() );

      // initial value for areas of a vertex
      std::map<osg::Vec3,double>::iterator mapIter = areaMap.find( p[i] );
      if( mapIter == areaMap.end() )
        areaMap.insert( std::make_pair( p[i], 0. ) );

      // initial value for contribution of a vertex
      std::map<osg::Vec3,osg::Vec3>::iterator mapIter2 = contrMap.find( p[i] );
      if( mapIter2 == contrMap.end() )
        contrMap.insert( std::make_pair( p[i], osg::Vec3( 0., 0., 0. ) ) );

      // initial value for gaussian curvature of a vertex
      std::map<osg::Vec3,double>::iterator mapIter3 = _gaussianCurvature.find( p[i] );
      if( mapIter3 == _gaussianCurvature.end() )
        _gaussianCurvature.insert( std::make_pair( p[i], 2. * M_PI ) );

      // initial value for mean curvature of a vertex
      std::map<osg::Vec3,double>::iterator mapIter4 = _meanCurvature.find( p[i] );
      if( mapIter4 == _meanCurvature.end() )
        _meanCurvature.insert( std::make_pair( p[i], 0. ) );
    }
  }

  // traverse over all facets for computation
  // of the vornoi areas of obtuse and non-obtuse triangles
  for( Alpha_Shape_3::Alpha_shape_facets_iterator iter = _as.Alpha_shape_facets_begin();
       iter != _as.Alpha_shape_facets_end(); ++iter )
  {
    osg::Vec3 p[3];
    Triangle tr = _as.triangle( *iter );

    // get the three points of the triangle
    for( int i=0;i<3;i++ )
    {
      PointK po = tr.vertex(i);
      p[i] = osg::Vec3( po.x(), po.y(), po.z() );
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
          areaMap.at(p[i]) += this->computeArea(p) / 2.;
          areaMap.at(p[(i+1)%3]) += this->computeArea(p) / 4.;
          areaMap.at(p[(i+2)%3]) += this->computeArea(p) / 4.;

          break;
        }
      }
    }
  }

  // second traversal over all facets to compute the gaussian curvature
  for( Alpha_Shape_3::Alpha_shape_facets_iterator iter = _as.Alpha_shape_facets_begin();
       iter != _as.Alpha_shape_facets_end(); ++iter )
  {
    osg::Vec3 p[3];
    Triangle tr = _as.triangle( *iter );

    // get the three points of the triangle
    for( int i=0;i<3;i++ )
    {
      PointK po = tr.vertex(i);
      p[i] = osg::Vec3( po.x(), po.y(), po.z() );
    }

    // compute angles
    angle[0] = this->computeAngle( p[1]-p[0], p[2]-p[0] );
    angle[1] = this->computeAngle( p[0]-p[1], p[2]-p[1] );
    angle[2] = M_PI - (angle[0] + angle[1]);

    if( angle[0] == 0. || angle[1] == 0. || angle[2] == 0. )
    //if( angle[0] <= 0.035 || angle[1] <= 0.035 || angle[2] <= 0.035 )
      continue;

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

  // last traversal over all vertices since all start types have the same size
  // and the same content of key values which are the vertices, we perform
  // loop over the areaMap
  for( std::map<osg::Vec3,double>::const_iterator iter = areaMap.begin();
       iter != areaMap.end(); ++iter )
  {
    osg::Vec3 p = iter->first;

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

const osg::Vec3 SAlphaShape::getVertexNormal( const osg::Vec3 &vertex )
{
  // initialize vertex normal
  osg::Vec3 nv( 0., 0., 0. );

  for( std::list<Triangle>::const_iterator tri = _triangleList.begin();
       tri != _triangleList.end(); ++tri )
  {
    PointK pk0 = tri->vertex(0);
    PointK pk1 = tri->vertex(1);
    PointK pk2 = tri->vertex(2);
    osg::Vec3 p0( pk0[0], pk0[1], pk0[2] );
    osg::Vec3 p1( pk1[0], pk1[1], pk1[2] );
    osg::Vec3 p2( pk2[0], pk2[1], pk2[2] );

    if( p0 == vertex ||
        p1 == vertex ||
        p2 == vertex )
    {
      osg::Vec3 dir1,dir2;

      // depending on the current vertex, determine
      // the two direction vectors defining the plane
      // and compute the normal pointing outwards
      // Note that if vertex == p[1], the directions
      // are swapped
      if( p0 == vertex )
      {
        dir1 = p1 - vertex; dir2 = p2 - vertex;
      }
      else if( p1 == vertex )
      {
        dir2 = p0 - vertex; dir1 = p2 - vertex;
      }
      else if( p2 == vertex )
      {
        dir1 = p0 - vertex; dir2 = p1 - vertex;
      }

      // compute normal and add it
      osg::Vec3 locNormal = dir1 ^ dir2;
      locNormal.normalize();
      nv += locNormal;
    }
  }

  // at last normalize normal
  nv.normalize();

  if( _renderNormals )
    this->renderSurfaceNormal( vertex, nv, _triangles, osg::Vec4( 0., 0., 1., 1. ) );

  return nv;
}

// ---------------------------------------------------------------------

const osg::Vec3 SAlphaShape::getSurfaceNormal( const osg::Vec3 &vertex,
                                               osg::Vec3 &closestPoint )
{
  AlphaTree tree( _triangleList.begin(), _triangleList.end() );
  PointK point_query( vertex[0], vertex[1], vertex[2] );
  // compute closest point and primitive id
  AlphaPoint_and_primitive_id pp = tree.closest_point_and_primitive( point_query );

  // get closest triangle
  Triangle tri = *(pp.second);

  // store closest point
  closestPoint = osg::Vec3( pp.first.x(), pp.first.y(), pp.first.z() );

  // compute normal vector
  PointK pk0 = tri.vertex(0);
  PointK pk1 = tri.vertex(1);
  PointK pk2 = tri.vertex(2);
  osg::Vec3 p0( pk0[0], pk0[1], pk0[2] );
  osg::Vec3 p1( pk1[0], pk1[1], pk1[2] );
  osg::Vec3 p2( pk2[0], pk2[1], pk2[2] );
  osg::Vec3 normal = ( p1 - p0 ) ^ ( p2 - p1 );
  normal.normalize();

  return normal;
}

// ---------------------------------------------------------------------

const std::vector<osg::Vec3> SAlphaShape::getOuterBoundary( const int viewingType,
                                                            const osg::Matrix &rotMat )
{
  std::list<Point_2> points;
  for( std::size_t v = 0; v < _vertices->size(); v++ )
  {
    osg::Vec3 pos = _vertices->at(v);
    pos = pos * rotMat;

    switch( viewingType )
    {
      case 0: points.push_back( Point_2( pos[0], pos[1] ) ); break;
      case 1: points.push_back( Point_2( pos[0], pos[2] ) ); break;
      case 2: points.push_back( Point_2( pos[1], pos[2] ) ); break;
    }
  }

  Alpha_shape_2 as( points.begin(), points.end(), FT(10000),
                    Alpha_shape_2::REGULARIZED );

  // find optimal alpha value
  Alpha_shape_2::Alpha_iterator opt = as.find_optimal_alpha( 1 );
  as.set_alpha( *opt * 4. );

  std::vector<osg::Vec3> boundary;
  for( Alpha_shape_edges_iterator it = as.alpha_shape_edges_begin();
       it != as.alpha_shape_edges_end(); ++it )
  {
    if( as.classify( *it ) == Alpha_shape_2::REGULAR )
    {
      Segment_2 seg = as.segment(*it);
      switch( viewingType )
      {
      case 0:
        boundary.push_back( osg::Vec3( seg.source()[0], seg.source()[1], 0. ) );
        boundary.push_back( osg::Vec3( seg.target()[0], seg.target()[1], 0. ) );
        break;
      case 1:
        boundary.push_back( osg::Vec3( seg.source()[0], 0., seg.source()[1] ) );
        boundary.push_back( osg::Vec3( seg.target()[0], 0., seg.target()[1] ) );
        break;
      case 2:
        boundary.push_back( osg::Vec3( 0., seg.source()[0], seg.source()[1] ) );
        boundary.push_back( osg::Vec3( 0., seg.target()[0], seg.target()[1] ) );
        break;
      }
    }
  }

  return boundary;
}

// ---------------------------------------------------------------------
