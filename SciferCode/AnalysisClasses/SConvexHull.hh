#ifndef SCONVEXHULL_HH
#define SCONVEXHULL_HH

/**
  @file   SConvexHull.hh
  @brief  Contains class for generating a 3D convex hull
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SSurfaceGenerator.hh"

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
// TODO: use this for later release of CGAL (>4.5) because
// AABB_polyhedron_triangle_primitive is deprectated
#include <CGAL/AABB_face_graph_triangle_primitive.h>
//#include <CGAL/AABB_polyhedron_triangle_primitive.h>

typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>                 Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                   Delaunay;
typedef Delaunay::Point                                          Point;
typedef Polyhedron::Facet::Halfedge_around_facet_circulator      FHIterator;
// TODO: use this for later release of CGAL (>4.5) because
// AABB_polyhedron_triangle_primitive is deprectated
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron>     Primitive;
//typedef CGAL::AABB_polyhedron_triangle_primitive<K,Polyhedron>     Primitive;
typedef CGAL::AABB_traits<K, Primitive>                          Traits;
typedef CGAL::AABB_tree<Traits>                                  Tree;
typedef Tree::Point_and_primitive_id                             Point_and_primitive_id;

class SConvexHull : public SSurfaceGenerator
{
public:
  SConvexHull();

  SConvexHull( const osg::ref_ptr<osg::Vec3Array> vertices,
               const std::size_t timeStep,
               const osg::Vec4 &color,
               const bool render,
               const bool wireframe = false );

  SConvexHull( const osg::ref_ptr<osg::Vec3Array> vertices,
               const std::size_t timeStep,
               const curvatureType cType,
               std::vector<std::pair<double,double> > &minMaxCurvatures,
               const std::map<osg::Vec3,std::pair<std::size_t,std::size_t> > &cells,
               const bool wireframe = false );

  void triangulate();

  const double computeVolume() const;

  void renderSurface();

  void computeDiscreteMeanAndGaussianCurvature();

  const osg::Vec3 getVertexNormal( const osg::Vec3 &vertex );

  const osg::Vec3 getSurfaceNormal( const osg::Vec3 &vertex,
                                    osg::Vec3 &closestPoint );

  const std::vector<osg::Vec3> getOuterBoundary( const int viewingType,
                                                 const osg::Matrix &rotMat );

private:

  /// instance of CGAL delaunay
  Delaunay _d;

  /// convex hull as polytype
  Polyhedron _poly;
};

#endif // SCONVEXHULL_HH
