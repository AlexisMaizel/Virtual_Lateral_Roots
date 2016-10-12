#ifndef SALPHASHAPE_HH
#define SALPHASHAPE_HH

/**
  @file   SAlphaShape.hh
  @brief  Contains class for generating a 3D alpha shape
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SSurfaceGenerator.hh"

#include <boost/optional.hpp>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::Alpha_shape_vertex_base_3<K>                    AlphaVb;
typedef CGAL::Alpha_shape_cell_base_3<K>                      AlphaFb;
typedef CGAL::Triangulation_data_structure_3<AlphaVb,AlphaFb> AlphaTds;
typedef CGAL::Delaunay_triangulation_3<K, AlphaTds>           AlphaDelaunay;
typedef CGAL::Alpha_shape_3<AlphaDelaunay>                    Alpha_Shape_3;
typedef Alpha_Shape_3::Alpha_iterator                         Alpha_iterator;
typedef K::Point_3                                            PointK;
typedef K::Segment_3                                          SegmentK;
typedef K::Vector_3                                           VectorK;
typedef std::list<Triangle>::iterator                         TriIterator;
typedef CGAL::AABB_triangle_primitive<K,TriIterator>          AlphaPrimitive;
typedef CGAL::AABB_traits<K, AlphaPrimitive>                  AlphaTraits;
typedef CGAL::AABB_tree<AlphaTraits>                          AlphaTree;
typedef AlphaTree::Point_and_primitive_id                     AlphaPoint_and_primitive_id;

// typedefs for 2D alpha shapes
typedef K::FT FT;
typedef K::Point_2  Point_2;
typedef K::Segment_2  Segment_2;
typedef CGAL::Alpha_shape_vertex_base_2<K> Vb_2;
typedef CGAL::Alpha_shape_face_base_2<K>  Fb_2;
typedef CGAL::Triangulation_data_structure_2<Vb_2,Fb_2> Tds_2;
typedef CGAL::Delaunay_triangulation_2<K,Tds_2> Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;

class SAlphaShape : public SSurfaceGenerator
{
public:
  SAlphaShape();

  SAlphaShape( const osg::ref_ptr<osg::Vec3Array> vertices,
               const std::size_t timeStep,
               const osg::Vec4 &color,
               const bool render,
               const bool wireframe = false,
               const double surfaceOffset = 0. );

  SAlphaShape( const osg::ref_ptr<osg::Vec3Array> vertices,
               const std::size_t timeStep,
               const curvatureType cType,
               std::vector<std::pair<double,double> > &minMaxCurvatures,
               const std::map<osg::Vec3,std::pair<std::size_t,std::size_t> > &cells,
               const bool wireframe = false,
               const double surfaceOffset = 0. );

  void storeTriangles();

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
  Alpha_Shape_3 _as;

  std::list<Triangle> _triangleList;

  bool _renderNormals;
};

#endif // SALPHASHAPE_HH
