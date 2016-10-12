#include "SRenderVoronoiCells.hh"

/**
  @file   SRenderVoronoiCells.cc
  @brief  Contains class for rendering voronoi cells
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include <fstream>

#include <osg/ShapeDrawable>
#include <osg/PositionAttitudeTransform>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <CGAL/squared_distance_3.h>

const double EPSILON = 1.e-10;

// ---------------------------------------------------------------------

SRenderVoronoiCells::SRenderVoronoiCells( const std::size_t timeStep,
                                          osg::ref_ptr<SSliderCallback> sliderCallback )
  : _numTimeSteps( timeStep )
{
  this->addUpdateCallback( sliderCallback );
}

// ---------------------------------------------------------------------

void SRenderVoronoiCells::drawEdges()
{
  // draw the voronoi edges from the stored class information
  for( std::size_t t = 0; t < (std::size_t)_numTimeSteps; t++ )
  {
    osg::ref_ptr<osg::Group> voronoiGroup = new osg::Group;

    for( std::vector<SEdgeInfo>::const_iterator edge = _voronoiEdges.at(t).begin();
         edge != _voronoiEdges.at(t).end(); ++edge )
    {
      osg::Vec4 color;

      if( edge->innerEdge )
        color = osg::Vec4( 0., 1., 0., 1. );
      else
        color = osg::Vec4( 0.8, 0.8, 0.8, 0.1 );

      std::vector<double> points = edge->points;

      osg::Vec3 p1( points[0], points[1], points[2] );
      osg::Vec3 p2( points[3], points[4], points[5] );

      SSimilarityMeasureGraphics::addCylinderBetweenPoints( p1, p2,
                                                        2., color, voronoiGroup );
    }

    voronoiGroup->getOrCreateStateSet()->setMode( GL_BLEND, osg::StateAttribute::ON );
    voronoiGroup->getOrCreateStateSet()->setRenderingHint( osg::StateSet::TRANSPARENT_BIN );
    _voronoiGroup[t+1] = voronoiGroup;
  }
}

// ---------------------------------------------------------------------

void SRenderVoronoiCells::renderVoronoiEdges( const std::vector<osg::BoundingBox> &bboxes,
                                              const VoronoiSet &voronoiCells )
{
  // render voronoi edges if they do not belong to the boundary
  for( std::size_t t = 0; t < (std::size_t)_numTimeSteps; t++ )
  {
    osg::BoundingBox bbox = bboxes.at(t);

    osg::ref_ptr<osg::Group> voronoiGroup = new osg::Group;

    // for each time step we store a set of all edges in order
    // to draw them only once since different voronoi cells
    // can share the same edge
    std::set< std::pair<osg::Vec3, osg::Vec3> > edgeSet;

    std::vector<SEdgeInfo> edgeInfos;
    bool inner = true;

    int edgeCounter = 0;

    // loop over all cells in this time step
    for( std::map<int, SVoronoiCell*>::iterator cellIter = voronoiCells->at(t).begin();
         cellIter != voronoiCells->at(t).end(); ++cellIter )
    {
      EdgeVector edges = cellIter->second->getEdgeVector();

      // loop over edges of voronoi cell
      for( std::size_t e = 0; e < edges.size(); e++ )
      {
        osg::Vec4 color( 0., 1., 0., 1. );

        osg::Vec3 pos1 = edges.at(e).first;
        osg::Vec3 pos2 = edges.at(e).second;

        // check if line formed by both positions belong to the boundary box
        // if so then do not render this boundary cylinder
        if( this->lineLiesOnBBox( bbox, pos1, pos2 ) )
        {
          cellIter->second->innerCell( false );
          color = osg::Vec4( 0.8, 0.8, 0.8, 0.1 );
          inner = false;
        }
        else
          inner = true;

        // if the edge was already rendererd then skip the
        // redundant one
        if( !edgeSet.insert( std::make_pair( pos1, pos2 ) ).second ||
            !edgeSet.insert( std::make_pair( pos2, pos1 ) ).second )
          continue;

        std::vector<double> points;
        points.resize( 6 );
        for( int i = 0; i < 3; i++ )
        {
          points[i] = pos1[i];
          points[i+3] = pos2[i];
        }

        SEdgeInfo edgeInfo;
        edgeInfo.innerEdge = inner;
        edgeInfo.points = points;

        edgeInfos.push_back( edgeInfo );

        edgeCounter++;

        SSimilarityMeasureGraphics::addCylinderBetweenPoints( pos1, pos2,
                                                          2., color, voronoiGroup );
      }
    }

    voronoiGroup->getOrCreateStateSet()->setMode( GL_BLEND, osg::StateAttribute::ON );
    voronoiGroup->getOrCreateStateSet()->setRenderingHint( osg::StateSet::TRANSPARENT_BIN );

    _voronoiGroup[t+1] = voronoiGroup;
    _voronoiEdges.push_back( edgeInfos );
    _numEdges.push_back( edgeCounter );
  }
}

// ---------------------------------------------------------------------

bool SRenderVoronoiCells::lineLiesOnBBox( const osg::BoundingBox &bbox,
                                          const osg::Vec3 &pos1,
                                          const osg::Vec3 &pos2 )
{
  // set line coordinates
  PointK l1( pos1[0], pos1[1], pos1[2] );
  PointK l2( pos2[0], pos2[1], pos2[2] );

  // generate bounday box corners
  std::vector<PointK> bb;
  // 0
  bb.push_back( PointK( bbox.xMin(), bbox.yMin(), bbox.zMin() ) );
  // 1
  bb.push_back( PointK( bbox.xMax(), bbox.yMin(), bbox.zMin() ) );
  // 2
  bb.push_back( PointK( bbox.xMin(), bbox.yMax(), bbox.zMin() ) );
  // 3
  bb.push_back( PointK( bbox.xMax(), bbox.yMax(), bbox.zMin() ) );
  // 4
  bb.push_back( PointK( bbox.xMin(), bbox.yMin(), bbox.zMax() ) );
  // 5
  bb.push_back( PointK( bbox.xMax(), bbox.yMin(), bbox.zMax() ) );
  // 6
  bb.push_back( PointK( bbox.xMin(), bbox.yMax(), bbox.zMax() ) );
  // 7
  bb.push_back( PointK( bbox.xMax(), bbox.yMax(), bbox.zMax() ) );

  // generate all planes of the boundary box
  std::vector<PlaneK> pl;
  pl.push_back( PlaneK( bb[0], bb[1], bb[3] ) );
  pl.push_back( PlaneK( bb[5], bb[4], bb[6] ) );
  pl.push_back( PlaneK( bb[4], bb[0], bb[2] ) );
  pl.push_back( PlaneK( bb[1], bb[5], bb[7] ) );
  pl.push_back( PlaneK( bb[0], bb[4], bb[5] ) );
  pl.push_back( PlaneK( bb[3], bb[6], bb[7] ) );

  // check for both positions if their distance to one of the plane
  // is zero
  // NOTE that we have chosen this kind of query since the 'has_on'
  // method of the plane class of CGAL does not work correctly; which
  // is maybe because of rounding errors
  for( int i = 0; i < 6; i++ )
  {
    //double dist1 = CGAL::squared_distance( pl.at(i), line );
    double dist1 = CGAL::squared_distance( pl.at(i), l1 );
    double dist2 = CGAL::squared_distance( pl.at(i), l2 );

    bool liesOnPlane;
    // if both positions are smaller than some epsilon
    // then this edge is a boundary edge
    if( dist1 <= EPSILON && dist2 <= EPSILON )
      liesOnPlane = true;
    else
      liesOnPlane = false;

    if( liesOnPlane )
      return liesOnPlane;
  }

  return false;
}

// ---------------------------------------------------------------------

void SRenderVoronoiCells::update( const int t )
{
  std::map<int, osg::ref_ptr<osg::Group> >::const_iterator iter =
      _voronoiGroup.find( t );
  if( iter != _voronoiGroup.end() )
  {
    this->removeChild( 0, 1 );
    this->addChild( iter->second );
    this->setNodeMask( ~0 );
  }
  else
    this->setNodeMask( 0 );
}

// ---------------------------------------------------------------------

void SRenderVoronoiCells::storeVoronoiEdges( const std::string &fileName )
{
  std::ofstream outfile( fileName.c_str(), std::ios::out | std::ios::binary );

  // save data to archive
  boost::archive::binary_oarchive oa( outfile );
  // write class instance to archive
  oa << *this;
}

// ---------------------------------------------------------------------

void SRenderVoronoiCells::loadVoronoiEdges( const std::string &fileName )
{
  std::ifstream infile( fileName.c_str(), std::ios::in | std::ios::binary );

  if( (infile.rdstate() & std::ifstream::failbit ) != 0 )
  {
    std::cerr << "Error opening " << fileName << std::endl;
    std::cerr << strerror( errno ) << std::endl;
    return;
  }

  boost::archive::binary_iarchive ia( infile );
  ia >> *this;
}

// ---------------------------------------------------------------------
