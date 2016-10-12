#include "SVoronoiCell.hh"

/**
  @file   SVoronoiCell.cc
  @brief  Contains class for handling a voronoi cell
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

// ---------------------------------------------------------------------

SVoronoiCell::SVoronoiCell( const int cellId,
                            const int timeStep )
  : _cellId( cellId ),
    _timeStep( timeStep ),
    _volume( 0. ),
    _innerCell( true )
{
}

// ---------------------------------------------------------------------

void SVoronoiCell::setEdges( const EdgeVector &edges )
{
  _edges = edges;

  // after setting the edges, the vertices are extracted
  this->generateVertices();

  // after setting the vertices, compute the centroid
  this->computeCentroid();
}

// ---------------------------------------------------------------------

void SVoronoiCell::generateVertices()
{
  for( EdgeVector::iterator edgeIter = _edges.begin();
       edgeIter != _edges.end(); ++edgeIter )
  {
    // insert unique vertices
    _vertices.insert( edgeIter->first );
    _vertices.insert( edgeIter->second );
  }
}

// ---------------------------------------------------------------------

void SVoronoiCell::computeCentroid()
{
  osg::Vec3 sum( 0., 0., 0. );

  for( std::set<osg::Vec3>::const_iterator setIter = _vertices.begin();
       setIter != _vertices.end(); ++setIter )
  {
    sum += *setIter;
  }

  sum /= _vertices.size();
  _centroid = sum;
}

// ---------------------------------------------------------------------
