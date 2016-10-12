#ifndef SVORONOICELL_HH
#define SVORONOICELL_HH

/**
  @file   SVoronoiCell.hh
  @brief  Contains class for handling a voronoi cell
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "dataSet/SLineageTree.hh"

#include <vector>
#include <osg/Group>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/set.hpp>

typedef std::vector< std::pair<osg::Vec3,osg::Vec3> > EdgeVector;

class SVoronoiCell
{
public:
  SVoronoiCell(){}

  SVoronoiCell( const int cellId,
                const int timeStep );

  void setVolume( const double volume )
  { _volume = volume; }

  double volume() const
  { return _volume; }

  void innerCell( const bool status )
  { _innerCell = status; }

  bool innerCell() const
  { return _innerCell; }

  void generateVertices();

  void computeCentroid();

  osg::Vec3 getCentroid() const
  { return _centroid; }

  void setEdges( const EdgeVector &edges );

  EdgeVector getEdgeVector() const
  { return _edges; }

  std::set<osg::Vec3> getVertices() const
  { return _vertices; }

  void setNeighbors( const std::set<const SLineageTree*> &neighbors )
  { _neighbors = neighbors; }

  std::set<const SLineageTree*> getNeighbors() const
  { return _neighbors; }

  void setNeighborIds( const std::set<int> &neighborIds )
  { _neighborIds = neighborIds; }

  std::set<int> getNeighborIds() const
  { return _neighborIds; }

  void setNode( SLineageTree *node )
  { _node = node; }

  SLineageTree* getNode() const
  { return _node; }

  int getCellId() const
  { return _cellId; }

  // archiving of all class parameters
  friend class boost::serialization::access;

  template<class Archive>
  void save( Archive &ar, const unsigned int version ) const
  {
    ar << _cellId;
    ar << _timeStep;
    ar << _volume;
    ar << _innerCell;
    ar << _neighborIds;
  }

  template<class Archive>
  void load( Archive &ar, const unsigned int version )
  {
    ar >> _cellId;
    ar >> _timeStep;
    ar >> _volume;
    ar >> _innerCell;
    ar >> _neighborIds;
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()

private:

  int _cellId;
  int _timeStep;
  double _volume;
  bool _innerCell;
  EdgeVector _edges;
  std::set<osg::Vec3> _vertices;
  osg::Vec3 _centroid;
  SLineageTree *_node;
  std::set<const SLineageTree*> _neighbors;
  std::set<int> _neighborIds;
};

#endif // SVORONOICELL_HH
