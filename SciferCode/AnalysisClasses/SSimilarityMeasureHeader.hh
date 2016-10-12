#ifndef SSIMILARITYMEASUREHEADER_HH
#define SSIMILARITYMEASUREHEADER_HH

#include <dataSet/SAdvancedLineageTreeCollection.hh>

#include <boost/signals2.hpp>

#include <osg/Node>

#include <set>

struct compareSet
{
  bool operator()( const SLineageTree *t1, const SLineageTree *t2 ) const
  {
    if( t1->cellId != t2->cellId )
      return ( t1->cellId < t2->cellId );
    else
      return ( t1->timeStep < t2->timeStep );
  }
};

struct comparePairSet
{
  bool operator()( const std::pair<int,int> &p1,
                   const std::pair<int,int> &p2 ) const
  {
    if( p1.first != p2.first )
      return (p1.first < p2.first);
    else
      return ( p1.second < p2.second );
  }
};

struct SphereC
{
  double theta;
  double phi;
  double r;
};

namespace NodeType
{
  enum NodeType{ NODE2D, NODE3D };
}

namespace DivisionType
{
  // order: 0, 1, 2 when exporting this info
  enum DivisionType{ ANTICLINAL, PERICLINAL, RADIAL };
}

typedef std::set<const SLineageTree*, compareSet> cellSet;
// alternative cellSet with unique identifictation by id and timestep
typedef std::set<std::pair<int,int>, comparePairSet> cellPairSet;
typedef std::set<unsigned int> cellHistory;
typedef std::vector< std::vector<cellSet> > cellLayerVector;
typedef std::vector<cellPairSet> cellPairVector;
typedef std::vector<cellSet> cellTimeVector;
typedef std::map<const SLineageTree*, cellHistory> cellHistoryMap;
typedef boost::signals2::signal< void (osg::Node*) > signalHighlight;
typedef signalHighlight::slot_type slotHighlight;
typedef boost::signals2::signal< void (osg::Node*) > signalLayerAssignment;
typedef signalLayerAssignment::slot_type slotLayerAssignment;
typedef std::vector< std::map<int, DivisionType::DivisionType> > divisionAssignment;


#endif // SSIMILARITYMEASUREHEADER_HH
