#ifndef SRENDERVORONOICELLS_HH
#define SRENDERVORONOICELLS_HH

/**
  @file   SRenderVoronoiCells.hh
  @brief  Contains class for rendering voronoi cells
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SVoronoiCell.hh"

#include "SSimilarityMeasureGraphics.hh"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

// loop over all time steps and over all vornoi edges
// with start and end point
typedef boost::shared_ptr< std::vector< std::map<int, SVoronoiCell*> > > VoronoiSet;

struct SEdgeInfo
{
public:

  SEdgeInfo()
    : innerEdge( true )
  {}

  // archiving of all class parameters
  friend class boost::serialization::access;

  template<class Archive>
  void save( Archive &ar, const unsigned int version ) const
  {
    ar << points;
    ar << innerEdge;
  }

  template<class Archive>
  void load( Archive &ar, const unsigned int version )
  {
    ar >> points;
    ar >> innerEdge;
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()

  std::vector<double> points;
  bool innerEdge;
};

typedef std::vector< std::vector<SEdgeInfo> > EdgeInfoVector;

class SRenderVoronoiCells : public osg::Group
{
public:

  SRenderVoronoiCells( const std::size_t timeStep,
                       osg::ref_ptr<SSliderCallback> sliderCallback );

  void renderVoronoiEdges( const std::vector<osg::BoundingBox> &bboxes,
                           const VoronoiSet &voronoiCells );

  bool lineLiesOnBBox( const osg::BoundingBox &bbox,
                       const osg::Vec3 &pos1,
                       const osg::Vec3 &pos2 );

  void update( const int t );

  void drawEdges();

  void storeVoronoiEdges( const std::string &fileName );

  void loadVoronoiEdges( const std::string &fileName );

  // archiving of all class parameters
  friend class boost::serialization::access;

  template<class Archive>
  void save( Archive &ar, const unsigned int version ) const
  {
    ar << _numTimeSteps;
    ar << _numEdges;

    for( std::size_t t = 0; t < _voronoiEdges.size(); t++ )
    {
      for( std::size_t e = 0; e < _voronoiEdges.at(t).size(); e++ )
        _voronoiEdges.at(t).at(e).save( ar, version );;
    }
  }

  template<class Archive>
  void load( Archive &ar, const unsigned int version )
  {
    ar >> _numTimeSteps;
    ar >> _numEdges;

    _voronoiEdges.resize( _numTimeSteps );

    for( std::size_t t = 0; t < _voronoiEdges.size(); t++ )
    {
      _voronoiEdges.at(t).resize( _numEdges.at(t) );

      for( std::size_t e = 0; e < _voronoiEdges.at(t).size(); e++ )
        _voronoiEdges.at(t).at(e).load( ar, version );
    }
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()

private:

  /// voronoi osg group for each time step switched by the time slider
  std::map<int, osg::ref_ptr<osg::Group> > _voronoiGroup;

  EdgeInfoVector _voronoiEdges;

  /// complete number of time steps
  /// TODO: should be std::size_t but not done at the moment due to boost serialization
  int _numTimeSteps;

  /// complete number of edges per time step
  std::vector<int> _numEdges;
};

#endif // SRENDERVORONOICELLS_HH
