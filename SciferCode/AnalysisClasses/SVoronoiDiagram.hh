#ifndef SVORONOIDIAGRAM_HH
#define SVORONOIDIAGRAM_HH

/**
  @file   SVoronoiDiagram.hh
  @brief  Contains class for handling a voronoi diagram of voronoi cells
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SRenderVoronoiCells.hh"

#include <vector>

#include <osg/BoundingBox>

class SVoronoiDiagram
{
public:
  SVoronoiDiagram( const std::string voronoiDirectory,
                   boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                   osg::ref_ptr<SSliderCallback> sliderCallback,
                   boost::shared_ptr<cellTimeVector> cTS );

  osg::Vec3 getCentroid( const int cellId, const int timeStep ) const
  { return _voronoiCells->at(timeStep-1).at(cellId)->getCentroid(); }

  std::set<osg::Vec3> getVertices( const int cellId, const int timeStep ) const
  { return _voronoiCells->at(timeStep-1).at(cellId)->getVertices(); }

  std::set<int> getNeighborIds( const int cellId, const int timeStep ) const
  { return _voronoiCells->at(timeStep-1).at(cellId)->getNeighborIds(); }

  std::set<const SLineageTree*> getNeighbors( const int cellId, const int timeStep ) const
  { return _voronoiCells->at(timeStep-1).at(cellId)->getNeighbors(); }

  osg::ref_ptr<osg::Group> getVoronoi() const
  { return _renderer; }

  void readVoronoiEdges( const std::string &filename,
                         const int timeStep );

  void readBBoxInformation( const std::string &filename,
                            const int timeStep );

  void readNeighborAndVolumeInformation( const std::string &filename,
                                         const int timeStep );

  void storeVoronoiCells( const std::string &fileName );

  void loadVoronoiCells( const std::string &fileName );

  // archiving of all class parameters
  friend class boost::serialization::access;

  template<class Archive>
  void save( Archive &ar, const unsigned int version ) const
  {
    ar << _numTimeSteps;
    ar << _numVornoiCells;

    for( std::size_t t = 0; t < _voronoiCells->size(); t++ )
    {
      for( std::map<int, SVoronoiCell*>::const_iterator mapIter = _voronoiCells->at(t).begin();
           mapIter != _voronoiCells->at(t).end(); ++mapIter )
      {
        mapIter->second->save( ar, version );
      }
    }

  }

  template<class Archive>
  void load( Archive &ar, const unsigned int version )
  {
    ar >> _numTimeSteps;
    ar >> _numVornoiCells;

    _voronoiCells = boost::shared_ptr< std::vector< std::map<int, SVoronoiCell*> > >(
          new std::vector< std::map<int, SVoronoiCell*> > );

    _voronoiCells->resize( _numTimeSteps );

    for( std::size_t t = 0; t < _voronoiCells->size(); t++ )
    {
      for( int c = 0; c < _numVornoiCells.at(t); c++ )
      {
        SVoronoiCell *vC = new SVoronoiCell();

        vC->load( ar, version );

        _voronoiCells->at(t).insert( std::make_pair( vC->getCellId(), vC ) );
      }
    }
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()

private:
  /// complete number of time steps
  int _numTimeSteps;

  /// number of voronoi cells per time step
  std::vector<int> _numVornoiCells;

  /// voronoi cells for each time step and cell id
  VoronoiSet _voronoiCells;

  /// bounding box information for each time step
  std::vector<osg::BoundingBox> _bboxes;

  /// voronoi renderer
  osg::ref_ptr<SRenderVoronoiCells> _renderer;
};

#endif // SVORONOIDIAGRAM_HH
