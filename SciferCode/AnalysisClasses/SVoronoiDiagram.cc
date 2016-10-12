#include "SVoronoiDiagram.hh"

/**
  @file   SVoronoiDiagram.cc
  @brief  Contains class for handling a voronoi diagram of voronoi cells
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "kernel/SKernel.hh"

#include <fstream>

#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>

// ---------------------------------------------------------------------

SVoronoiDiagram::SVoronoiDiagram( const std::string voronoiDirectory,
                                  boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                                  osg::ref_ptr<SSliderCallback> sliderCallback,
                                  boost::shared_ptr<cellTimeVector> cTS )
{
  _numTimeSteps = lineages->getOrComputeNumberOfTimesteps();

  _voronoiCells = boost::shared_ptr< std::vector< std::map<int, SVoronoiCell*> > >(
        new std::vector< std::map<int, SVoronoiCell*> > );
  _voronoiCells->resize( _numTimeSteps );
  _bboxes.resize( _numTimeSteps );

  bool loadClassData = true;

  /* read voronoi information from external files */

  // renderer of the edges of the voronoi diagram
  _renderer = new SRenderVoronoiCells( _numTimeSteps, sliderCallback );

  SKernel::get()->setProgress( 0.25, "Reading voronoi info" );

  if( !loadClassData )
  {
    for( int i = 0; i < _numTimeSteps; i++ )
    {
      // read edges
      std::string filename = voronoiDirectory + "/voronoiCells/vT" +
          boost::lexical_cast<std::string>(i) + ".txt_v.pov";
      this->readVoronoiEdges( filename, i );

      // read bbox information
      std::string filenameBbox = voronoiDirectory + "/bbox/vT" +
          boost::lexical_cast<std::string>(i) + ".txt";
      this->readBBoxInformation( filenameBbox, i );

      // read volume information
      std::string filenameVolume = voronoiDirectory + "/voronoiCells/vT" +
          boost::lexical_cast<std::string>(i) + ".txt.vol";
      this->readNeighborAndVolumeInformation( filenameVolume, i );
    }

    _renderer->renderVoronoiEdges( _bboxes, _voronoiCells );
    _renderer->storeVoronoiEdges( voronoiDirectory + "/edgeInfo.bin" );

    this->storeVoronoiCells( voronoiDirectory + "/voronoiCellInfo.bin" );
  }
  else
  {
    _renderer->loadVoronoiEdges( voronoiDirectory + "/edgeInfo.bin" );
    _renderer->drawEdges();


    this->loadVoronoiCells( voronoiDirectory + "/voronoiCellInfo.bin" );
  }

  _renderer->setName( "Voronoi Diagram" );

  SKernel::get()->setProgress( 0.75, "Generating neighborhood info" );

  // loop over all cells in order to set the corresponding node pointer
  // in the voronoi cells
  for( SAdvancedLineageTreeCollection::const_iterator iter = lineages->begin();
       iter != lineages->end(); ++iter )
  {
    SLineageTree *tree = iter->second;

    // set for all trees the layer information which is later used
    // in the coloring step
    for( SLineageTree::iterator treeIter = tree->begin();
         treeIter != tree->end(); ++treeIter )
    {
      int timeStep = treeIter->timeStep-1;
      int cellId = treeIter->cellId;

      _voronoiCells->at(timeStep).at(cellId)->setNode( *treeIter );

      // for each cell and its set of neighbor ids generate a corresponding
      // set of neighbor SLineageTree node pointers which makes the later
      // access simpler
      std::set<int> neighborIds =
          _voronoiCells->at(timeStep).at(cellId)->getNeighborIds();

      // new set for node pointers
      std::set<const SLineageTree*> neighbors;

      // for each neighbor id
      for( std::set<int>::iterator setIter = neighborIds.begin();
           setIter != neighborIds.end(); ++setIter )
      {
        // for each node at the current time step
        for( cellSet::const_iterator iter = cTS->at(timeStep).begin();
             iter != cTS->at(timeStep).end(); ++iter )
        {
          // find the node with the same cell id
          if( (*iter)->cellId == *setIter )
          {
            neighbors.insert( *iter );
            break;
          }
        }
      }

      _voronoiCells->at(timeStep).at(cellId)->setNeighbors( neighbors );
    }
  }

  SKernel::get()->setProgress( 1., "Done" );
}

// ---------------------------------------------------------------------

void SVoronoiDiagram::readVoronoiEdges( const std::string &filename,
                                        const int timeStep )
{
  // open file and only allow reading
  std::ifstream in( filename.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to read data "
                     + filename );
    return;
  }

  // map for all voronoi cells in this time step
  std::map<int, SVoronoiCell*> cells;

  // current line of file
  std::string line;
  getline( in, line );

  while( in.good() )
  {
    // stop loop if the content is empty since we have
    // reached the end of the data
    if( line == "" )
      break;

    // for each new voronoi cell get the cell id and assign
    // the corresponding voronoi edges
    std::vector<std::string> tmpCellLine;
    boost::split( tmpCellLine, line, boost::is_any_of(" ") );
    int cellId = boost::lexical_cast<int>( tmpCellLine[2] );

    SVoronoiCell *vCell = new SVoronoiCell( cellId, timeStep );
    EdgeVector edges;

    // next line
    getline( in, line );

    while( line.at(0) != '/' )
    {
      // ignore sphere coordinates
      // just store the cylinder coordinates
      if( line.at(0) == 'c' )
      {
        // first split between the brackets
        std::vector<std::string> tmpParam;
        boost::split( tmpParam, line, boost::is_any_of("{}") );

        // then split between the commas
        std::vector<std::string> tmpCoord;
        boost::split( tmpCoord, tmpParam[1], boost::is_any_of(",") );

        osg::Vec3 pos1, pos2;

        // store the first position, erasing the first bracket
        pos1[0] = boost::lexical_cast<double>( tmpCoord[0].erase(0,1) );
        // the second entry just needs a conversion
        pos1[1] = boost::lexical_cast<double>( tmpCoord[1] );
        // for the last entry, erase the last bracket
        pos1[2] = boost::lexical_cast<double>( tmpCoord[2].erase(tmpCoord[2].size()-1,1) );

        // store the first position, erasing the first bracket
        pos2[0] = boost::lexical_cast<double>( tmpCoord[3].erase(0,1) );
        // the second entry just needs a conversion
        pos2[1] = boost::lexical_cast<double>( tmpCoord[4] );
        // for the last entry, erase the last bracket
        pos2[2] = boost::lexical_cast<double>( tmpCoord[5].erase(tmpCoord[5].size()-1,1) );

        edges.push_back( std::make_pair( pos1, pos2 ) );
      }

      getline( in, line );

      // stop loop if the content is empty since we have
      // reached the end of the data
      if( line == "" )
        break;
    }

    // set the voronoi edges
    vCell->setEdges( edges );

    // store this voronoi cell in the map
    cells[cellId] = vCell;
  }

  // store number of voronoi cells in the current time step
  _numVornoiCells.push_back( cells.size() );

  // at last store the complete map for this time step
  _voronoiCells->at(timeStep) = cells;

  in.close();
}

// ---------------------------------------------------------------------

void SVoronoiDiagram::readBBoxInformation( const std::string &filename,
                                           const int timeStep )
{
  // open bbox file and only allow reading
  std::ifstream in( filename.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to read data "
                     + filename );
    return;
  }

  std::string line;

  getline( in, line );

  std::vector<std::string> bb;
  std::vector<double> bbV;
  boost::split( bb, line, boost::is_any_of(" ") );
  for( int i = 0; i < 6; i++ )
    bbV.push_back( boost::lexical_cast<double>(bb[i]) );

  osg::BoundingBox bbox( bbV[0], bbV[2], bbV[4], bbV[1], bbV[3], bbV[5] );

  _bboxes.at(timeStep) = bbox;

  in.close();
}

// ---------------------------------------------------------------------

void SVoronoiDiagram::readNeighborAndVolumeInformation( const std::string &filename,
                                                        const int timeStep )
{
  // open bbox file and only allow reading
  std::ifstream in( filename.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to read data "
                     + filename );
    return;
  }

  std::string line;

  while( in.good() )
  {
    getline( in, line );

    if( line == "" )
      break;

    std::set<int> neighborIds;
    // first get the cell id and the volume of its
    // voronoi cell, then store all neighbor cell ids
    // the different categories are separated by the
    // '_' sign else a blank is used

    std::vector<std::string> cellInfos;
    boost::split( cellInfos, line, boost::is_any_of(" _") );

    // i == 0 -> cellid
    int cellId = boost::lexical_cast<int>( cellInfos[0] );

    // i == 1 -> volume of voronoi cell
    double volume = boost::lexical_cast<double>( cellInfos[1] );

    for( std::size_t i = 2; i < cellInfos.size(); i++ )
      neighborIds.insert( boost::lexical_cast<int>( cellInfos[i] ) );

    _voronoiCells->at(timeStep).at(cellId)->setVolume( volume );
    _voronoiCells->at(timeStep).at(cellId)->setNeighborIds( neighborIds );
  }

  in.close();
}

// ---------------------------------------------------------------------

void SVoronoiDiagram::storeVoronoiCells( const std::string &fileName )
{
  std::ofstream outfile( fileName.c_str(), std::ios::out | std::ios::binary );

  // save data to archive
  boost::archive::binary_oarchive oa( outfile );
  // write class instance to archive
  oa << *this;
}

// ---------------------------------------------------------------------

void SVoronoiDiagram::loadVoronoiCells( const std::string &fileName )
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
