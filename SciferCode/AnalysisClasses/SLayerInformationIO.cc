#include "SLayerInformationIO.hh"

/**
  @file   SLayerInformationIO.cc
  @brief  Contains class for I/O operations of layer information
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "common/VException.hh"

#include "kernel/SKernel.hh"

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/limits.hpp>
#include <boost/numeric/conversion/bounds.hpp>

namespace SLayerInformationIO
{
// ---------------------------------------------------------------------

void writeVolumes( const std::string &filename,
                   const boost::shared_ptr<cellLayerVector> &cellLayersPerTimeStep,
                   const std::vector< std::vector<double> > &layerVolumes,
                   const std::size_t numLayers )
{
  // open file at the end of the stream and allow writing
  std::ofstream out( filename.c_str(), std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to write into data "
                     + filename );
    return;
  }

  out << "Time ";

  for( std::size_t l = 0; l < numLayers; l++ )
    out << "Layer" << boost::lexical_cast<std::string>(l+1)
        << " Cells ";

  out << "MaxCells\n";

  // for each layer output the temporal development of
  // the volume sizes
  for( std::size_t t = 0; t < layerVolumes.size(); t++ )
  {
    std::size_t maxCells = 0;

    out << t+1 << " ";

    for( std::size_t l = 0; l < numLayers; l++ )
    {
      if( l < layerVolumes.at(t).size() )
      {
        out << layerVolumes.at(t).at(l) << " "
            << cellLayersPerTimeStep->at(t).at(l).size() << " ";

        maxCells += cellLayersPerTimeStep->at(t).at(l).size();
      }
      else
        out << "0 0 ";
    }
    out << maxCells << "\n";
  }

  out.close();

  std::cout << "Volume info is stored in " << filename << "." << std::endl;
}

// ---------------------------------------------------------------------

void writeVolume( const std::string &filename,
                  const std::vector< std::pair<double,std::size_t> > &volumes )
{
  // open file at the end of the stream and allow writing
  std::ofstream out( filename.c_str(), std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to write into data "
                     + filename );
    return;
  }

  // for each time step output the temporal development of
  // the volume size and the number of total cells at this time step
  for( std::size_t t = 0; t < volumes.size(); t++ )
    out << volumes.at(t).first << " " << volumes.at(t).second << "\n";

  out.close();

  std::cout << "Volume info is stored in " << filename << "." << std::endl;
}

// ---------------------------------------------------------------------

void writeCellPositions( const std::string &filename,
                         const boost::shared_ptr<SCellLayers> cellLayers )
{
  // this function is used for exporting the bbox and cell position
  // information such that an external program computes the voronoi diagram
  boost::shared_ptr<cellTimeVector> cellsPerTimeStep =
      cellLayers->getCellsPerTimeStep();

  for( std::size_t t = 0; t < cellsPerTimeStep->size(); t++ )
  {
    std::string bboxFN = filename + "/bbox/vT" + boost::lexical_cast<std::string>(t) + ".txt";
    std::string FN = filename + "/voronoiCells/vT" + boost::lexical_cast<std::string>(t) + ".txt";
    // open file at the end of the stream and allow writing
    std::ofstream out( FN.c_str(), std::ofstream::out );

    if( !out.is_open() )
    {
      THROW_EXCEPTION( VInvalidFilenameException, "Unable to write into data "
                       + FN );
      return;
    }

    std::vector<double> min;
    std::vector<double> max;
    min.resize( 3, boost::numeric::bounds<double>::highest() );
    max.resize( 3, boost::numeric::bounds<double>::lowest() );

    for( cellSet::const_iterator iter = cellsPerTimeStep->at(t).begin();
         iter != cellsPerTimeStep->at(t).end(); ++iter )
    {
      int cellId = (*iter)->cellId;
      osg::Vec3 pos( (*iter)->getX(), (*iter)->getY(), (*iter)->getZ() );
      out << cellId << " " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";

      for( int i = 0; i < 3; i++ )
      {
        if( pos[i] < min[i] )
          min[i] = pos[i];

        if( pos[i] >= max[i] )
          max[i] = pos[i];
      }
    }

    out.close();

    // open bbox file
    std::ofstream outbb( bboxFN.c_str(), std::ofstream::out );

    if( !outbb.is_open() )
    {
      THROW_EXCEPTION( VInvalidFilenameException, "Unable to write into data "
                       + bboxFN );
      return;
    }

    // before the min/max values are written, the min values are slightly
    // decreased while the max values are slightly increased such that
    // the external computation of the voronoi cells also includes the
    // cells located at the boundary
    for( int i = 0; i < 3; i++ )
    {
      min[i] -= 5.;
      max[i] += 5.;
    }

    // output the bounding box info
    outbb << min[0] << " " << max[0] << " "
          << min[1] << " " << max[1] << " "
          << min[2] << " " << max[2] << "\n";

    outbb.close();
  }
}

// ---------------------------------------------------------------------

// WARNING: This method should only be called if no cell layer
// information was already appended in the data as a column else it will
// erase the division type column or anything else at the back
void writeLayerInformation( const std::string &filename,
                            const boost::shared_ptr<cellLayerVector> cellLayers )
{
  // open file
  std::ifstream in( filename.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data "
                     + filename );
    return;
  }

  // current line of file
  std::string line;

  // store complete file contents in vector of strings
  std::vector<std::string> oData;

  while( in.good() )
  {
    getline( in, line );

    // stop loop if the content is empty since we have
    // reached the end of the data
    if( line == "" )
      break;

    oData.push_back( line );
  }

  in.close();

  // after reading create new data file with appending layer
  // assignment in each column
  std::fstream out( filename.c_str(), std::ofstream::in | std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data "
                     + filename );
    return;
  }

  // skip the first three header lines
  // and get the forth line with the descriptions
  std::size_t begin = 3;
  for( std::size_t i = 0; i < begin; i++ )
    getline( out, line );

  // for each cell there is always a pair of entries
  // describing the cell movement, thus the cell id
  // is the same for two lines which is checked here
  // because those two lines for one cell is assigned
  // the same layer value
  std::size_t lastId = 0;
  std::size_t currentId = 0;
  std::size_t lastLayer = 0;

  std::size_t currentLine = begin;
  while( currentLine < oData.size() )
  {
    // get only cell id and time information
    std::stringstream oLineStream( oData.at(currentLine),
                                   std::ios_base::app | std::ios_base::out );
    std::stringstream readLineStream( oData.at(currentLine) );

    // the first line of data is the descripton information
    // which is appended by the layer column
    if( currentLine == begin )
    {
      oLineStream << " " << "Layer\n";
      out << oLineStream.str();
      currentLine++;
      continue;
    }

    int id, time;
    double dummy;
    int layer = -1;

    readLineStream >> id >> dummy >> dummy >> dummy >> time;

    currentId = id;
    bool sameCell = false;
    if( currentId == lastId )
    {
      sameCell = true;
    }
    else
    {
      lastId = currentId;
      sameCell = false;
    }

    // if we are situated at the destiny position of a cell
    // then we already know the layer value
    if( !sameCell )
    {
      // get layer value of current cell id and time step
      // loop over all layers in the current time step
      for( std::size_t l = 0; l < cellLayers->at(time-1).size(); l++ )
      {
        // loop over all cells in the layer and time step
        for( cellSet::iterator iter = cellLayers->at(time-1).at(l).begin();
             iter != cellLayers->at(time-1).at(l).end(); ++iter )
        {
          if( (*iter)->cellId == id )
          {
            layer = l;
            break;
          }
        }

        if( layer != -1 )
          break;
      }

      // update last layer
      lastLayer = layer;
    }
    // else set current layer to the last layer value
    else
      layer = lastLayer;

    // append layer information to original stream
    oLineStream << " " << layer << "\n";
    out << oLineStream.str();

    currentLine++;
  }

  out.close();

  std::cout << "Data is extented by layer information as column." << std::endl;
}

// ---------------------------------------------------------------------

// WARNING: This method should only be called if some other column
// information is appended after the layer information (which is
// for example the division type column)
void overwriteLayerInformation( const std::string &filename,
                                const boost::shared_ptr<cellLayerVector> cellLayers )
{
  // open file
  std::ifstream in( filename.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data "
                     + filename );
    return;
  }

  // current line of file
  std::string line;

  // store complete file contents in vector of strings
  std::vector<std::string> oData;
  // since the layer information is not the last
  // column store everything in the line after the
  // layer information
  std::vector<std::string> backData;

  // get the initial four lines with header info only
  for( std::size_t i = 0; i < 4; i++ )
  {
    getline( in, line );
    oData.push_back( line );

    // for the header info just push an empty string
    // since oData already includes everything
    backData.push_back( "" );
  }

  while( in.good() )
  {
    getline( in, line );

    // stop loop if the content is empty since we have
    // reached the end of the data
    if( line == "" )
      break;

    // get last entry of line since and store it;
    // erase the next entry
    int lastEntryLength = 0;
    std::string backString;
    std::string::reverse_iterator sIter = line.rbegin();
    int positionIndex = 0;
    for( ;sIter != line.rend(); sIter++ )
    {
      if( *sIter == ' ' )
      {
        positionIndex++;

        if( positionIndex == 1 )
          break;
      }

      backString += *sIter;
      lastEntryLength++;
    }

    // swap content
    backString = std::string( backString.rbegin(), backString.rend() );

    line.erase( line.end() - lastEntryLength - 3, line.end() );

    oData.push_back( line );
    backData.push_back( backString );
  }

  in.close();

  // after reading create new data file with appending layer
  // assignment in each column
  std::fstream out( filename.c_str(), std::ofstream::in | std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data "
                     + filename );
    return;
  }

  // skip the first three header lines
  // and get the forth line with the descriptions
  std::size_t begin = 4;
  for( std::size_t i = 0; i < begin; i++ )
    getline( out, line );

  // for each cell there is always a pair of entries
  // describing the cell movement, thus the cell id
  // is the same for two lines which is checked here
  // because those two lines for one cell is assigned
  // the same layer value
  std::size_t lastId = 0;
  std::size_t currentId = 0;
  std::size_t lastLayer = 0;

  std::size_t currentLine = begin;
  while( currentLine < oData.size() )
  {
    // get only cell id and time information
    std::stringstream oLineStream( oData.at(currentLine),
                                   std::ios_base::app | std::ios_base::out );
    std::stringstream readLineStream( oData.at(currentLine) );

    int id;
    std::size_t time;
    double dummy;
    int layer = -1;

    readLineStream >> id >> dummy >> dummy >> dummy >> time;

    currentId = id;
    bool sameCell = false;
    if( currentId == lastId )
    {
      sameCell = true;
    }
    else
    {
      lastId = currentId;
      sameCell = false;
    }

    // if we are situated at the destiny position of a cell
    // then we already know the layer value
    if( !sameCell )
    {
      // get layer value of current cell id and time step
      // loop over all layers in the current time step
      for( std::size_t l = 0; l < cellLayers->at(time-1).size(); l++ )
      {
        // loop over all cells in the layer and time step
        for( cellSet::iterator iter = cellLayers->at(time-1).at(l).begin();
             iter != cellLayers->at(time-1).at(l).end(); ++iter )
        {
          if( (*iter)->cellId == id )
          {
            layer = l;
            break;
          }
        }

        if( layer != -1 )
          break;
      }

      // update last layer
      lastLayer = layer;
    }
    // else set current layer to the last layer value
    else
      layer = lastLayer;

    // append layer information to original stream
    oLineStream << " " << layer << " " << backData.at(currentLine) << "\n";
    out << oLineStream.str();

    currentLine++;
  }

  out.close();

  std::cout << "Layer information saved." << std::endl;
}

// ---------------------------------------------------------------------

// This method is only valid to be called if the raw data includes
// the information of each individual time step. If only start and
// end positions are given then the method 'readLayerInformation'
// have to be called
void readLayerInformationPerTimeStep( const std::string &filename,
                                      boost::shared_ptr<SCellLayers> cellLayers,
                                      const int layerColumnIndex )
{
  // open file and only allow reading
  std::ifstream in( filename.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to read data "
                     + filename );
    return;
  }

  std::string line;

  // skip the first forth header lines
  for( std::size_t i = 0; i < 4; i++ )
    getline( in, line );

  boost::shared_ptr<cellTimeVector> ctv = cellLayers->getCellsPerTimeStep();
  boost::shared_ptr<cellLayerVector> clv = cellLayers->getCellLayersPerTimeStep();

  // resize vector depending on number of time steps
  clv->resize( ctv->size() );

  // get max number of layers
  std::size_t maxNumberOfLayers = 1;

  int step = 0;

  // Now read the file line per line
  while( in.good() )
  {
    SKernel::get()->setProgress( (double)step/((double)ctv->size() - 1), "Reading layer information" );

    getline( in, line );

    if( line == "" )
      break;

    std::stringstream lineStream( line );

    int id, time;
    double dummy;

    // get cell id and time info
    lineStream >> id >> dummy >> dummy >> dummy >> time;

    // do nothing for the first time step since
    // this is already done in the initialization
    // of the cell layers
    if( time != 1 )
    {
      std::size_t layer = assignLayerValue( id, time, line, ctv, clv, layerColumnIndex );

      if( maxNumberOfLayers < layer + 1 )
        maxNumberOfLayers = layer + 1;
    }

    step++;
  }

  // at last update maximum number of layers
  cellLayers->setMaxNumLayers( maxNumberOfLayers );

  in.close();

  std::cout << "Layer information loaded." << std::endl;
}

// ---------------------------------------------------------------------

void readLayerInformation( const std::string &filename,
                           boost::shared_ptr<SCellLayers> cellLayers,
                           const int layerColumnIndex )
{
  // open file and only allow reading
  std::ifstream in( filename.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to read data "
                     + filename );
    return;
  }

  std::string line, line1, line2;

  // skip the first forth header lines
  for( std::size_t i = 0; i < 4; i++ )
    getline( in, line );

  boost::shared_ptr<cellTimeVector> ctv = cellLayers->getCellsPerTimeStep();
  boost::shared_ptr<cellLayerVector> clv = cellLayers->getCellLayersPerTimeStep();

  // resize vector depending on number of time steps
  clv->resize( ctv->size() );

  // get max number of layers
  std::size_t maxNumberOfLayers = 1;

  int step = 0;

  // Now read the file line per line
  while( in.good() )
  {
    SKernel::get()->setProgress( (double)step/((double)ctv->size() - 1), "Reading layer information" );

    // always read a pair of lines
    getline( in, line1 );
    getline( in, line2 );

    if( line1 == "" )
      break;

    std::stringstream lineStream1( line1 );
    std::stringstream lineStream2( line2 );

    // id and layer are identical for the pair of lines
    // if it is a cell movement; if the cell divides immediately
    // then there is only one entry for that cell
    // but the time steps differ
    int id1, id2, time1, time2;
    double dummy;

    // get cell id and time info
    lineStream1 >> id1 >> dummy >> dummy >> dummy >> time1;
    lineStream2 >> id2 >> dummy >> dummy >> dummy >> time2;

    std::size_t layer;

    // check if the current cell is a cell that only
    // exists for one time step.
    // if so, then assign it its layer value (if it in time step != 1)
    // and read the next line
    // since this division node was already assigned
    // to layer 0 in time step 0
    if( id1 != id2 )
    {
      while( id1 != id2 )
      {
        // process the first line
        // get layer value and assign it to the cell that exists
        // only in one time step
        layer = assignLayerValue( time1, time1, id1, line1, ctv, clv, layerColumnIndex );

        // read next line and update id and time step values
        time1 = time2;
        id1 = id2;
        line1 = line2;
        getline( in, line2 );
        std::stringstream lineStream( line2 );

        lineStream >> id2 >> dummy >> dummy >> dummy >> time2;
      }
    }

    layer = assignLayerValue( time1, time2, id1, line1, ctv, clv, layerColumnIndex );

    if( maxNumberOfLayers < layer + 1 )
      maxNumberOfLayers = layer + 1;

    step++;
  }

  // at last update maximum number of layers
  cellLayers->setMaxNumLayers( maxNumberOfLayers );

  in.close();

  std::cout << "Layer information loaded." << std::endl;
}

// ---------------------------------------------------------------------

int assignLayerValue( const std::size_t tStart, const std::size_t tEnd,
                      const int id, const std::string &line,
                      const boost::shared_ptr<cellTimeVector> ctv,
                      boost::shared_ptr<cellLayerVector> clv,
                      const int layerColumnIndex )
{
  // get layer info only for one line since they are identical
  // The value is the pre last entry
  std::string layerString;
  std::string::const_reverse_iterator sIter = line.rbegin();
  // the position index denotes the word (layer) number counting
  // from the back of each line
  int positionIndex = 0;
  for( ;sIter != line.rend(); sIter++ )
  {
    if( *sIter == ' ' )
    {
      positionIndex++;
      sIter++;

      if( positionIndex == layerColumnIndex+1 )
        break;
    }

    // store the entry after the first word
    if( positionIndex == layerColumnIndex )
      layerString += *sIter;
  }

  // swap content
  layerString = std::string( layerString.rbegin(), layerString.rend() );

  std::size_t layer = boost::lexical_cast<std::size_t>( layerString );

  // loop over all time steps in which the cells move
  for( std::size_t i = tStart; i <= tEnd; i++ )
  {
    // do nothing for the first time step since
    // this is already done in the initialization
    // of the cell layers
    if( i == 1 )
      continue;

    // before setting the new layer value we check if we have to
    // resize the cell layer vector
    if( clv->at( i-1 ).size() <= layer )
      clv->at( i-1 ).resize( layer+1 );

    for( cellSet::const_iterator treeIter = ctv->at( i-1 ).begin();
         treeIter != ctv->at( i-1 ).end(); treeIter++ )
    {
      if( (*treeIter)->cellId == id && (*treeIter)->timeStep == i )
      {
        clv->at( i-1 ).at( layer ).insert( *treeIter );
        break;
      }
    }
  }

  return layer;
}

// ---------------------------------------------------------------------

int assignLayerValue( const int id,
                      const std::size_t timeStep,
                      const std::string &line,
                      const boost::shared_ptr<cellTimeVector> ctv,
                      boost::shared_ptr<cellLayerVector> clv,
                      const int layerColumnIndex )
{
  // get layer info only for one line since they are identical
  // The value is the pre last entry
  std::string layerString;
  std::string::const_reverse_iterator sIter = line.rbegin();
  // the position index denotes the word (layer) number counting
  // from the back of each line
  int positionIndex = 0;
  for( ;sIter != line.rend(); sIter++ )
  {
    if( *sIter == ' ' )
    {
      positionIndex++;
      sIter++;

      if( positionIndex == layerColumnIndex+1 )
        break;
    }

    // store the entry after the first word
    if( positionIndex == layerColumnIndex )
      layerString += *sIter;
  }

  // swap content
  layerString = std::string( layerString.rbegin(), layerString.rend() );

  std::size_t layer = boost::lexical_cast<std::size_t>( layerString );

  // before setting the new layer value we check if we have to
  // resize the cell layer vector
  if( clv->at( timeStep-1 ).size() <= layer )
    clv->at( timeStep-1 ).resize( layer+1 );

  for( cellSet::const_iterator treeIter = ctv->at( timeStep-1 ).begin();
       treeIter != ctv->at( timeStep-1 ).end(); treeIter++ )
  {
    if( (*treeIter)->cellId == id && (*treeIter)->timeStep == timeStep )
    {
      clv->at( timeStep-1 ).at( layer ).insert( *treeIter );
      break;
    }
  }

  return layer;
}

// ---------------------------------------------------------------------

std::string readDivisionType( const std::string &line )
{
  // get division type info which is the last column entry
  std::string divTypeString;
  std::string::const_reverse_iterator sIter = line.rbegin();
  for( ;sIter != line.rend(); sIter++ )
  {
    if( *sIter == ' ' )
      break;

    divTypeString += *sIter;
  }

  // swap content
  divTypeString = std::string( divTypeString.rbegin(), divTypeString.rend() );

  return divTypeString;
}

// ---------------------------------------------------------------------

void writeDivisionScheme( const std::string &filename,
                          const boost::shared_ptr<cellTimeVector> cTV,
                          const NodeFeatureInfo &divisionScheme )
{
  // open file
  std::ifstream in( filename.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data "
                     + filename );
    return;
  }

  // current line of file
  std::string line;

  // store complete file contents in vector of strings
  std::vector<std::string> oData;

  while( in.good() )
  {
    getline( in, line );

    // stop loop if the content is empty since we have
    // reached the end of the data
    if( line == "" )
      break;

    oData.push_back( line );
  }

  in.close();

  // after reading create new data file with appending layer
  // assignment in each column
  std::fstream out( filename.c_str(), std::ofstream::in | std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data "
                     + filename );
    return;
  }

  // skip the first three header lines
  // and get the forth line with the descriptions
  std::size_t begin = 3;
  for( std::size_t i = 0; i < begin; i++ )
    getline( out, line );

  std::size_t currentLine = begin;
  while( currentLine < oData.size() )
  {
    // get only cell id and time information
    std::stringstream oLineStream( oData.at(currentLine),
                                   std::ios_base::app | std::ios_base::out );
    std::stringstream readLineStream( oData.at(currentLine) );

    // the first line of data is the descripton information
    // which is appended by the division type column
    if( currentLine == begin )
    {
      oLineStream << " " << "DivisionType\n";
      out << oLineStream.str();
      currentLine++;
      continue;
    }

    int id;
    std::size_t time;
    double dummy;
    bool divFound = false;

    readLineStream >> id >> dummy >> dummy >> dummy >> time;

    // check all nodes in the current time step
    for( cellSet::const_iterator iter = cTV->at(time-1).begin();
         iter != cTV->at(time-1).end(); ++iter )
    {
      // check if the current node is a division node
      if( (*iter)->cellId == id && (*iter)->children.size() == 2 )
      {
        // append division type information to original stream
        // 0 -> anticlinal
        // 1 -> periclinal
        // 2 -> radial
        int divType = divisionScheme.at(time-1).at(id);
        oLineStream << " " << divType << "\n";

        divFound = true;
        break;
      }
    }

    if( !divFound )
       oLineStream << " -\n";

    out << oLineStream.str();

    currentLine++;
  }

  out.close();

  std::cout << "Data is extented by division type information as column." << std::endl;
}

// ---------------------------------------------------------------------

void overwriteDivisionScheme( const std::string &filename,
                              const NodeFeatureInfo &divisionScheme )
{
  // open file
  std::ifstream in( filename.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data "
                     + filename );
    return;
  }

  // current line of file
  std::string line;

  // store complete file contents in vector of strings
  std::vector<std::string> oData;

  // get the initial four lines with header info only
  for( int i = 0; i < 4; i++ )
  {
    getline( in, line );
    oData.push_back( line );
  }

  while( in.good() )
  {
    getline( in, line );

    // stop loop if the content is empty since we have
    // reached the end of the data
    if( line == "" )
      break;

    // erase last entry of line since we write
    // the new division type information into it
    int lastEntryLength = 0;
    std::string::reverse_iterator sIter = line.rbegin();
    for( ;sIter != line.rend(); sIter++ )
    {
      if( *sIter == ' ' )
        break;

      lastEntryLength++;
    }

    line.erase( line.end() - lastEntryLength - 1, line.end() );

    oData.push_back( line );
  }

  in.close();

  // after reading create new data file with appending
  // the division type in each row
  std::fstream out( filename.c_str(), std::ofstream::in | std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data "
                     + filename );
    return;
  }

  // skip the first three header lines
  // and get the forth line with the descriptions
  std::size_t begin = 4;
  for( std::size_t i = 0; i < begin; i++ )
    getline( out, line );

  std::size_t currentLine = begin;
  while( currentLine < oData.size() )
  {
    // get only cell id and time information
    std::stringstream oLineStream( oData.at(currentLine),
                                   std::ios_base::app | std::ios_base::out );
    std::stringstream readLineStream( oData.at(currentLine) );

    std::size_t id, time;
    double dummy;

    readLineStream >> id >> dummy >> dummy >> dummy >> time;

    std::map<std::size_t, int>::const_iterator mapIter =
        divisionScheme.at(time-1).find(id);

    if( mapIter != divisionScheme.at(time-1).end() )
    {
      int divType = mapIter->second;
      // append division type information to original stream
      oLineStream << " " << divType << "\n";
    }
    // else include a '-' sign
    else
      oLineStream << " -\n";

    out << oLineStream.str();

    currentLine++;
  }

  out.close();

  std::cout << "Division type information saved." << std::endl;
}

// ---------------------------------------------------------------------

void readDivisionScheme( const std::string &filename,
                         const std::size_t numTimeSteps,
                         NodeFeatureInfo &divisionScheme )
{
  // open file and only allow reading
  std::ifstream in( filename.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to read data "
                     + filename );
    return;
  }

  divisionScheme.resize( numTimeSteps );

  std::string line;

  // skip the first forth header lines
  for( int i = 0; i < 4; i++ )
    getline( in, line );

  // Now read the file line per line
  while( in.good() )
  {
    // always read a pair of lines
    getline( in, line );

    if( line == "" )
      break;

    std::stringstream lineStream( line );

    // id and layer are identical for the pair of lines
    // if it is a cell movement; if the cell divides immediately
    // then there is only one entry for that cell
    // but the time steps differ
    int id;
    std::size_t time;
    double dummy;

    // get cell id and time info
    lineStream >> id >> dummy >> dummy >> dummy >> time;

    std::string divTypeString = readDivisionType( line );

    if( divTypeString == "-" )
      continue;

    int divType = boost::lexical_cast<int>( divTypeString );
    divisionScheme.at(time-1).insert( std::make_pair( id, divType ) );
  }

  in.close();
}

// ---------------------------------------------------------------------

void generateDivisionScheme( const boost::shared_ptr<SAdvancedLineageTreeCollection> lineageTrees,
                             const NodeFeatureInfo &cellLayers,
                             NodeFeatureInfo &divisionScheme,
                             const osg::Vec3 &sideViewNormal,
                             const double radialAngleThreshold,
                             const osg::Matrix &rotMat )
{
  divisionScheme.resize( lineageTrees->getOrComputeNumberOfTimesteps() );

  for( SAdvancedLineageTreeCollection::const_iterator l = lineageTrees->begin();
       l != lineageTrees->end(); ++l )
  {
    SLineageTree *tree = l->second;

    // traverse the whole tree
    for( SLineageTree::iterator iter = tree->begin();
         iter != tree->end(); ++iter )
    {
      if( iter->children.size() == 2 )
      {
        // divPattern defines the type of division
        // 0 -> anticlinal
        // 1 -> periclinal
        // 2 -> radial
        int divPattern = 0;

        // get layer values of parent and its children nodes
        int pLayer = cellLayers.at(iter->timeStep-1).find(iter->cellId)->second;
        int cLayer0 = cellLayers.at(iter->children[0]->timeStep-1).find(iter->children[0]->cellId)->second;
        int cLayer1 = cellLayers.at(iter->children[1]->timeStep-1).find(iter->children[1]->cellId)->second;

        // periclinal division if the layer of parent is different to layer of child
        if( pLayer != cLayer0 || pLayer != cLayer1 )
          divPattern = 1;
        // else check if it is anticlinal or radial division by comparing the division direction
        // with the normal of the side view. Note that the side view is set manually by the rotation
        // parameters. TODO: Do this in future by computing the principal component of the master
        // cell file -> perhaps better approximation
        else
        {
          osg::Vec3 cPos1( iter->children[0]->getX(), iter->children[0]->getY(), iter->children[0]->getZ() );
          osg::Vec3 cPos2( iter->children[1]->getX(), iter->children[1]->getY(), iter->children[1]->getZ() );

          // check which direction should be chosen: cPos1 -> cPos2 or cPos2 -> cPos1?
          osg::Vec3 c1 = cPos1;
          osg::Vec3 c2 = cPos2;
          c1 = c1 * rotMat;
          c2 = c2 * rotMat;

          // get direction vector between both cells
          // set the smaller y value to be the cell
          osg::Vec3 childrenDir;
          if( c1[1] < c2[1] )
            childrenDir = cPos1 - cPos2;
          else
            childrenDir = cPos2 - cPos1;

          childrenDir.normalize();

          // compute angle
          double angle = acos( sideViewNormal * childrenDir );
          angle = angle * 180./M_PI;

          // if smaller than threshold then we set a radial division
          if( angle < radialAngleThreshold )
            divPattern = 2;
        }

        // store division pattern
        divisionScheme.at( iter->timeStep-1 ).insert( std::make_pair( iter->cellId, divPattern ) );
      }
    }
  }
}

// ---------------------------------------------------------------------

void generateDivisionScheme( const boost::shared_ptr<SAdvancedLineageTreeCollection> lineageTrees,
                             const boost::shared_ptr<SCellLayers> &cellLayers,
                             NodeFeatureInfo &divisionScheme,
                             const osg::Vec3 &sideViewNormal,
                             const double radialAngleThreshold,
                             const osg::Matrix &rotMat )
{
  divisionScheme.resize( lineageTrees->getOrComputeNumberOfTimesteps() );

  boost::shared_ptr< std::vector<cellHistoryMap> > cH = cellLayers->getCellHistories();

  for( SAdvancedLineageTreeCollection::const_iterator l = lineageTrees->begin();
       l != lineageTrees->end(); ++l )
  {
    SLineageTree *tree = l->second;

    // traverse the whole tree
    for( SLineageTree::iterator iter = tree->begin();
         iter != tree->end(); ++iter )
    {
      if( iter->children.size() == 2 )
      {
        // divPattern defines the type of division
        // 0 -> anticlinal
        // 1 -> periclinal
        // 2 -> radial
        int divPattern = 0;

        // get layer values of parent and its children nodes
        std::size_t pLayer = cH->at(iter->timeStep-1).find(*iter)->second.size();
        std::size_t cLayer0 = cH->at(iter->children[0]->timeStep-1).find(iter->children[0])->second.size();
        std::size_t cLayer1 = cH->at(iter->children[1]->timeStep-1).find(iter->children[1])->second.size();

        // periclinal division if the layer of parent is different to layer of child
        if( pLayer != cLayer0 || pLayer != cLayer1 )
          divPattern = 1;
        // else check if it is anticlinal or radial division by comparing the division direction
        // with the normal of the side view. Note that the side view is set manually by the rotation
        // parameters. TODO: Do this in future by computing the principal component of the master
        // cell file -> perhaps better approximation
        else
        {
          osg::Vec3 cPos1( iter->children[0]->getX(), iter->children[0]->getY(), iter->children[0]->getZ() );
          osg::Vec3 cPos2( iter->children[1]->getX(), iter->children[1]->getY(), iter->children[1]->getZ() );

          // check which direction should be chosen: cPos1 -> cPos2 or cPos2 -> cPos1?
          osg::Vec3 c1 = cPos1;
          osg::Vec3 c2 = cPos2;
          c1 = c1 * rotMat;
          c2 = c2 * rotMat;

          // get direction vector between both cells
          // set the smaller y value to be the cell
          osg::Vec3 childrenDir;
          if( c1[1] < c2[1] )
            childrenDir = cPos1 - cPos2;
          else
            childrenDir = cPos2 - cPos1;

          childrenDir.normalize();

          // compute angle
          double angle = acos( sideViewNormal * childrenDir );
          angle = angle * 180./M_PI;

          // if smaller than threshold then we set a radial division
          if( angle < radialAngleThreshold )
            divPattern = 2;
        }

        // store division pattern
        divisionScheme.at( iter->timeStep-1 ).insert( std::make_pair( iter->cellId, divPattern ) );
      }
    }
  }
}

// ---------------------------------------------------------------------

void generatePartialDivisionScheme( const boost::shared_ptr<SAdvancedLineageTreeCollection> lineageTrees,
                                    const NodeFeatureInfo &cellLayers,
                                    NodeFeatureInfo &divisionScheme )
{
  divisionScheme.resize( lineageTrees->getOrComputeNumberOfTimesteps() );

  for( SAdvancedLineageTreeCollection::const_iterator l = lineageTrees->begin();
       l != lineageTrees->end(); ++l )
  {
    SLineageTree *tree = l->second;

    // traverse the whole tree
    for( SLineageTree::iterator iter = tree->begin();
         iter != tree->end(); ++iter )
    {
      if( iter->children.size() == 2 )
      {
        // divPattern defines the type of division
        // 0 -> anticlinal
        // 1 -> periclinal
        int divPattern = 0;

        // get layer values of parent and its children nodes
        int pLayer = cellLayers.at(iter->timeStep-1).find(iter->cellId)->second;
        int cLayer0 = cellLayers.at(iter->children[0]->timeStep-1).find(iter->children[0]->cellId)->second;
        int cLayer1 = cellLayers.at(iter->children[1]->timeStep-1).find(iter->children[1]->cellId)->second;

        // periclinal division if the layer of parent is different to layer of child
        if( pLayer != cLayer0 || pLayer != cLayer1 )
          divPattern = 1;

        // store division pattern
        divisionScheme.at( iter->timeStep-1 ).insert( std::make_pair( iter->cellId, divPattern ) );
      }
    }
  }
}

// ---------------------------------------------------------------------

void readCellWallInformation( const std::string &filename,
                              CellWalls &cellWalls )
{
  // open file and only allow reading
  std::ifstream in( filename.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to read data "
                     + filename );
    return;
  }

  std::string line;

  // skip header information
  getline( in, line );

  std::size_t oldTime = 1;
  cellWalls.resize(1);

  // Now read the file line per line
  while( in.good() )
  {
    // always read a pair of lines
    getline( in, line );

    if( line == "" )
      break;

    std::stringstream lineStream( line );
    std::size_t id;
    std::size_t time;
    std::string neighbors;

    lineStream >> id >> time >> neighbors;

    std::vector<std::string> coords;
    boost::split( coords, neighbors, boost::is_any_of( "," ) );
    std::vector<osg::Vec3> array;
    for( std::size_t c=0; c<coords.size(); c+=3 )
    {
      osg::Vec3 temp( boost::lexical_cast<double>( coords.at(c) ),
                      boost::lexical_cast<double>( coords.at(c+1) ),
                      boost::lexical_cast<double>( coords.at(c+2) ) );

      array.push_back( temp );
    }

    if( time == oldTime )
      cellWalls.back().insert( std::pair<std::size_t,std::vector<osg::Vec3> >( id, array ) );
    else
    {
      std::map<std::size_t, std::vector<osg::Vec3> > map;
      map.insert( std::pair<std::size_t,std::vector<osg::Vec3> >( id, array ) );
      cellWalls.push_back( map );
      oldTime = time;
    }
  }

  in.close();
}

// ---------------------------------------------------------------------

}
