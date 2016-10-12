#include "SLoadArabidopsisLineages.hh"

#include "common/VException.hh"

#include "kernel/SKernel.hh"
#include "kernel/SAlgorithmProfile.hh"
#include "dataSet/SDataManager.hh"

#include <string>
#include <fstream>

#include <osg/MatrixTransform>

#include <boost/numeric/conversion/bounds.hpp>

using namespace std;

//----------------------------------------------------------------------------//

typedef SAlgorithm *maker_t();
extern std::map<std::string, maker_t*> algorithmFactory;

//----------------------------------------------------------------------------//

extern "C" 
{
  SAlgorithm *maker()
  {
    return new SLoadArabidopsisLineages();
  }

  class proxy 
  { 
  public:
    proxy()
    {
      algorithmFactory[SLoadArabidopsisLineages::algoName()] = maker;
    }
  };

  proxy p;
}

//----------------------------------------------------------------------------//

SLoadArabidopsisLineages::SLoadArabidopsisLineages() :
  _want2Dcopy( true ),
  _applyRotations( true ),
  _interpolateMovements( true )
{
  _data = boost::shared_ptr<SAdvancedLineageTreeCollection>(
        new SAdvancedLineageTreeCollection );
}

//----------------------------------------------------------------------------//

SLoadArabidopsisLineages::~SLoadArabidopsisLineages()
{

}

//----------------------------------------------------------------------------//

void SLoadArabidopsisLineages::initProfile()
{
  _profile->addFilename( _filename, "Input file" );
  _profile->addCheckBox( _want2Dcopy, "Create an abstract 2D copy" );
  _profile->addCheckBox( _applyRotations, "Apply rotation information" );
  _profile->addCheckBox( _interpolateMovements, "Apply movement interpolation" );
}

//----------------------------------------------------------------------------//

bool SLoadArabidopsisLineages::canUndo()
{
  return false;
}

//----------------------------------------------------------------------------//

void SLoadArabidopsisLineages::undo()
{

}

//----------------------------------------------------------------------------//

std::string SLoadArabidopsisLineages::getMenuEntry() const
{
  return std::string( "Algorithms/Load/" + getName() );
}

//----------------------------------------------------------------------------//

std::string SLoadArabidopsisLineages::algoName()
{
  return "Load Arabidopsis Lineages";
}

//----------------------------------------------------------------------------//

std::string SLoadArabidopsisLineages::getName() const
{
  return algoName();
}

//----------------------------------------------------------------------------//

void SLoadArabidopsisLineages::execute( void* )
{
  try {
    print();

//    filePreprocessing();
    loadLineages();

    // Get some statistics
//    int mints = _data->getMinimumTimestep();
//    int maxts = mints + _data->getOrComputeNumberOfTimesteps() - 1;
//    std::cout << "min time step: " << mints << " | max time step: " << maxts << std::endl;
//    std::cout << "num trees: " << _data->size() << std::endl;

    /* Convert the _filename string back to a boost::filesystem::path to be able
     * to use the stem() method later. */
    boost::filesystem::path p( _filename );

    // Use filename without path and extension as data name
    _data->setName( p.stem().string() );

    std::cout << "name of data set: " << _data->getName() << std::endl;

    // and store full path name for later writing and reading from the data
    _data->setFullPath( _filename );

     SDataManager::get()->add( _data );

     /* If desired make a 2D copy of the lineage trees e.g. for use in tracking
      * editor */
     if( _want2Dcopy )
     {
       boost::shared_ptr<SLineageTreeCollection> copy2D( new
         SLineageTreeCollection( *_data ) );
       copy2D->setName( p.stem().string() + "-abstr" );
       copy2D->adaptTreeLayouts( false );
       SDataManager::get()->add( copy2D );

       //_data->setName( _data->getName() + "-geom" );
     }
  }
  CATCH_N_RETHROW( VException );
}

//----------------------------------------------------------------------------//

void SLoadArabidopsisLineages::loadLineages()
{
  ifstream in( _filename.c_str() );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to read input file "
                     + _filename );
    return;
  }

  std::string line; // This will always contain one line of data from the file

  // data types for storing the cell layers and cell files
  NodeFeatureInfo cellLayers, cellDivisionTypes;
  NodeDoubleFeatureInfo cellAreas;
  std::map<int,int> cellFiles;
  std::pair<double, double> minMaxArea;
  minMaxArea.first = boost::numeric::bounds<double>::highest();
  minMaxArea.second = boost::numeric::bounds<double>::lowest();

  // store information if data is interpolated or not
  _data->interpolatedData( _interpolateMovements );

  // first line is the rotation
  getline( in, line );
  stringstream ssRotation( line );
  SArray rotation(3);
  ssRotation >> rotation[0] >> rotation[1] >> rotation[2];
  _data->setRotation( rotation );
  _data->applyRotation( _applyRotations );

  // second line is the center of the root
  getline( in, line );
  stringstream ssCenter( line );
  SArray center(3);
  ssCenter >> center[0] >> center[1] >> center[2];
  _data->setCenter( center );

  // third line is the dimension information of the raw data set
  getline( in, line );
  stringstream ssDim( line );
  SArray dimension(3);
  ssDim >> dimension[0] >> dimension[1] >> dimension[2];
  _data->setDimension( dimension );

  // Read forth line -> column description
  getline( in, line );

  // check if trackGroup, layer or division type column are stored in the data
  this->checkDataContent( line );

  // Remember previous node
  int prevId = 0; // Data in text file begins with 1
  int prevTime = 0;
  double prevX, prevY, prevZ;
  SLineageTree* prevNode = NULL;

  unsigned int lineNumber = 1; // To hint towards faulty data
  unsigned int tooManyChildrenCases = 0; // Usefull during development
  while( in.good() ) // Now read the file line per line
  {
    ++lineNumber;

    // Temporary variables for reading
    string dummyString;
    std::size_t time;
    int id, precursor;
    double x, y, z, rad;

    SLineageTree* node = new SLineageTree;
    std::list<int> precursors;
    int treeId = -1;


    /* Read the node information off the line */

    getline( in, line );

    if( line.size() == 0 )
      continue;

    stringstream lineStream( line );

    // ObjectID, X, Y, Z, Timepoint and Radius are easy
    lineStream >> id >> x >> y >> z >> time >> rad;

    // always z flip the data set
    z = -1.*z;

    // If there is more than one precursor, the format changes from {} or {p} to
    // "{p1, p2, ..., pn}" (note the quotation marks and the space).
    // So because the >> operator reads till it meets a space, we need to read
    // once more for each precursor.
    lineStream >> dummyString; // The critical piece of string
    // If dummyString == "{}" precursors list remains empty
    if( dummyString.at(0) == '{' && dummyString.at(1) != '}' ) // One precursor
    {
      stringstream ss( dummyString.substr( 1, dummyString.length()-2 ) );
      ss >> precursor;
      precursors.push_back( precursor );
    }

    if( dummyString.at(0) == '\"' ) // Multiple precursors
    {
      stringstream ss( dummyString.substr( 2, dummyString.length()-3 ) );
      ss >> precursor;
      precursors.push_back( precursor );

      while( *dummyString.rbegin() != '\"' )
      {
        lineStream >> precursor;
        precursors.push_back( precursor );
        lineStream >> dummyString; // Either , oder }"
      }
    }

    // Ignore the RGBColor
    lineStream >> dummyString >> dummyString >> dummyString;


    // Control output
//    if( lineNumber % 1000 == 0 )
//    {
//      cout << "Line " << lineNumber << ": " << id << " "
//           << x << " " << y << " " << z << " " << time << " " << rad << " { ";
//      for( list<int>::const_iterator it = precursors.begin();
//           it != precursors.end(); ++it ) cout << *it << " ";
//      cout << "}" << endl;
//    }

    // for the raw data of the Arabidopsis the rad variable holds a constant
    // value for the radius of cells (10) but in the case of the modelling
    // we use this column to store the area of the tissue of a cell for a
    // specific time step; also determine the min and max values
    if( time > cellAreas.size() )
      cellAreas.resize( time );

    cellAreas.at( time-1 ).insert( std::make_pair( id, log(rad) ) );
    if( log(rad) <= minMaxArea.first )
      minMaxArea.first = log(rad);

    if( log(rad) > minMaxArea.second )
      minMaxArea.second = log(rad);

    // store division type for each division node
    if( _data->divisionInfoIsIncluded() )
    {
      int divType = this->getDivisionType( line );

      // do always increase the size of the division type
      // data structure to ensure always a valid access
      // to the outer vector (=number of time steps)
      if( time > cellDivisionTypes.size() )
        cellDivisionTypes.resize( time );

      if( divType != -1 )
        cellDivisionTypes.at( time-1 ).insert( std::make_pair( id, divType ) );
    }

    /* Build the tree structure */

    if( id != prevId ) // We have a new "object", i.e. (sub)tree
    {
      if( _data->layerInfoIsIncluded() )
      {
        int layerValue = this->getSequenceLayerValue( line );

        // resize data types if required
        if( time > cellLayers.size() )
          cellLayers.resize( time );

        cellLayers.at(time-1).insert( std::make_pair( id, layerValue ) );
      }

      if( _data->trackGroupInfoIncluded() && precursors.empty() )
      {
        // get and check the track group column in order to get information
        // about the number of cells in each file and to know the master file
        // NOTE that this is only applied for the first entry
        // of a cell migration (always given as pair of lines)
        int cellFile = this->getCellFile( line );
        cellFiles.insert( std::make_pair( id, cellFile ) );
      }

      // Check (and correct if necessary) the last object
      if( prevNode != NULL ) // Only NULL for the first line
        correctObjectHistory( prevNode, lineNumber-1 );

      if( precursors.empty() ) // Insert root nodes into the collection
      {
        treeId = id;
        _data->insert( treeId, node );
      }
      else // Look for the right tree and the right path there
      {
        SLineageTree* root = (*_data)[precursors.front()];
        if( root == NULL )
          THROW_EXCEPTION( VInvalidFileException, "Unable to find root" );

        treeId = root->treeId;

        for( list<int>::const_iterator it = precursors.begin();
             it != precursors.end(); ++it )
        {
          root = root->findCell( *it );
          if( root == NULL )
            THROW_EXCEPTION( VInvalidFileException, "Unable to find subRoot" );
        }

        root->addSubtree( node ); // TODO Maybe search right inclusion place between root and root->subRoot()

        // Some hints in case of corrupt data
        if( root->children.size() > 2 )
        {
          ++tooManyChildrenCases;
          cout << "line " << lineNumber << ": The following node has "
               << root->children.size() << " children: " << "cell id: "
               << root->cellId << ", time step: " << root->timeStep << endl
               << "Just added this node: " << "cell id: " << id
               << ", time step: " << time << endl;
          std::cout << "children of " <<  root->cellId << " t = " << root->timeStep << " children: ";
          for( std::size_t i = 0; i < root->children.size(); ++i )
            std::cout << root->children[i]->cellId << " t = " << root->children[i]->timeStep << " ";
          std::cout << endl;
        }
      }
    }
    else // No division, just move
    {
      if( _data->layerInfoIsIncluded() )
      {
        int layerValue = this->getSequenceLayerValue( line );

        // resize data types if required
        if( time > cellLayers.size() )
          cellLayers.resize( time );

        for( std::size_t t = prevTime; t <= time; t++ )
          cellLayers.at(t-1).insert( std::make_pair( id, layerValue ) );
      }

      if( _data->trackGroupInfoIncluded() && precursors.empty() )
      {
        // probably not required
        int cellFile = this->getCellFile( line );
        cellFiles.insert( std::make_pair( id, cellFile ) );
      }

      // handle case of performing an interpolation for cell movements
      if( _interpolateMovements )
      {
        treeId = prevNode->treeId;
        SArray p1( prevX, prevY, prevZ*2. );
        SArray p2( x, y, z*2. );
        this->interpolateCellMovements( p1, p2,
                                        prevTime, time,
                                        prevNode, node,
                                        id, treeId );
      }
      else
      {
        treeId = prevNode->treeId;
        prevNode->addSubtree( node );
      }
    }

    // Set the renderable tree attributes
    node->setX( x );
    node->setY( y );
    node->setZ( z*2. );

    // Set the lineage tree attributes
    node->cellId = id;
    node->timeStep = time;
    node->treeId = treeId;

    // Update prev data
    prevId = id;
    prevNode = node;
    prevTime = time;

    prevX = x;
    prevY = y;
    prevZ = z;

    if( node->parent )
      sortChildren( node->parent );
  }

  in.close();

  if( _data->trackGroupInfoIncluded() )
  {
    // set master file id
    _data->setMasterFileId( 0 );

    // insert file, layer and division information into data
    _data->setCellFiles( cellFiles );
  }

  if( _data->layerInfoIsIncluded() )
    _data->setCellLayers( cellLayers );

  if( _data->divisionInfoIsIncluded() )
    _data->setDivisionType( cellDivisionTypes );

  // set information about tissue area
  _data->setAreaData( cellAreas );
  _data->setMinMaxArea( minMaxArea );

  SKernel::get()->setProgress( 1.0, "Lineages Built" );

  cout << "Accumulated cases of too many children: " << tooManyChildrenCases << endl;
  unsigned int problematicChildren = 0;
  for( SAdvancedLineageTreeCollection::const_iterator tree = _data->begin();
       tree != _data->end(); ++tree )
  {
    for( SLineageTree::const_iterator node = tree->second->begin();
         node != tree->second->end(); ++node )
    {
      if( node->children.size() > 2 ) ++problematicChildren;
    }
  }
  cout << "Cases of too many children: " << problematicChildren << endl;

}

//----------------------------------------------------------------------------//

void SLoadArabidopsisLineages::sortChildren( SLineageTree *root )
{
  if( root->children.size() == 2 )
  {
    osg::Vec3 p1( root->children[0]->getX(), root->children[0]->getY(),
                 root->children[0]->getZ() );
    osg::Vec3 p2( root->children[1]->getX(), root->children[1]->getY(),
                 root->children[1]->getZ() );

    if( _applyRotations )
    {
      SArray rotation = _data->getRotation();

      osg::MatrixTransform *trafo = new osg::MatrixTransform(
                       osg::Matrix::rotate( rotation[0]/180.*M_PI,
                                            osg::Vec3(1,0,0) ) );
      trafo->postMult( osg::Matrix::rotate( rotation[1]/180.*M_PI,
                                            osg::Vec3(0,1,0) ));
      trafo->postMult( osg::Matrix::rotate( rotation[2]/180.*M_PI,
                                            osg::Vec3(0,0,1)));

      p1 = p1 * trafo->getMatrix();
      p2 = p2 * trafo->getMatrix();
    }
    
    //std::cout << root->children[0]->cellId << " " << root->children[1]->cellId << " " << p1[0] << " " << p2[0] << std::endl;

    // switch children if in wrong order
    if( p1[0] > p2[0] )
    {
      SLineageTree *dummy = root->children[0];
      root->children[0] = root->children[1];
      root->children[1] = dummy;
    }
  }
}

//----------------------------------------------------------------------------

void SLoadArabidopsisLineages::checkDataContent( const std::string &line )
{
  // the optimal setting would be that the last three columns
  // include the following columns:
  // TrackGroup Layer DivisionType
  //
  // check which one of them exists and only allow access to existing data

  // initialize vector for the three last titles of the columns
  std::vector<std::string> entries;
  entries.resize(3, "");
  int positionIndex = 0;
  for( std::string::const_reverse_iterator sIter = line.rbegin();
       sIter != line.rend(); sIter++ )
  {
    if( *sIter == ' ' )
    {
      positionIndex++;
      sIter++;

      // break loop after the third word
      if( positionIndex == 3 )
        break;
    }

    // add the chars to the corresponding string content
    entries.at(positionIndex) += *sIter;
  }

  // swap content for all words
  for( std::size_t s = 0; s < entries.size(); s++ )
    entries.at(s) = std::string( entries.at(s).rbegin(), entries.at(s).rend() );

  for( std::size_t s = 0; s < entries.size(); s++ )
  {
    if( entries.at(s) == "Layer" )
    {
      _data->layerInfoIsIncluded( true );
      _data->layerColumnIndex( s );
      _indexLayer = s;
    }
  }

  for( std::size_t s = 0; s < entries.size(); s++ )
  {
    if( entries.at(s) == "DivisionType" )
    {
      _data->divisionInfoIsIncluded( true );
      _data->divisionColumIndex( s );
      _indexDivisionType = s;
    }
  }

  for( std::size_t s = 0; s < entries.size(); s++ )
  {
    if( entries.at(s) == "TrackGroup" )
    {
      _data->trackGroupInfoIncluded( true );
      _data->trackGroupColumnIndex( s );
      _indexTrackGroup = s;
    }
  }
}

//----------------------------------------------------------------------------

int SLoadArabidopsisLineages::getDivisionType( const std::string &line )
{
  std::string cellDivStr;

  int positionIndex = 0;
  for( std::string::const_reverse_iterator sIter = line.rbegin();
       sIter != line.rend(); sIter++ )
  {
    if( *sIter == ' ' )
    {
      positionIndex++;
      sIter++;

      if( positionIndex == _indexDivisionType+1 )
        break;
    }

    if( positionIndex == _indexDivisionType )
      cellDivStr += *sIter;
  }

  // swap content and cast to std::size_t
  cellDivStr = std::string( cellDivStr.rbegin(), cellDivStr.rend() );

  if( cellDivStr != "-" )
    return boost::lexical_cast<int>( cellDivStr );
  else
    return -1;
}

//----------------------------------------------------------------------------

int SLoadArabidopsisLineages::getCellFile( const std::string &line )
{
  std::string cellFileStr;

  int positionIndex = 0;
  for( std::string::const_reverse_iterator sIter = line.rbegin();
       sIter != line.rend(); sIter++ )
  {
    if( *sIter == ' ' )
    {
      positionIndex++;
      sIter++;

      if( positionIndex == _indexTrackGroup+1 )
        break;
    }

    if( positionIndex == _indexTrackGroup )
      cellFileStr += *sIter;
  }

  // swap content
  cellFileStr = std::string( cellFileStr.rbegin(), cellFileStr.rend() );
  return boost::lexical_cast<int>( cellFileStr );
}

//----------------------------------------------------------------------------

int SLoadArabidopsisLineages::getSequenceLayerValue( const std::string &line )
{
  std::string::const_reverse_iterator sIter = line.rbegin();
  for( ;sIter != line.rend(); sIter++ )
  {
    // start of reversed layer string
    if( *sIter == '}' )
      break;
  }

  sIter++;
  // store the layer sequence
  std::string layerSeq = "";

  while( *sIter != '{' )
  {
    layerSeq += *sIter;
    sIter++;
  }

  // swap content
  layerSeq = std::string( layerSeq.rbegin(), layerSeq.rend() );

  // compute layer value
  int layerValue = 1;
  for( std::string::iterator iter = layerSeq.begin();
       iter != layerSeq.end(); iter++ )
  {
    if( *iter == '1' )
      layerValue *= 2;

    if( *iter == '2' )
    {
      layerValue *= 2;
      layerValue += 1;
    }
  }

  //std::cout << "val: " << layerValue << std::endl;
  return layerValue;
}

//----------------------------------------------------------------------------

int SLoadArabidopsisLineages::getLayerValue( const std::string &line )
{
  std::string cellLayerStr;

  int positionIndex = 0;
  for( std::string::const_reverse_iterator sIter = line.rbegin();
       sIter != line.rend(); sIter++ )
  {
    if( *sIter == ' ' )
    {
      positionIndex++;
      sIter++;

      if( positionIndex == _indexLayer+1 )
        break;
    }

    if( positionIndex == _indexLayer )
      cellLayerStr += *sIter;
  }

  // swap content
  cellLayerStr = std::string( cellLayerStr.rbegin(), cellLayerStr.rend() );
  return boost::lexical_cast<int>( cellLayerStr );
}

//----------------------------------------------------------------------------

void SLoadArabidopsisLineages::interpolateCellMovements( const SArray &p1, const SArray &p2,
                                                    const int prevTime, const int time,
                                                    SLineageTree *prNode, SLineageTree *node,
                                                    const int id, const int treeId )
{
  int deltaTime = time - prevTime;
  int startTime;
  SArray p,q;

  // if we have no time difference at all which is
  // the case of the full segmented data set then
  // leave at once the method
  if( deltaTime == 0 )
    return;

  // if delta time == 1 then no interpolation
  // is needed and we add the only node as a subtree
  // and return
  if( deltaTime == 1 )
  {
    prNode->addSubtree( node );
    return;
  }

  // swap cell events if they occur (what they sometimes do in the data set)
  if( deltaTime < 0 )
  {
    deltaTime *= -1.;
    startTime = time + 1;
    p = p2;
    q = p1;
  }
  else
  {
    startTime = prevTime + 1;
    p = p1;
    q = p2;
  }

  int curTime = startTime;

  double steps = 1./deltaTime;

  SLineageTree *prevNode = NULL;

  double k = steps;
  for( int s = 0; s < deltaTime-1; s++, k += steps )
  {
    SLineageTree *node = new SLineageTree;

    double x = (1-k) * p[0] + k* q[0];
    double y = (1-k) * p[1] + k* q[1];
    double z = (1-k) * p[2] + k* q[2];

    node->setX( x );
    node->setY( y );
    node->setZ( z );
    node->cellId = id;
    node->timeStep = curTime;
    node->treeId = treeId;

//    std::cout << "info: " << treeId << " " << id
//              << " " << x << " " << y << " " << z
//              << " " << curTime << "s: " << s << std::endl;

    if( startTime == curTime )
    {
      prNode->addSubtree( node );
      prevNode = node;
    }
    else
    {
      prevNode->addSubtree( node );
      prevNode = node;
    }

    curTime++;
  }

  // at last append last node
  prevNode->addSubtree( node );
}

//----------------------------------------------------------------------------//

void SLoadArabidopsisLineages::correctObjectHistory( SLineageTree* bottom,
                                                unsigned int lnum )
{
  // The data is partially pretty messed up ...

  SLineageTree* up = bottom->parent;

  while( up != NULL && up->cellId == bottom->cellId )
  {
    bool duplicateTimeStep = false;
    SLineageTree* duplicateBackup = NULL;
    bool swapped = false;

    if( bottom->timeStep == up->timeStep )
    /* This case actually occurs in the data! There are a few lines which are
     * exactly equal to the one before. We just ignore such nodes here. */
    {
      duplicateTimeStep = true;
      duplicateBackup = bottom;
//      cout << "Line " << lnum << ": Two instances of one object at the same "
//           << "time! Ignoring this line" << endl;
    }

    if( bottom->timeStep < up->timeStep )
    /* This occurs quite often and in various different situations.
     * Generally we remove this node from the chain and put it back in between
     * up->parent and up. */
    {
//      cout << "Line " << lnum << ": time step is smaller than in "
//           << "the line before! But don't worry - it is being delt with."
//           << endl;

      if( up->parent == NULL ) // time steps like 2 1 3 4 5 ...
      {
        // Overwrite false root
        (_data->iExceptionallyWantToEditTheLineages())[bottom->treeId] = bottom;
        bottom->removeThisNodeFromTree();
        bottom->addSubtree( up );
        return;
      }
      else // time steps like 1 3 2 4 5 ...
      {
        // There is even the case 216-240 then 217-240 (same id!)
        // For the later group the following is true
        if( up->getPatriarch()->findNode( bottom->cellId, bottom->timeStep )
            != bottom )
        {
          duplicateTimeStep = true;
          duplicateBackup = bottom;
          cout << "-> Attention! The data is severely corrupted. Found "
               << "reocurring cell, will ignore it." << endl;
        }
        else
        {
          // Possibly more than one wrong time step directly after each other
          // NOTE this is all very inefficent if the whole object was given in
          // reverse order :-(
          swapped = true;
          bottom->removeThisNodeFromTree();
          up->parent->addSubtree( bottom );
          up->cropThisSubtree();
          bottom->addSubtree( up );
        }
      }
    }


    // Usually traverse bottom up
    if( !swapped )
    {
      bottom = bottom->parent;
      up = up->parent;
      --lnum;
    }
    else // Handling (multiple) swaps
    {
      if( !up->children.empty() )
      {
        // bottom is now up's parent, so we put it back down
//        for( int gap = 1; gap < 7; ++gap )
//        {
//          bottom = up->findNode( up->cellId, up->timeStep+gap );
//          if( bottom != NULL ) break;
//        }
//        if( bottom == NULL )
//        {
//          THROW_EXCEPTION( VInvalidOrderException,
//                           "Could not find adequate child node for traversal" );
//          return;
//        }
        bottom = up->children.back();
        if( bottom->parent != up )
        {
          THROW_EXCEPTION( VInvalidOrderException,
                           "Node is not its parent's child" );
          return;
        }
        ++lnum;
      }
      else
      {
        SLineageTree* newBottom = up;
        up = bottom;
        bottom = newBottom;
      }
    }

    if( duplicateTimeStep )
    {
      duplicateBackup->removeThisNodeFromTree();
      delete duplicateBackup;
    }
//    cout << "bottom cell id: " << bottom->cellId << " up cell id: " << up->cellId << endl;
  }

//  cout << "Finished correcting cell " << bottom->cellId << endl;
}
