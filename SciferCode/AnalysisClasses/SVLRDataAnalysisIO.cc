#include "SVLRDataAnalysisIO.hh"

#include "SSimilarityMeasureUtil.hh"

#include <fstream>

namespace SVLRDataAnalysisIO
{

// ---------------------------------------------------------------------

void printDivisionInformation( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                               const std::size_t numTotalTimesteps,
                               const std::vector<osg::Vec3> &centresOfMass )
{
  std::pair<std::size_t, std::size_t> range;
  range.first = 18;
  range.second = 143;
  std::size_t epsilon = 1;
  std::size_t numDivisions = 0;
  std::size_t numAnticlinalDivisions = 0;
  std::size_t numPericlinalDivisions = 0;
  std::size_t numRadialDivisions = 0;
  std::size_t numSyncAnticlinalDivisions = 0;
  std::size_t numSyncPericlinalDivisions = 0;
  std::size_t numSyncRadialDivisions = 0;
  std::size_t numCellsInMasterFileStart = 0;
  std::size_t numCellsInMasterFileEnd = 0;

  std::size_t numDivisionsInMaster = 0;
  std::size_t numAnticlinalDivisionsInMaster = 0;
  std::size_t numPericlinalDivisionsInMaster = 0;
  std::size_t numRadialDivisionsInMaster = 0;
  std::size_t numSyncAnticlinalDivisionsInMaster = 0;
  std::size_t numSyncPericlinalDivisionsInMaster = 0;
  std::size_t numSyncRadialDivisionsInMaster = 0;
  std::size_t maxCells = 0;

  SArray rotation = lineages->getRotation();
  osg::Matrix rotMat;
  rotMat.postMult( osg::Matrix::rotate( rotation[2]/180.*M_PI, osg::Vec3(0,0,1) ) );
  rotMat.postMult( osg::Matrix::rotate( rotation[1]/180.*M_PI, osg::Vec3(0,1,0) ) );
  rotMat.postMult( osg::Matrix::rotate( rotation[0]/180.*M_PI, osg::Vec3(1,0,0) ) );

  osg::Vec3 lowBound( -60., -300., -300. );
  osg::Vec3 upBound( 60., 300., 300. );
  std::size_t numCellsInBoundary = 0;
  std::size_t registeredTimeStep = 0;

  if( lineages->getName() == "120830" )
    registeredTimeStep = 269;
  else if( lineages->getName() == "121204" )
    registeredTimeStep = 277;
  else if( lineages->getName() == "121211" )
    registeredTimeStep = 230;
  else if( lineages->getName() == "130508" )
    registeredTimeStep = 344;
  else if( lineages->getName() == "130607" )
    registeredTimeStep = 213;

  // first determine the number of cells for each time step
  std::vector<std::size_t> numCellsPerTimeStep;
  numCellsPerTimeStep.resize( numTotalTimesteps, 0 );
  for( auto l = lineages->begin(); l != lineages->end(); ++l )
  {
    SLineageTree *tree = l->second;
    for( auto nodeIt = tree->begin(); nodeIt != tree->end(); ++nodeIt )
      numCellsPerTimeStep.at(nodeIt->timeStep-1)++;
  }

  NodeFeatureInfo divisionScheme = lineages->getDivisionType();
  std::map<int,int> cellFiles = lineages->getCellFiles();

  for( SAdvancedLineageTreeCollection::const_iterator l = lineages->begin();
       l != lineages->end(); ++l )
  {
    SLineageTree *tree = l->second;

    int cellFile = cellFiles.at( tree->treeId );
    bool master = cellFile == 0;

    // store number of cells at the last time step for all trees located in the master cell file
    if( master )
    {
      numCellsInMasterFileStart += 1;
      numCellsInMasterFileEnd += tree->numDivisions()+1;
    }

    maxCells += tree->nodes();

    for( SLineageTree::const_iterator treeIter = tree->begin();
         treeIter != tree->end(); ++treeIter )
    {
      int cellId = treeIter->cellId;
      std::size_t timeStep = treeIter->timeStep;
      std::size_t numCells = numCellsPerTimeStep.at(timeStep-1);

      // if the current time steps is the registered time step
      // check if the cell lies within the boundary zone
      if( timeStep == registeredTimeStep )
      {
        osg::Vec3 cen = centresOfMass.at( timeStep-1 );
        cen = cen * rotMat;
        osg::Vec3 pos( treeIter->getX(), treeIter->getY(), treeIter->getZ() );
        pos = pos * rotMat;
        pos = pos - cen;

        if( pos[0] >= lowBound[0] && pos[0] < upBound[0] &&
            pos[1] >= lowBound[1] && pos[1] < upBound[1] &&
            pos[2] >= lowBound[2] && pos[2] < upBound[2] )
          numCellsInBoundary++;
      }

      if( treeIter->children.size() == 2 )
      {
        numDivisions++;

        if( master )
          numDivisionsInMaster++;

        int divType = divisionScheme.at(timeStep-1).at(cellId);

        switch(divType)
        {
        case 0:
          numAnticlinalDivisions++;
          if( master )
            numAnticlinalDivisionsInMaster++;

          if( (numCells >= range.first - epsilon) &&
              (numCells < range.second + epsilon) )
          {
            numSyncAnticlinalDivisions++;
            if( master )
              numSyncAnticlinalDivisionsInMaster++;
          }
          break;
        case 1:
          numPericlinalDivisions++;
          if( master )
            numPericlinalDivisionsInMaster++;

          if( (numCells >= range.first - epsilon) &&
              (numCells < range.second + epsilon) )
          {
            numSyncPericlinalDivisions++;
            if( master )
              numSyncPericlinalDivisionsInMaster++;
          }
          break;
        case 2:
          numRadialDivisions++;
          if( master )
            numRadialDivisionsInMaster++;

          if( (numCells >= range.first - epsilon) &&
              (numCells < range.second + epsilon) )
          {
            numSyncRadialDivisions++;
            if( master )
              numSyncRadialDivisionsInMaster++;
          }
          break;
        }
      }
    }
  }

  std::cout << std::endl;
  std::cout << "=============================================" << std::endl;
  std::cout << "Division Properties for data: " << lineages->getName() << std::endl;
  std::cout << "=============================================" << std::endl;
  std::cout << "Max Divisions (in master): " << numDivisions << " ("
            << numDivisionsInMaster << ")" << std::endl;
  std::cout << "Anticlinal Divisions (in master): " << numAnticlinalDivisions << " ("
            << numAnticlinalDivisionsInMaster << ")" << std::endl;
  std::cout << "Periclinal Divisions (in master): " << numPericlinalDivisions << " ("
            << numPericlinalDivisionsInMaster << ")" << std::endl;
  std::cout << "Radial Divisions (in master): " << numRadialDivisions << " ("
            << numRadialDivisionsInMaster << ")" << std::endl;
  std::cout << "Number of cells in master cell file at start: " << numCellsInMasterFileStart << std::endl;
  std::cout << "Number of cells in master cell file at end: " << numCellsInMasterFileEnd << std::endl;
  std::cout << "=============================================" << std::endl;
  std::cout << "Synchronized cell range: " << range.first << ", " << range.second << std::endl;
  //std::cout << "Synchronized cell range in master: " << masterRange.first << ", " << masterRange.second << std::endl;
  std::cout << "Anticlinal Divisions in range (in master): " << numSyncAnticlinalDivisions << " ("
            << numSyncAnticlinalDivisionsInMaster << ")" << std::endl;
  std::cout << "Periclinal Divisions in range (in master): " << numSyncPericlinalDivisions << " ("
            << numSyncPericlinalDivisionsInMaster << ")" << std::endl;
  std::cout << "Radial Divisions in range (in master): " << numSyncRadialDivisions << " ("
            << numSyncRadialDivisionsInMaster << ")" << std::endl;
  //std::cout << "Number of cells in master cell file at synced start: " << numSyncCellsInMasterFileStart << std::endl;
  //std::cout << "Number of cells in master cell file at synced end: " << numSyncCellsInMasterFileEnd << std::endl;
  std::cout << "Max number of cells over all time steps: " << maxCells << std::endl;
  std::cout << "Number of cells in the boundary volume: " << numCellsInBoundary << std::endl;
  std::cout << "=============================================" << std::endl;
  std::cout << std::endl;
}

// ---------------------------------------------------------------------

void exportCompressedDataset( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                              const NodeFeatureInfo &layerValues )
{
  std::string origData = "/home/necrolyte/Uni/LateralRootGrowth/170815/";
  origData += lineages->getName() + ".csv";
  std::string line;
  std::ifstream in( origData.c_str() );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to read input file " + origData );
    return;
  }

  std::string fileName = "/tmp/" + lineages->getName() + ".csv";

  // open file at the end of the stream and allow writing
  std::ofstream out( fileName.c_str(), std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to write into data " + fileName );
    return;
  }

  SArray rotation = lineages->getRotation();
  SArray center = lineages->getCenter();
  SArray dimension = lineages->getDimension();

  out << "DataRotation: " << rotation[0] << " " << rotation[1] << " " << rotation[2] << "\n";
  out << "RootCenter: " << center[0] << " " << center[1] << " " << center[2] << "\n";
  out << "DataDimension: " << dimension[0] << " " << dimension[1] << " " << dimension[2] << "\n";
  out << "Id X Y Z T Precursors Lineage CellFile Layer DivType\n";

  for( std::size_t i=0;i<4;i++ )
    getline( in, line );

  while( in.good() ) // Now read the file line per line
  {
    std::string dummyString, divType, layerStr;
    std::size_t time, id, lineage;
    std::string precursor = "";
    double x, y, z;
    int cellFile;

    getline( in, line );
    if( line == "" )
      break;

    std::stringstream lineStream( line );
    lineStream >> id >> x >> y >> z >> time >> dummyString >> dummyString;

    // multiple precursors
    if( dummyString.at(0) == '\"' )
    {
      precursor += dummyString;
      while( *dummyString.rbegin() != '\"' )
      {
        lineStream >> dummyString;
        precursor += dummyString;
      }
    }
    else
      precursor += dummyString;

    lineStream >> dummyString >> dummyString >> dummyString
               >> lineage;

    lineStream >> dummyString;

    // ignore multiple trackgroup stuff
    if( dummyString.at(0) == '\"' )
    {
      while( *dummyString.rbegin() != '\"' )
        lineStream >> dummyString;
    }

    lineStream >> dummyString >> dummyString >> dummyString;

    lineStream >> cellFile >> layerStr;

    if( layerStr.at(0) == '\"' )
    {
      while( *layerStr.rbegin() != '\"' )
        lineStream >> layerStr;
    }

    lineStream >> divType;

    // get new layering in which both daughter cell were updated
    int newLayer = layerValues.at(time-1).at(id);

    // write new data
    out << id << " " << x << " " << y << " " << z << " " << time << " "
        << precursor << " " << lineage << " "
        << cellFile << " " << newLayer << " " << divType << "\n";
  }

  in.close();
  out.close();
}

// ---------------------------------------------------------------------


void exportCellDivisions( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                          const std::string &fileName )
{
  // open file at the end of the stream and allow writing
  std::ofstream out( fileName.c_str(), std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to write into data "
                     + fileName );
    return;
  }

  out << "Id Parent X Y Z DivX DivY DivZ DivType Time DivisionAngle\n";

  NodeFeatureInfo divisionScheme = lineages->getDivisionType();

  for( SAdvancedLineageTreeCollection::const_iterator l = lineages->begin();
       l != lineages->end(); ++l )
  {
    // we create a new copy of this lineage tree in order
    // to change its layout which should not happen
    // to the original data since for them we store the
    // pointers to identify them
    SLineageTree *tree = l->second;

    for( SLineageTree::const_iterator treeIter = tree->begin();
         treeIter != tree->end(); ++treeIter )
    {
      // write root properties
      if( !treeIter->parent )
      {
        int cellId = treeIter->cellId;
        std::size_t timeStep = treeIter->timeStep;

        out << cellId << " NA "
            << treeIter->getX() << " " << treeIter->getY() << " "
            << treeIter->getZ() << " NA NA NA NA " << timeStep << " NA\n";
      }
      // only process division nodes
      if( treeIter->children.size() == 2 )
      {
        int cellId = treeIter->cellId;
        std::size_t timeStep = treeIter->timeStep;

        int parentId;
        if( !treeIter->getLastDivision() )
          parentId = cellId;
        else
          parentId = treeIter->getLastDivision()->cellId;

        osg::Vec3 posC1( treeIter->children[0]->getX(),
            treeIter->children[0]->getY(), treeIter->children[0]->getZ() );

        osg::Vec3 posC2( treeIter->children[1]->getX(),
            treeIter->children[1]->getY(), treeIter->children[1]->getZ() );

        osg::Vec3 divDir = posC2 - posC1;
        divDir.normalize();

        int divType = divisionScheme.at(timeStep-1).at(cellId);

        out << cellId << " " << parentId << " "
            << treeIter->getX() << " " << treeIter->getY() << " "
            << treeIter->getZ() << " " << divDir[0] << " "
            << divDir[1] << " " << divDir[2] << " "
            << divType << " " << timeStep << " NA\n";

        for( std::size_t c=0; c < 2; c++ )
        {
          out << treeIter->children[c]->cellId << " " << cellId << " "
              << treeIter->children[c]->getX() << " "
              << treeIter->children[c]->getY() << " "
              << treeIter->children[c]->getZ() << " "
              << "NA NA NA NA " << treeIter->children[c]->timeStep << " ";
          double angle = SSimilarityMeasureUtil::computeDivisionAngle( treeIter->children[c] );
          if( angle != -1 )
            out << angle << "\n";
          else
            out << "NA\n";
        }
      }
      else if( treeIter->children.size() == 0 )
      {
        int cellId = treeIter->cellId;
        std::size_t timeStep = treeIter->timeStep;

        int parentId;
        if( !treeIter->getLastDivision() )
          parentId = cellId;
        else
          parentId = treeIter->getLastDivision()->cellId;

        out << cellId << " " << parentId << " "
            << treeIter->getX() << " " << treeIter->getY() << " "
            << treeIter->getZ() << " "
            << "NA NA NA NA " << timeStep << " NA\n";
      }
    }
  }

  out.close();
}

// ---------------------------------------------------------------------

void appendDivisionProperties( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                               const osg::Matrix &rotMat )
{
  std::string fileName = "/home/necrolyte/Uni/LateralRootGrowth/261115/DivisionOrientations_LateTimeStep_29102015_";
  fileName += lineages->getName();
  fileName += ".csv";
  std::ifstream in( fileName.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data " + fileName );
    return;
  }

  std::map< std::pair<std::size_t, std::size_t>, Divs > divs;

  // open file at the end of the stream and allow writing
  std::string outFileName = "/tmp/DivisionAnalysis_" + lineages->getName() + ".csv";
  std::ofstream out( outFileName.c_str(), std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to write into data "
                     + outFileName );
    return;
  }

  bool firstDivTime = false;

  for( SAdvancedLineageTreeCollection::const_iterator l = lineages->begin();
       l != lineages->end(); ++l )
  {
    for( SLineageTree::const_iterator iter = l->second->begin();
         iter != l->second->end(); ++iter )
    {
      if( iter->children.size() == 2 )
      {
        Divs div;
        std::size_t id = iter->cellId;
        std::size_t t = iter->timeStep;

        osg::Vec3 divPos( iter->getX(), iter->getY(), iter->getZ() );
        div.divPosOrig = divPos;
        divPos = divPos * rotMat;
        div.divPos = divPos;

        osg::Vec3 c1, c2;
        if( !firstDivTime )
        {
          // first child
          SLineageTree::const_iterator iter1 = iter->children[0];
          while( iter1->children.size() != 2 )
          {
            if( iter1->children.size() == 0 )
              break;
            else
              ++iter1;
          }

          c1 = osg::Vec3( iter1->getX(), iter1->getY(), iter1->getZ() );

          // second child
          SLineageTree::const_iterator iter2 = iter->children[1];
          while( iter2->children.size() != 2 )
          {
            if( iter2->children.size() == 0 )
              break;
            else
              ++iter2;
          }

          c2 = osg::Vec3( iter2->getX(), iter2->getY(), iter2->getZ() );
        }
        else
        {
          c1 = osg::Vec3( iter->children[0]->getX(),
              iter->children[0]->getY(),
              iter->children[0]->getZ() );

          c2 = osg::Vec3( iter->children[1]->getX(),
              iter->children[1]->getY(),
              iter->children[1]->getZ() );
        }

        div.dauPos1Orig = c1;
        c1 = c1 * rotMat;
        div.dauPos1 = c1;
        div.dauPos2Orig = c2;
        c2 = c2 * rotMat;
        div.dauPos2 = c2;

        divs.insert( std::make_pair( std::make_pair( id, t ), div ) );
      }
    }
  }

  std::string line;
  getline( in, line );
  line += " DivPosXPOrig DivPosYOrig DivPosZOrig DivPosX DivPosY DivPosZ";
  line += " DauPos1XOrig DauPos1YOrig DauPos1ZOrig";
  line += " DauPos2XOrig DauPos2YOrig DauPos2ZOrig";
  out << line << "\n";

  while( in.good() )
  {
    getline( in, line );
    if( line == "" )
      break;

    std::string lineC = line;
    std::stringstream lineSC( lineC );

    std::size_t id, t;
    lineSC >> id >> t;

    Divs div = divs.at( std::make_pair( id, t ) );

    line += " " + std::to_string( div.divPosOrig[0] ) + " "
            + std::to_string( div.divPosOrig[1] ) + " "
            + std::to_string( div.divPosOrig[2] );
    line += " " + std::to_string( div.divPos[0] ) + " "
        + std::to_string( div.divPos[1] ) + " "
        + std::to_string( div.divPos[2] );
    line += " " + std::to_string( div.dauPos1Orig[0] ) + " "
        + std::to_string( div.dauPos1Orig[1] ) + " "
        + std::to_string( div.dauPos1Orig[2] );
    line += " " + std::to_string( div.dauPos2Orig[0] ) + " "
        + std::to_string( div.dauPos2Orig[1] ) + " "
        + std::to_string( div.dauPos2Orig[2] ) + "\n";

    out << line;
  }

  in.close();
  out.close();
}

// ---------------------------------------------------------------------

void readAngleData( const std::string &fileName,
                    std::map<std::size_t, FileAngles> &angles )
{
  std::ifstream in( fileName.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data " + fileName );
    return;
  }

  std::string line;
  getline( in, line );

  while( in.good() )
  {
    getline( in, line );
    if( line == "" )
      break;

    std::string dummy;
    std::size_t id, t, numCells, divType;
    int cellFile;
    double theta, phi;
    std::stringstream lineStream( line );
    lineStream >> id >> t;

    for( std::size_t e = 0; e < 16; e++ )
      lineStream >> dummy;

    // read cell number
    lineStream >> numCells >> dummy >> cellFile >> dummy;

    if( dummy.at(0) == '\"' )
    {
      while( *dummy.rbegin() != '\"' )
        lineStream >> dummy;
    }

    lineStream >> divType >> theta >> phi;

//    std::size_t divType = 0;
//    if( phi > 60. )
//    {
//      if( theta > 60. )
//        divType = 2;
//      else if( theta < 30. )
//        divType = 0;
////      else
////        std::cout << "Type not set" << std::endl;
//    }
//    else if( phi < 30. && theta > 60. )
//      divType = 1;
////    else
////      std::cout << "Type not set" << std::endl;

    FileAngles fA;
    fA.cellFile = cellFile;
    fA.theta = theta;
    fA.phi = phi;
    fA.divType = divType;

    angles.insert( std::make_pair( numCells, fA ) );
  }

  in.close();
}

// ---------------------------------------------------------------------

void readModelAngleData( const std::string &fileName,
                         std::vector< std::map<std::size_t, FileAngles> > &angles )
{
  std::ifstream in( fileName.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data " + fileName );
    return;
  }

  std::string line;
  getline( in, line );

  std::map<std::size_t, FileAngles> Temp;
  std::size_t lastT = 1;

  while( in.good() )
  {
    getline( in, line );
    if( line == "" )
      break;

    std::size_t id, t, numCells, divType;
    int cellFile;
    double theta;
    std::stringstream lineStream( line );
    lineStream >> id >> t >> numCells >> cellFile >> divType >> theta;

    // new model run
    if( t < lastT )
    {
      angles.push_back( Temp );
      Temp.clear();
    }

    FileAngles fA;
    fA.cellFile = cellFile;
    fA.theta = theta;
    fA.phi = 0.;

    Temp.insert( std::make_pair( numCells, fA ) );
    lastT = t;
  }

  in.close();
}

// ---------------------------------------------------------------------

void exportAngleAnalysis( const std::string &fileName,
                          const std::vector< std::map<int, double> > &averagedValues )
{
  std::ofstream out( fileName.c_str(), std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to write into data "
                     + fileName );
    return;
  }

  out << "NumCells -2 -1 0 1 2\n";

  std::size_t numCells = 0;
  for( auto iter = averagedValues.begin(); iter != averagedValues.end(); ++iter, numCells++ )
  {
    if( iter->size() == 0 )
      continue;

    std::vector<double> values;
    values.resize(5);

    for( auto fileIter = iter->begin(); fileIter != iter->end(); ++fileIter )
      values.at( fileIter->first + 2 ) = fileIter->second;

    out << numCells << " ";
    for( std::size_t v=0; v < values.size(); v++ )
      out << values.at(v) << " ";

    out << "\n";
  }

  out.close();
}

// ---------------------------------------------------------------------

void readCellShapeVertices( const std::string &fileName,
                            std::vector< std::vector<osg::Vec3> > &vertices )
{
  std::ifstream in( fileName.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data " + fileName );
    return;
  }

  std::string line;
  unsigned int cellCounter = 0;

  while( in.good() )
  {
    getline( in, line );
    if( line == "" )
      break;

    unsigned int numVertices;
    std::stringstream lineStr( line );
    lineStr >> numVertices;
    //std::cout << "numVert: " << numVertices << std::endl;
    double x,y,z;
    std::vector<osg::Vec3> vert;
    for( std::size_t v=0; v < numVertices; v++ )
    {
      getline( in, line );
      std::stringstream lineStr( line );
      lineStr >> x >> y >> z;
      //std::cout << "pos: " << x << " " << y << " " << z << std::endl;
      vert.push_back( osg::Vec3( x, y, z ) );
    }
    vertices.push_back( vert );

    cellCounter++;
  }

  in.close();
}

// ---------------------------------------------------------------------

}
