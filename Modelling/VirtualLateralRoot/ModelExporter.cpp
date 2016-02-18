#include "ModelExporter.h"

/**
  @file   ModelExporter.cpp
  @brief  Namespace for exporting modelling results
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

namespace ModelExporter
{

//----------------------------------------------------------------

void exportLineageInformation( const std::string &filename,
                               const cell& c,
                               const MyTissue& T,
                               const bool init )
{
  if( init )
  {
    std::ofstream out( filename.c_str(), std::ofstream::out );
  
    // rotation information
    out << "0 0 0\n";
    // center of root
    out << "0 -2 0\n";
    // dimension
    out << "3 2 0\n";
    // header
    out << "ObjectID X Y Z Timepoint Radius Precursors Color Lineage TrackID TrackColor TrackGroup Layer DivisionType\n";
    out.close();
    return;
  }
 
  std::ofstream out( filename.c_str(), std::ofstream::out | std::ofstream::app );
  
  out << c->id << " "
      << c->center.i() << " "
      << c->center.j() << " "
      << c->center.k() << " "
      << c->timeStep << " "
      << c->area << " "
      << "\"{"; 
  
  std::size_t counter = 1;
  for( std::set<std::size_t>::iterator setIter = c->precursors.begin();
        setIter != c->precursors.end(); setIter++, counter++ )
  {
    if( counter != c->precursors.size() )
      out << *setIter << ", ";
    else
      out << *setIter << "}\" ";
  }
  
  if( c->precursors.size() == 0 )
    out << "}\" ";
      
  out << "Color "
      << c->treeId << " "
      << "TrackId "
      << "TrackColor "
      << c->cellFile << " "
      << c->layerValue << " ";
      
  if( c->divType == DivisionType::ANTICLINAL )
    out << "0 \n";
  else if( c->divType == DivisionType::PERICLINAL )
    out << "1 \n";
  else
    out << "2 \n";
  
  out.close();
}

// ---------------------------------------------------------------------

void exportTimeAgainstCells( const std::string &filename,
                             const double dT,
                             const std::size_t numCells,
                             const bool init )
{
  std::ofstream out;
  
  // header
  if( init )
  {
    out.open( filename.c_str(), std::ofstream::out );
    out << "dT Cells\n";
  }
  else
  {
    out.open( filename.c_str(), std::ofstream::out | std::ofstream::app );
    out << dT << " " << numCells << "\n";
  }
  
  out.close();
}

// ---------------------------------------------------------------------

void exportCellWalls( const std::string &filename,
                      const cell& c,
                      const MyTissue &T,
                      const bool init )
{
  // header
  if( init )
  {
    std::ofstream out( filename.c_str(), std::ofstream::out );
    out << "ID Time CellCorners\n";
    out.close();
    return;
  }
  
  std::ofstream out( filename.c_str(), std::ofstream::out | std::ofstream::app );
  out << c->id << " " << c->timeStep << " ";
  bool first = true;
  forall(const junction& j, T.S.neighbors(c))
  {
    if( first )
      first = false;
    else
      out << ",";
  
    out << j->getPos().i() << "," 
        << j->getPos().j() << "," 
        << j->getPos().k();
  }
  
  out << "\n";
  
  out.close();
}

// ---------------------------------------------------------------------

void exportDivisionAngles( const std::string &filename,
                           const cell& c,
                           const bool init )
{
  // header
  if( init )
  {
    std::ofstream out( filename.c_str(), std::ofstream::out );
    out << "CellID TimeStep CellNumber File DivisionType Theta\n";
    out.close();
    return;
  }
  
  std::ofstream out( filename.c_str(), std::ofstream::out | std::ofstream::app );
  out << c->id << " "
      << c->timeStep << " "
      << c->cellNumber << " "
      << c->cellFile << " ";

  if( c->divType == DivisionType::ANTICLINAL )
    out << "0 ";
  else if( c->divType == DivisionType::PERICLINAL )
    out << "1 ";
  else
    out << "2 ";
  
  out << c->angle << "\n";
  
  out.close();
}

// ---------------------------------------------------------------------

void exportNumberOfCellsInFiles( const std::string &filename,
                                 const MyTissue &T,
                                 const int run )
{ 
  std::vector<std::size_t> numCells;
  numCells.resize(5);
  
  forall(const cell& c, T.C)
    numCells.at( c->cellFile+2 )++;
  
  std::ofstream out( filename.c_str(), std::ofstream::out | std::ofstream::app );
  out << run << " ";
  
  for( std::size_t f=0; f < numCells.size(); f++ )
    out << numCells.at(f) << " ";
  
  out << "\n";
  
  out.close(); 
}

// ---------------------------------------------------------------------

void exportDivisionDaughterProperties( const std::string &filename,
                                       const std::vector<cell> &dC,
                                       const DivisionType::type divType,
                                       const double angleThreshold,
                                       std::pair<std::size_t, std::size_t> &divOccurrences,
                                       const bool init )
{
  if( init )
  {
    std::ofstream out( filename.c_str(), std::ofstream::out );
    // header
    out << "'Cell ID' 'Parent ID' 'Lineage ID' 'Timestep' 'Cell Cycle' " //'Area' 'LongestCellWall'
        << "'Dir of last Division X' 'Dir of last Division Y' 'Dir of last Division Z' "
        << "'Dir of this Division X' 'Dir of this Division Y' 'Dir of this Division Z' "
        << "'Angle' 'Cell File' 'Cell File Sequence' 'Cell Layer' 'Division Type' 'Division Sequence' 'XPos' 'YPos' 'ZPos'\n";
    out.close();
    return;
  }
  
  std::ofstream out( filename.c_str(), std::ofstream::out | std::ofstream::app );
  
  for( std::size_t i = 0; i < dC.size(); i++ )
  {
    cell c = dC.at(i);
    
    out << c->id << " " << c->parentId << " "
        << c->treeId << " " << c->timeStep << " "
        << c->cellCycle << " " << c->previousDivDir.i() << " "
        << c->previousDivDir.j() << " " << c->previousDivDir.k() << " "
        << c->divDir.i() << " " << c->divDir.j() << " "
        << c->divDir.k() << " ";
      
    if( c->cellCycle != 1 )
    {
      Point3d lastDir = c->previousDivDir;
      Point3d curDir = c->divDir;
      
      lastDir.normalize();
      curDir.normalize();
      
      out << acos( lastDir*curDir )*180./M_PI << " ";
    }
    else
      out << "NA ";
    
    out << c->cellFile << " \"" << c->cellFileSequence << "\" " << c->layerValue << " ";
    
    if( divType == DivisionType::ANTICLINAL ||
        divType == DivisionType::RADIAL )
    {
      if( i == 0 )
        divOccurrences.first++;
      
      if( divType == DivisionType::ANTICLINAL )
        out << "0 ";
      else
        out << "2 ";
    }
    else
    {
      if( i == 0 )
        divOccurrences.second++;
      
      out << "1 ";
    }
    
    out << "\"" << c->divisionLetterSequence << "\" ";
    
    // at last position of nucleus
    out << c->center.i() << " " << c->center.j() << " " << c->center.k() << "\n";
  }
  
  out.close();
}

// ---------------------------------------------------------------------

void exportModelProperties( const std::string &filename,
                            const std::size_t loopCounter,
                            const layerMap &firstLayerAppearances,
                            const std::map<std::string, std::size_t> &totalLayerCount,
                            const std::pair<std::size_t, std::size_t> &divOccurrences,
                            const std::size_t numCells,
                            const std::size_t numSamples,
                            const bool init )
{
  // generate the division sequences automatically
  std::size_t maxLength = 6;
  std::vector<std::string> seqVecTemp;
  std::vector<std::string> totalSeqVec;
  seqVecTemp.push_back( "0" );
  totalSeqVec.push_back( "0" );
  for( std::size_t l = 1; l < maxLength; l++ )
  {
    // double entries of seqVecTemp
    std::size_t entries = seqVecTemp.size();
    for( std::size_t m = 0; m < entries; m++ )
      seqVecTemp.push_back( seqVecTemp.at(m) );
    
    std::size_t newL = std::pow(2,l);
    for( std::size_t s = 0; s < newL; s++ )
    {
      if( s < (std::size_t)(newL/2.) )
        seqVecTemp.at(s).insert( 1, "1" );
      else
        seqVecTemp.at(s).insert( 1, "2" );
      
      totalSeqVec.push_back( seqVecTemp.at(s) );
    }
  }
  
  std::ofstream out;
  if( init )
  {
    out.open( filename.c_str(), std::ofstream::out );
    // header
    out << numSamples << " " << maxLength << "\n";
    out << "'Sample' 'Cells' 'Anticlinal' 'Periclinal' ";
      
    // export the sequences twice, one for the time steps
    // and the other one for the counts
    for( std::size_t i = 0; i < 2; i++ )
      for( std::size_t l = 0; l < totalSeqVec.size(); l++ )
        out << "'" << totalSeqVec.at(l) << "' ";

    out << "\n";
  }
  else
    out.open( filename.c_str(), std::ofstream::out | std::ofstream::app );
  
  out << loopCounter << " ";
  out << numCells << " ";
  out << divOccurrences.first << " ";
  out << divOccurrences.second << " ";
  
  // first export the time steps
  for( std::size_t l = 0; l < totalSeqVec.size(); l++ )
  {
    bool found = false;
    for( auto iter = firstLayerAppearances.begin();
         iter != firstLayerAppearances.end(); ++iter )
    {
      if( iter->second.second.compare( totalSeqVec.at(l) ) == 0 )
      {
        // write the time step of the corresponding sequence
        out << iter->second.first << " ";
        found = true;
        break;
      }
    }
    
    if( !found )
      out << "-1 ";
  }
  
  // then export the counts of each sequence
  for( std::size_t l = 0; l < totalSeqVec.size(); l++ )
  {
    bool found = false;
    for( auto iter = totalLayerCount.begin();
         iter != totalLayerCount.end(); ++iter )
    {
      if( iter->first.compare( totalSeqVec.at(l) ) == 0 )
      {
        // write the time step of the corresponding sequence
        out << iter->second << " ";
        found = true;
        break;
      }
    }
    
    if( !found )
      out << "-1 ";
  }
  
  out << "\n";
  out.close();
}

// ---------------------------------------------------------------------

void exportDivisionSequences( const std::string &filename,
                              const std::map< std::string, std::vector<std::size_t> > &divisionSequences )
{
  std::ofstream out;
  out.open( filename.c_str(), std::ofstream::out );
  
  for( auto iter = divisionSequences.begin();
       iter != divisionSequences.end(); iter++ )
  {
    out << iter->first << " ";
    for( std::size_t f = 0; f < iter->second.size(); f++ )
      out << iter->second.at(f) << " ";
    
    out << "\n";
  }
  
  out.close();
}

// ---------------------------------------------------------------------

void exportPropabilityDistribution( const std::string &filename,
                                    const std::vector<std::vector<double> > &probValues,
                                    const std::vector<std::vector<double> > &lengths,
                                    const std::vector<std::size_t> choice )
{ 
  std::ofstream out;
  out.open( filename.c_str(), std::ofstream::out | std::ofstream::app );
  
  if( probValues.size() == 0 || lengths.size() == 0 )
    return;
  
  // for each length and probability value
  for( std::size_t l=0; l<probValues.at(0).size(); l++ )
  {
    // for each division event
    for( std::size_t d=0; d<probValues.size(); d++ )
    {
      out << lengths.at(d).at(l) << " "
          << probValues.at(d).at(l) << " ";
    }
    out << "\n";
  }
  
  for( std::size_t c=0; c<choice.size(); c++ )
    out << choice.at(c) << " " << "-" << " ";
    
  out << "\n";
  out.close();
}

// ---------------------------------------------------------------------

}