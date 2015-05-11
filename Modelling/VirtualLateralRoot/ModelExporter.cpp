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
    out << "ObjectID X Y Z Timepoint Radius Precursors Color Lineage TrackID TrackColor TrackGroup Layer\n";
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
      << "0 "
      << c->layerValue << "\n";
  
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

void exportDivisionDaughterProperties( const std::string &filename,
                                       const cell& cl,
                                       const cell& cr,
                                       const MyTissue::division_data& ddata,
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
        << "'Angle' 'Cell File' 'Cell Layer' 'Division Type' 'XPos' 'YPos' 'ZPos'\n";
    out.close();
    return;
  }
  
  std::ofstream out( filename.c_str(), std::ofstream::out | std::ofstream::app );
  // left daughter cell
  out << cl->id << " "
      << cl->parentId << " "
      << cl->treeId << " " 
      << cl->timeStep << " "
      << cl->cellCycle << " "
      //<< cl->area << " "
      //<< cl->longestWallLength << " "
      << cl->previousDivDir.i() << " "
      << cl->previousDivDir.j() << " "
      << cl->previousDivDir.k() << " "
      << cl->divDir.i() << " "
      << cl->divDir.j() << " "
      << cl->divDir.k() << " ";
      
  if( cl->cellCycle != 1 )
  {
    Point3d lastDir = cl->previousDivDir;
    Point3d curDir = cl->divDir;
    
    lastDir.normalize();
    curDir.normalize();
    
    out << acos( lastDir*curDir )*180./M_PI << " ";
  }
  else
    out << "NA ";
  
  // master cell file
  out << "0 "
      << cl->layerValue << " ";
  
  DivisionType::type divType = ModelUtils::determineDivisionType( ddata, angleThreshold );
  if( divType == DivisionType::ANTICLINAL )
  {
    divOccurrences.first++;
    out << "0 ";
  }
  else
  {
    divOccurrences.second++;
    out << "1 ";
  }
  
  // at last position of nucleus
  out << cl->center.i() << " "
      << cl->center.j() << " "
      << cl->center.k() << "\n";
  
  // right daughter cell
  out << cr->id << " "
      << cr->parentId << " "
      << cr->treeId << " " 
      << cr->timeStep << " "
      << cr->cellCycle << " "
      //<< cr->area << " "
      //<< cr->longestWallLength << " "
      << cr->previousDivDir.i() << " "
      << cr->previousDivDir.j() << " "
      << cr->previousDivDir.k() << " "
      << cr->divDir.i() << " "
      << cr->divDir.j() << " "
      << cr->divDir.k() << " ";
      
  if( cr->cellCycle != 1 )
  {
    Point3d lastDir = cr->previousDivDir;
    Point3d curDir = cr->divDir;
    
    lastDir.normalize();
    curDir.normalize();
    
    out << acos( lastDir*curDir )*180./M_PI << " ";
  }
  else
    out << "NA ";
  
  // master cell file
  out << "0 "
      << cr->layerValue << " ";
  
  if( divType == DivisionType::ANTICLINAL )
    out << "0 ";
  else
    out << "1 ";
  
  // at last position of nucleus
  out << cr->center.i() << " "
      << cr->center.j() << " "
      << cr->center.k() << "\n";
  
  out.close();
}

// ---------------------------------------------------------------------

void exportModelProperties( const std::string &filename,
                            const std::size_t loopCounter,
                            const layerMap &firstLayerAppearances,
                            const std::pair<std::size_t, std::size_t> &divOccurrences,
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
    out << "'Sample' 'Anticlinal' 'Periclinal' ";
      
    for( std::size_t l = 0; l < totalSeqVec.size(); l++ )
      out << "'" << totalSeqVec.at(l) << "' ";

    out << "\n";
  }
  else
  {
    out.open( filename.c_str(), std::ofstream::out | std::ofstream::app );  
  }
  
  out << loopCounter << " ";
  out << divOccurrences.first << " ";
  out << divOccurrences.second << " ";
  
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
      out << "0 ";
  }
  out << "\n";
  out.close();
}

// ---------------------------------------------------------------------

void exportModelProperties( const std::string &filename,
                            const std::size_t loopCounter,
                            const std::map<std::string, std::size_t> &totalLayerCount,
                            const std::pair<std::size_t, std::size_t> &divOccurrences,
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
    out << "'Sample' 'Anticlinal' 'Periclinal' ";
      
    for( std::size_t l = 0; l < totalSeqVec.size(); l++ )
      out << "'" << totalSeqVec.at(l) << "' ";

    out << "\n";
  }
  else
  {
    out.open( filename.c_str(), std::ofstream::out | std::ofstream::app );  
  }
  
  out << loopCounter << " ";
  out << divOccurrences.first << " ";
  out << divOccurrences.second << " ";
  
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
      out << "0 ";
  }
  out << "\n";
  out.close();
}

// ---------------------------------------------------------------------

void exportDivisionProperties( const std::string &filename,
                               const cell& c,
                               const MyTissue::division_data& ddata,
                               const double angleThreshold,
                               const bool init )
{
  if( init )
  {
    std::ofstream out( filename.c_str(), std::ofstream::out );
    // header
    out << "'ID' 'ParentID' 'X' 'Y' 'Z' 'Time' 'CellCycle' 'Area' 'LongestCellWall' 'Lineage' 'Layer' 'DivisionAngle' 'AngleBetweenCurrentAndPreviousDivision' 'DivisionType'\n";
    out.close();
    return;
  }
  
  std::ofstream out( filename.c_str(), std::ofstream::out | std::ofstream::app );
  out << c->id << " ";
  
  if( c->cellCycle != 0 )
    out << c->parentId << " ";
  else
    out << "NA "; 
  
  out << c->center.i() << " "
      << c->center.j() << " "
      << c->center.k() << " "
      << c->timeStep << " "
      << c->cellCycle << " "
      << c->area << " "
      << c->longestWallLength << " "
      << c->treeId << " "
      << c->layerValue << " "
      << c->angle << " ";
  if( c->cellCycle != 0 )
  {
    double a, pa;
    a = c->angle;
    pa = c->previousAngle;
    if( a > 90. )
      a = 180. - a;
    
    if( pa > 90. )
      pa = 180. - pa;
    
    out << fabs( a - pa ) << " ";
  }
  else
    out << "NA ";
  
  DivisionType::type divType = ModelUtils::determineDivisionType( ddata, angleThreshold );
  if( divType == DivisionType::ANTICLINAL )
    out << "0\n";
  else
    out << "1\n";
      
  out.close();
}

// ---------------------------------------------------------------------

}