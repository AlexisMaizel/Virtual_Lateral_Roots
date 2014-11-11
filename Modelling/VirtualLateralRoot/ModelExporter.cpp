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
                               const std::size_t timeStep,
                               const bool init,
                               const std::size_t surfaceType )
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
  
  std::vector<Point2d> polygon;
  Point3d center = Point3d( 0., 0., 0. );
  forall(const junction& j, T.S.neighbors(c))
  {
    Point3d pos;
    
    if( surfaceType == 0 )
      pos = j->sp.Pos();
    else
      pos = Point3d( j->tp.Pos().i(), j->tp.Pos().j(), 0. );
    
    polygon.push_back( Point2d( pos.i(), pos.j() ) );
    center += pos;
  }
  
  center /= polygon.size();
  c->center = center;
  c->area = geometry::polygonArea(polygon);
  c->timeStep = timeStep;
  
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
  std::ofstream out( filename.c_str(), std::ofstream::out | std::ofstream::app );
  
  // header
  if( init )
    out << "dT Cells\n";
  else
    out << dT << " " << numCells << "\n";
  
  out.close();
}

// ---------------------------------------------------------------------

void exportCellWalls( const std::string &filename,
                      const cell& c,
                      const MyTissue &T,
                      const bool init,
                      const std::size_t surfaceType )
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
  
    if( surfaceType == 0 )
    {
      out << j->sp.Pos().i() << "," 
          << j->sp.Pos().j() << "," 
          << j->sp.Pos().k();
    }
    else
    {
      out << j->tp.Pos().i() << "," 
          << j->tp.Pos().j() << "," 
          << 0.;
    }
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
    out << "'ID' 'ParentID' 'Lineage' 'Time' 'CellCycle' 'XPos' 'YPos' 'ZPos' 'Area' 'LongestCellWall' "
        << "'Dir of last Division X' 'Dir of last Division Y' 'Dir of last Division Z' "
        << "'Dir of this Division X' 'Dir of this Division Y' 'Dir of this Division Z' "
        << "'Angle' 'Layer' 'DivisionType'\n";
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
      << cl->center.i() << " "
      << cl->center.j() << " "
      << cl->center.k() << " "
      << cl->area << " "
      << cl->longestWallLength << " "
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
  
  out << cl->layerValue << " ";
  
  DivisionType::type divType = ModelUtils::determineDivisionType( ddata, angleThreshold );
  if( divType == DivisionType::ANTICLINAL )
  {
    divOccurrences.first++;
    out << "0\n";
  }
  else
  {
    divOccurrences.second++;
    out << "1\n";
  }
      
  
  // right daughter cell
  out << cr->id << " "
      << cr->parentId << " "
      << cr->treeId << " " 
      << cr->timeStep << " "
      << cr->cellCycle << " "
      << cr->center.i() << " "
      << cr->center.j() << " "
      << cr->center.k() << " "
      << cr->area << " "
      << cr->longestWallLength << " "
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
  
  out << cr->layerValue << " ";
  
  if( divType == DivisionType::ANTICLINAL )
    out << "0\n";
  else
    out << "1\n";
  
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