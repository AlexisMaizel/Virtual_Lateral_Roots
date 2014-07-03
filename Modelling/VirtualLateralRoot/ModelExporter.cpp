#include "ModelExporter.h"

/**
  @file   ModelExporter.cpp
  @brief  Namespace for exporting modelling results
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

namespace ModelExporter
{

// ---------------------------------------------------------------------

void initExportFile( const std::string &filename )
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
}

//----------------------------------------------------------------

void exportLineageInformation( const std::string &filename,
                               const cell& c,
                               const MyTissue& T,
                               const std::size_t timeStep )
{
  std::ofstream out( filename.c_str(), std::ofstream::out | std::ofstream::app );
  
  std::vector<Point3d> polygon;
  Point3d center;
  forall(const junction& j, T.S.neighbors(c))
  {
    polygon.push_back(j->sp.Pos());
    center += j->sp.Pos();
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

void initCellWallFile( const std::string &filename )
{
  std::ofstream out( filename.c_str(), std::ofstream::out );

  // header
  out << "ID Time CellCorners\n";
}

// ---------------------------------------------------------------------

void exportCellWalls( const std::string &filename,
                      const cell& c,
                      const MyTissue &T )
{
  std::ofstream out( filename.c_str(), std::ofstream::out | std::ofstream::app );
  
  out << c->id << " " << c->timeStep << " ";
  bool first = true;
  forall(const junction& j, T.S.neighbors(c))
  {
    if( first )
      first = false;
    else
      out << ",";
    
    out << j->sp.Pos().i() << "," 
        << j->sp.Pos().j() << "," 
        << j->sp.Pos().k();
  }
  
  out << "\n";
}

// ---------------------------------------------------------------------

void initDivisionDaughterFile( const std::string &filename )
{
  std::ofstream out( filename.c_str(), std::ofstream::out );

  // header
  out << "'ID' 'ParentID' 'Lineage' 'Time' 'CellCycle' 'XPos' 'YPos' 'ZPos' 'Area' 'LongestCellWall' "
      << "'Dir of last Division X' 'Dir of last Division Y' 'Dir of last Division Z' "
      << "'Dir of this Division X' 'Dir of this Division Y' 'Dir of this Division Z' "
      << "'Angle' 'Layer' 'DivisionType'\n";
}

// ---------------------------------------------------------------------

void exportDivisionDaughterProperties( const std::string &filename,
                                       const cell& cl,
                                       const cell& cr,
                                       const MyTissue::division_data& ddata,
                                       const double angleThreshold,
																			 std::pair<std::size_t, std::size_t> &divOccurrences )
{
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

void initDivisionFile( const std::string &filename )
{
  std::ofstream out( filename.c_str(), std::ofstream::out );

  // header
  out << "'ID' 'ParentID' 'X' 'Y' 'Z' 'Time' 'CellCycle' 'Area' 'LongestCellWall' 'Lineage' 'Layer' 'DivisionAngle' 'AngleBetweenCurrentAndPreviousDivision' 'DivisionType'\n";
}

// ---------------------------------------------------------------------

void exportDivisionProperties( const std::string &filename,
                               const cell& c,
                               const MyTissue::division_data& ddata,
                               const double angleThreshold )
{
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