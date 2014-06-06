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

void initDivisionFile( const std::string &filename )
{
  std::ofstream out( filename.c_str(), std::ofstream::out );

  // header
  out << "ID X Y Z Time CellCycle Area LongestCellWall Lineage Layer DivisionAngle DivisionType\n";
}

// ---------------------------------------------------------------------

void exportDivisionProperties( const std::string &filename,
                               const cell& c,
                               const MyTissue::division_data& ddata,
                               const double angleThreshold )
{
  std::ofstream out( filename.c_str(), std::ofstream::out | std::ofstream::app );
  
  out << c->id << " "
      << c->center.i() << " "
      << c->center.j() << " "
      << c->center.k() << " "
      << c->timeStep << " "
      << c->cellCycle << " "
      << c->area << " "
      << c->longestWallLength << " "
      << c->treeId << " "
      << c->layerValue << " "
      << ModelUtils::determineDivisionAngle( ddata ) << " ";
      DivisionType::type divType = ModelUtils::determineDivisionType( ddata, angleThreshold );
      if( divType == DivisionType::ANTICLINAL )
        out << "0\n";
      else
        out << "1\n";
      
  out.close();
}

// ---------------------------------------------------------------------

}