#ifndef ModelExporter_HH
#define ModelExporter_HH

#include "ModelHeader.h"
#include "ModelUtils.h"

/**
  @file   ModelExporter.h
  @brief  Namespace for exporting modelling results
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

namespace ModelExporter
{
  void exportLineageInformation( const std::string &filename,
                                 const cell& c,
                                 const MyTissue& T,
                                 const std::size_t timeStep,
                                 const bool init );
  
  void exportTimeAgainstCells( const std::string &filename,
                               const double dT,
                               const std::size_t numCells,
                               const bool init );
  
  void exportDivisionDaughterProperties( const std::string &filename,
                                         const cell& cl,
                                         const cell& cr,
                                         const MyTissue::division_data& ddata,
                                         const double angleThreshold,
                                         std::pair<std::size_t, std::size_t> &divOccurrences,
                                         const bool init );
  
  void exportDivisionProperties( const std::string &filename,
                                 const cell& c,
                                 const MyTissue::division_data& ddata,
                                 const double angleThreshold,
                                 const bool init );
  
  void exportCellWalls( const std::string &filename,
                        const cell& c,
                        const MyTissue &T,
                        const bool init );
}

#endif // ModelExporter_HH