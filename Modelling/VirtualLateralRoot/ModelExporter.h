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
  void initExportFile( const std::string &filename );
  
  void exportLineageInformation( const std::string &filename,
                                 const cell& c,
                                 const MyTissue& T,
                                 const std::size_t timeStep );
  
  void initDivisionDaughterFile( const std::string &filename );
  
  void initDivisionFile( const std::string &filename );
  
  void exportDivisionDaughterProperties( const std::string &filename,
                                         const cell& cl,
                                         const cell& cr,
                                         const MyTissue::division_data& ddata,
                                         const double angleThreshold );
  
  void exportDivisionProperties( const std::string &filename,
                                 const cell& c,
                                 const MyTissue::division_data& ddata,
                                 const double angleThreshold );
}

#endif // ModelExporter_HH