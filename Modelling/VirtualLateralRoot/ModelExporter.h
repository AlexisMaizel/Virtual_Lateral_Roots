#ifndef ModelExporter_HH
#define ModelExporter_HH

#include "ModelHeader.h"

/**
  @file   ModelExporter.h
  @brief  Contains namespace for similarity measure tools
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

namespace ModelExporter
{
  void initExportFile( const std::string &filename );
  void exportLineageInformation( const std::string &filename,
                                 const cell& c,
                                 const MyTissue& T,
                                 const std::size_t timeStep );
}

#endif // ModelExporter_HH