#ifndef ModelUtils_HH
#define ModelUtils_HH

#include "ModelHeader.h"

/**
  @file   ModelUtils.h
  @brief  Contains namespace for util methods used in the model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

namespace ModelUtils
{
  double determineLongestWallLength( const cell& c,
                                     const MyTissue& T );
  
  DivisionType::type determineDivisionType( const MyTissue::division_data& ddata,
                                            const double angleThreshold );
  
  double getDivisionAngle( const MyTissue::division_data& ddata );
  
  double determineDivisionAngle( const MyTissue::division_data& ddata );
}

#endif // ModelUtils_HH