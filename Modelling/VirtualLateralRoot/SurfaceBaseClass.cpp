#include "SurfaceBaseClass.h"

/**
  @file   SurfaceBaseClass.cpp
  @brief  Base class for handling surface properties
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

//----------------------------------------------------------------

void SurfaceBaseClass::initSurfaceProperties()
{
}

//----------------------------------------------------------------

void SurfaceBaseClass::initPos( SurfacePoint &sp )
{
}

//----------------------------------------------------------------

void SurfaceBaseClass::calcPos( SurfacePoint &sp )
{
}

//----------------------------------------------------------------

void SurfaceBaseClass::calcNormal( SurfacePoint &sp )
{
}

//----------------------------------------------------------------

void SurfaceBaseClass::setPos( SurfacePoint &sp, const Point3d &cp )
{
}

//----------------------------------------------------------------

void SurfaceBaseClass::getPos( SurfacePoint &sp )
{
}

//----------------------------------------------------------------

void SurfaceBaseClass::growStep( const double dt )
{
}

//----------------------------------------------------------------

void SurfaceBaseClass::growStep( const double dt, std::vector<SurfacePoint> &sps )
{
}

//----------------------------------------------------------------