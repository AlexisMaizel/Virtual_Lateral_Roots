#include "RealSurface.h"

/**
  @file   RealSurface.cpp
  @brief  Class for handling triangulation of real data points
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

//----------------------------------------------------------------

RealSurface::RealSurface( util::Parms &parms, 
                          const std::string &section )
  : _time( 0. )
{
}

//----------------------------------------------------------------

void RealSurface::init( const double surfaceScale,
                        const std::string &fileName,
                        const bool useAutomaticContourPoints )
{
  _time = 0.;
  std::string dataset = fileName;
  if( dataset.compare( 0, 7, "Average") == 0 )
    dataset = "Average";
  
  std::string name = "triangulation-";
  name += dataset;
  if( useAutomaticContourPoints )
    name += "_auto.txt";
  else
    name += ".txt";
  _curSurface.readTriangulation( name, surfaceScale );
}

//----------------------------------------------------------------

void RealSurface::growStep( const double dt,
                            std::vector<SurfacePoint> &sps )
{ 
  _time += dt;
  _curSurface.interpolate( _time , sps );
}

//----------------------------------------------------------------

void RealSurface::initPos( SurfacePoint &sp )
{
  this->getPos( sp );
}

//----------------------------------------------------------------

void RealSurface::getPos( SurfacePoint &sp )
{
  this->calcPos( sp );
  this->calcNormal( sp );
}

//----------------------------------------------------------------

void RealSurface::setPos( SurfacePoint &sp, const Point3d &p )
{
  _curSurface.determinePosProperties( sp, p );
  this->calcNormal( sp );
}

//----------------------------------------------------------------

void RealSurface::calcPos( SurfacePoint &sp )
{
  sp.boundParameters();
  _curSurface.getCoord( sp );
}

//----------------------------------------------------------------

void RealSurface::calcNormal( SurfacePoint &sp )
{
  _curSurface.determineNormal( sp );
}

//----------------------------------------------------------------
