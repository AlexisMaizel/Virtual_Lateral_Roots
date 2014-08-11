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
  //parms(section.data(), "SurfTimeScale", _surfTimeScale);
  
  // depending on data set -> TODO
  _curSurface.readTriangulation( "/home/necrolyte/Data/VLR/Virtual_Lateral_Roots/FinalVLRForMatlab/triangulation-120830_raw.txt" );
  
  //this->growStep( 0.5 );
  //_curSurface.printTriangleProperties( 0 );
}

//----------------------------------------------------------------

void RealSurface::growStep( const double dt )
{ 
  _time += dt;
  _curSurface.interpolate( _time );
}

//----------------------------------------------------------------