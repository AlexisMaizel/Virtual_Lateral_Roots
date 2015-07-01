#ifndef GraphicsClass_HH
#define GraphicsClass_HH

#include "ModelHeader.h"
#include "ModelUtils.h"

/**
  @file   GraphicsClass.h
  @brief  Contains namespace for drawing methods used in the model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

namespace GraphicsClass
{
  void drawBezierCurve( const std::set<Point3d, lessXPos> &cps,
                        const util::Palette::Color &color,
                        GLUquadricObj *quadratic );
  
  void drawLine( const Point3d &pos1, const Point3d &pos2,
                 const util::Palette::Color &color,
                 const double lineWidth = 2.5 );
  
  void drawCylinder( const Point3d &pos1, const Point3d &pos2,
                     const util::Palette::Color &color,
                     GLUquadricObj *quadratic,
                     const double r, std::size_t slices,
                     std::size_t stacks );
  
  void drawSphere( const Point3d &pos, double r, std::size_t lats,
                   std::size_t longs, const util::Palette::Color &color,
                   GLUquadricObj *quadratic );
  
  void drawControlPoint( const Point3d &pos,
                         const util::Palette::Color &color );
  
  void drawBezierSurface( const conpoi &cps,
                          const util::Palette::Color &color );
}

#endif // GraphicsClass_HH