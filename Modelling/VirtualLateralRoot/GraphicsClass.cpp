#include "GraphicsClass.h"

/**
  @file   GraphicsClass.cpp
  @brief  Contains namespace for drawing methods used in the model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

namespace GraphicsClass
{

//----------------------------------------------------------------

void drawBezierCurve( const std::set<Point3d, lessXPos> &cps,
                      const util::Palette::Color &color,
                      GLUquadricObj *quadratic )
{
//   auto iter = cps.begin();
//   Point3d firstPos = *iter;
//   iter++;
//   for( ; iter != cps.end(); iter++ )
//   {
//     GraphicsClass::drawLine( firstPos, *iter, color, 3.5 );
//     firstPos = *iter;
//   }
  
  std::vector<Point3d> positions;
  for( auto iter = cps.begin(); iter != cps.end(); iter++ )
    positions.push_back( *iter );
  
  Point3d firstPos = ModelUtils::computeBezierPoint( positions, 0. );
  double deltaT = 0.02;
  for( double t = deltaT; t < 1.+deltaT; t += deltaT )
  {
    Point3d pos = ModelUtils::computeBezierPoint( positions, t );
    //GraphicsClass::drawLine( firstPos, pos, color, 3.5 );
    GraphicsClass::drawCylinder( firstPos, pos, color, quadratic, 1.5, 20, 20 );
    firstPos = pos;
  }
}

//----------------------------------------------------------------

void drawCylinder( const Point3d &pos1, const Point3d &pos2,
                   const util::Palette::Color &color,
                   GLUquadricObj *quadratic,
                   const double r, std::size_t slices,
                   std::size_t stacks )
{
  glColor4fv( color.c_data() );
  glShadeModel( GL_SMOOTH );
  GLfloat lmodel_ambient[] = { 0.01, 0.01, 0.01, 1.0 };
  glPushMatrix();
  glLightModelfv( GL_LIGHT_MODEL_AMBIENT, lmodel_ambient );
  Point3d mid = (pos2 + pos1)/2.;
  glTranslatef( mid.i(), mid.j(), mid.k() );
  glRotatef( 90., 1., 0., 0. );
  Point3d xaxisDir = Point3d( 1., 0., 0. );
  Point3d dir = pos2 - pos1;
  dir.normalize();
  double angle = 180./M_PI * acos( dir*xaxisDir );
  if( ModelUtils::determineSlope( dir ) < 0 )
    angle = 90. - angle;
  else
    angle += 90.;
  
  glRotatef( angle, 0., 1., 0. );
  gluCylinder( quadratic, r, r, 1.1*norm(pos2-pos1), slices, stacks );
  glPopMatrix();
}

//----------------------------------------------------------------

void drawLine( const Point3d &pos1, const Point3d &pos2,
               const util::Palette::Color &color,
               const double lineWidth )
{
  glLineWidth( lineWidth ); 
  glColor4fv( color.c_data() );
  glBegin( GL_LINES );
  glVertex3f( pos1.i(), pos1.j(), 2. );
  glVertex3f( pos2.i(), pos2.j(), 2. );
  glEnd();
}

//----------------------------------------------------------------

void drawControlPoint( const Point3d &pos,
                       const util::Palette::Color &color )
{
  double halfLength = 1.;
  glNormal3f( 0., 0., 1. );
  glPolygonMode( GL_FRONT, GL_FILL );
  glColor4fv( color.c_data() );
  glPushMatrix();
  glBegin( GL_QUADS );
  glVertex3f( pos.i()-halfLength, pos.j()+halfLength, 1. );
  glVertex3f( pos.i()+halfLength, pos.j()+halfLength, 1. );
  glVertex3f( pos.i()+halfLength, pos.j()-halfLength, 1. );
  glVertex3f( pos.i()-halfLength, pos.j()-halfLength, 1. );
  glEnd();
  glPopMatrix();
}

//----------------------------------------------------------------

void drawBezierSurface( const conpoi &cps,
                        const util::Palette::Color &color )
{
  glLineWidth(1.); 
  glColor4fv( color.c_data() );
  glShadeModel( GL_SMOOTH );
  double steps = 0.05;
  
  for( double u = 0.; u < 1.+steps; u+=steps )
  {
    glBegin( GL_LINE_STRIP );
    for( double v = 0.; v < 1.+steps; v+=steps )
    {
      Point3d pos = ModelUtils::computeBezierPoint( cps, u, v );
      glVertex3f( pos.i(), pos.j(), 0.9 );
    }
    glEnd();
  }
  
  for( double v = 0.; v < 1.+steps; v+=steps )
  {
    glBegin( GL_LINE_STRIP );
    for( double u = 0.; u < 1.+steps; u+=steps )
    {
      Point3d pos = ModelUtils::computeBezierPoint( cps, u, v );
      glVertex3f( pos.i(), pos.j(), 0.9 );
    }
    glEnd();
  }
}

//----------------------------------------------------------------

void drawSphere( const Point3d &pos, double r, std::size_t lats,
                 std::size_t longs, const util::Palette::Color &color,
                 GLUquadricObj *quadratic )
{
  glColor4fv( color.c_data() );
  glShadeModel( GL_SMOOTH );
  GLfloat lmodel_ambient[] = { 0.01, 0.01, 0.01, 1.0 };
  glPushMatrix();
  glLightModelfv( GL_LIGHT_MODEL_AMBIENT, lmodel_ambient );
  glTranslatef( pos.i(), pos.j(), 0. );
  gluSphere( quadratic, r, lats, longs );
  glPopMatrix();
}
 
//----------------------------------------------------------------
 
}