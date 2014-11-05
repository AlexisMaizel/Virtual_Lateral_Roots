// Implement surface and surface point class
#include "surface.h"
#include "bezier.h"

#include <string>
#include <iostream>

#include <cstdio>
#include <cmath>
#include <util/parms.h>

using std::cout;
using std::cerr;
using std::endl;

// Clip to bounds
static double InBounds(double s, double min, double max) 
{
  if(s < min)
    return(min);
  else if(s > max)
    return(max);
  else
    return(s);
}

//----------------------------------------------------------------

Point3d SurfacePoint::Pos()
{
  return(pos);
}

//----------------------------------------------------------------

Point3d SurfacePoint::Pos() const
{
  return(pos);
}

//----------------------------------------------------------------

Point3d SurfacePoint::Normal()
{
  return(normal);
}

//----------------------------------------------------------------

Surface::Surface(util::Parms &parms, string section)
{
  string surffile;

  time = 0;

  // testing
  surface[0].LoadBezierSurface( "/home/necrolyte/Uni/LateralRootGrowth/TIFFStacks/120830Bezier.mgxv" );
  
  // Load surfaces 
  parms(section.data(), "Surfaces", surfaces);
  parms(section.data(), "SurfTimeScale", surfTimeScale);
  parms(section.data(), "SurfMaxDist", surfMaxDist);
  for(int i = 0; i < surfaces; i++) {
    // Bezier surface
    ostringstream key;
    key << "Surface" << i;
    parms(section.data(), key.str().data(), surffile);
    surface[i].Load(surffile);

    // Scaling constant
    key.str("");
    key << "SurfaceScale" << i;
    parms(section.data(), key.str().data(), surfScale[i]);

    // Time constant
    key.str("");
    key << "SurfaceTime" << i;
    parms(section.data(), key.str().data(), surfTime[i]);
  }
}

//----------------------------------------------------------------

// Zero a surface Point
void Surface::Zero(SurfacePoint &p)
{
  p.u = p.v = 0.0;
  CalcPos(p);
  CalcNormal(p);
}

//----------------------------------------------------------------

// Initial surface Point
void Surface::InitPoint(SurfacePoint &p, double u, double v)
{
  p.u = InBounds(u, 0.0, 1.0);
  p.v = InBounds(v, 0.0, 1.0);
  
  CalcPos(p);
  CalcNormal(p);
}

//----------------------------------------------------------------

// Determine parametric coords u and v such that the distance between
// the corresponding point for (u,v) is minimal to the desired point cp
bool Surface::SetPoint(SurfacePoint &p, SurfacePoint sp, Point3d cp)
{
  // Initial guess for p
  p.u = sp.u;
  p.v = sp.v;
  CalcPos(p);
  double lastd = -(DX * 1000);
  double count = 0;
  double du, dv;

  while(fabs(lastd - norm(cp - p.pos)) > DX * 2 && count++ < MAXSEARCHSTEPS)
  {
    // Save previous distance
    lastd = norm(cp - p.pos);

    // Calc partials
    SurfacePoint u1 = p, u2 = p;
    u1.u -= DX;
    u2.u += DX;
    CalcPos(u1);
    CalcPos(u2);
    du = (norm(cp - u1.pos) - norm(cp - u2.pos))/fabs(u1.u - u2.u);

    SurfacePoint v1 = p, v2 = p;
    v1.v -= DX;
    v2.v += DX;
    CalcPos(v1);
    CalcPos(v2);
    dv = (norm(cp - v1.pos) - norm(cp - v2.pos))/fabs(v1.v - v2.v);

    // Do line minimization
    double step = .01;
    while(step > DX)
    {
      SurfacePoint t = p;
      t.u += step * du;
      t.v += step * dv;
      CalcPos(t);
      if(norm(cp - t.pos) < norm(cp - p.pos))
      {
        step *= 1.11;
        p = t;
      }
      else
        step /= 2.0;
    }
  }

  CalcNormal(p);
  if(count >= MAXSEARCHSTEPS || norm(p.pos - cp) > surfMaxDist) {
    cerr << "Surface::SetPoint:Error Failed, point " << cp << 
            " closest point " << p.pos << " du " << du << " dv " << dv << 
            " count " << count << " distance " <<  norm(p.pos - cp) << endl;
  }

  return(true);
}

//----------------------------------------------------------------

// Advance time
void Surface::GrowStep(double dt)
{
  time += dt * surfTimeScale;
  int surf = 1;
  double surftime = 0;
  //while(time > surftime + surfTime[surf] && surf < surfaces - 1)
  //  surftime += surfTime[surf++ - 1];

  // setting for considering three bezier surfaces
  /*
  double halfTime = surfTime[surfaces-1]/2.;
  double maxTime = surfTime[surfaces-1];
	
  if( time < halfTime )
  {
    surf = 1;
    surfCurr.Interpolate( surface[surf-1], surface[surf], 
                          surfScale[surf-1], surfScale[surf],
                          (time - surftime)/halfTime );
  }
  else
  {
    surf = 2;
    surfCurr.Interpolate( surface[surf-1], surface[surf], 
                          surfScale[surf-1], surfScale[surf],
                          (time - surftime - halfTime)/halfTime );
  }
  */
  
  surfCurr.Interpolate( surface[surf-1], surface[surf], 
                        surfScale[surf-1], surfScale[surf],
                        (time - surftime)/surfTime[surf]);
}

//----------------------------------------------------------------

double Surface::GetTime() const
{
	return time;
}

//----------------------------------------------------------------

// Get position at current time 
void Surface::GetPos(SurfacePoint &p)
{
  CalcPos(p);
  CalcNormal(p);
}

//----------------------------------------------------------------

// Return distance between u and v
double Surface::Distance(SurfacePoint &u, SurfacePoint &v)
{
  // For now use Euclidean distance
  return(norm(u.pos - v.pos));
}

//----------------------------------------------------------------

// Calculate xyz postition
void Surface::CalcPos(SurfacePoint &p)
{
  p.u = InBounds(p.u, 0.0, 1.0);
  p.v = InBounds(p.v, 0.0, 1.0);
  p.pos = surfCurr.EvalCoord(p.u, p.v);
}

//----------------------------------------------------------------

// Calculate normal
void Surface::CalcNormal(SurfacePoint &p)
{
  // Calc normal from surface.
  SurfacePoint n, s, e, w;
  n.u = p.u + DX; n.v = p.v; 
  s.u = p.u - DX; s.v = p.v;
  e.u = p.u; e.v = p.v + DX;
  w.u = p.u; w.v = p.v - DX;
  CalcPos(n); CalcPos(s);
  CalcPos(e); CalcPos(w);
    
  p.normal = s.Pos() - n.Pos();
  p.normal = p.normal ^ (w.Pos() - e.Pos());
  p.normal.normalize();
}
