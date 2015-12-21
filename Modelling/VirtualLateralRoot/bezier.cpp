// bezier.cpp - Load Bezier from files output by interactive editor.
//
// Handles multiple patch files (advanced editor) and combines into 
// one surface parameterized by u, v, both between 0 and 1.

#include "bezier.h"

#include <cctype>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#define USEONESURFACE

//----------------------------------------------------------------

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

// Return points (debugging)
float *Bezier::Points()
{
  return(ptsV);
}

//----------------------------------------------------------------
#ifndef USEONESURFACE
// Load from file written by patch editor
void Bezier::Load( const std::string bezFile )
{
  ifstream bIn(bezFile.c_str());
  if(!bIn) {
    cerr << "Bezier::Bezier:Error opening " << bezFile << endl;
    exit(1);
  }
  
  string buff("");

  // Find first patch
  while(bIn) {
    bIn >> buff;
    transform (buff.begin(), buff.end(), buff.begin(), ::tolower);
    if(buff.substr(0, 5) == "patch")
      break;
  }
  unsigned int patchcount = 0;
  //int pi = 0;
  while(bIn && patchcount < MAXPATCHES * MAXPATCHES) {
    if(!bIn || buff.substr(0, 5) != "patch")
      break;
    // skip color info  name
    bIn >>buff>>buff>>buff>>buff>>buff>>buff>>buff>>buff>>buff>>buff; 
    // Skip nhbd info
    bIn >>buff>>buff>>buff>>buff>>buff>>buff; 
    bIn >>buff>>buff>>buff>>buff; 
    bIn >>buff>>buff>>buff>>buff>>buff>>buff; 

    // Load points
    float *pts = ptsV + PATCHMAP[patchcount] * PATCHSIZE;
    for(unsigned int i = 0; i < PATCHPOINTS; i++)
      for(unsigned int j = 0; j < PATCHPOINTS; j++)
        for(int k = 0; k < 3; k++)
          bIn >> *pts++;
    patchcount++;
    bIn >> buff; // next patch name
    transform (buff.begin(), buff.end(), buff.begin(), ::tolower);
  }
  count = (unsigned int)pow((double)patchcount, 0.5);
  if(count*count != patchcount)
    cerr << "Bezier::Load:Error " << patchcount << 
                                 " patches loaded s/b power of 2" << endl; 
  else
    cerr << "Bezier::Load:Info " << patchcount << " patches loaded" << endl;
}

//----------------------------------------------------------------

void Bezier::LoadGrowthBezier( const std::string bezFile,
                               const std::size_t numPatches,
                               const std::size_t numControlPoints )
{
  ifstream bIn(bezFile.c_str());
  if(!bIn) {
    cerr << "Bezier::Bezier:Error opening " << bezFile << endl;
    exit(1);
  }
  
  for( std::size_t i = 0; i < MAXPATCHES * MAXPATCHES * PATCHSIZE; i++ )
    ptsV[i] = 0.;
  
  unsigned int patchcount = 0;
  while(bIn && patchcount < MAXPATCHES * MAXPATCHES)
  {
    if( !bIn || patchcount == numPatches )
      break;
    
    float *pts = ptsV + PATCHMAP[patchcount] * PATCHSIZE;
    for(unsigned int i = 0; i < PATCHPOINTS; i++)
      for(unsigned int j = 0; j < PATCHPOINTS; j++)
        for(int k = 0; k < 3; k++)
          bIn >> *pts++;
    
    patchcount++;
  }
  
  count = (unsigned int)pow((double)(patchcount), 0.5);
  if(count*count != patchcount)
    cerr << "Bezier::Load:Error " << patchcount << 
                                 " patches loaded s/b power of 2" << endl;
}

//----------------------------------------------------------------

// Copy and scale
void Bezier::Scale( Bezier &b, double scale )
{
  float *pts = ptsV;
  float *bpts = b.ptsV;
  count = b.count;
  unsigned int size = count * count * PATCHSIZE;
  for(unsigned int i = 0; i < size; i++)
    *pts++ = *bpts++ * scale;
}

//----------------------------------------------------------------

// Interpolate (plus scale) between two patches
void Bezier::Interpolate( Bezier &src1, Bezier &src2, 
                          double scale1, double scale2, double s )
{
  if(src1.count != src2.count)
    cerr << "Bezier::Interpolate Error must have matching patch count" << endl;
  count = src1.count;
  
  s = InBounds(s, 0.0, 1.0);
  for(unsigned int i = 0; i < src1.count; i++)
    for(unsigned int j = 0; j < src1.count; j++)
    {
      int offset = ((i * MAXPATCHES) + j) * PATCHSIZE;
      float *pts = ptsV + offset;
      float *s1pts = src1.ptsV + offset;
      float *s2pts = src2.ptsV + offset;
      for(unsigned int k = 0; k < PATCHSIZE; k++)
        *pts++ = ((1.0 - s) * scale1 * *s1pts++ + s * scale2 * *s2pts++);
    }
}

//----------------------------------------------------------------

// Print out points (debugging)
void Bezier::Print()
{
  cout << "Patch count " << count << endl;
  /*
  for(unsigned int k = 0; k < count * count; k++)
  {
    for(int i = 0; i < 16; i++)
    {
      for(int j = 0; j < 3; j++)
        cout << *(ptsV + k*PATCHSIZE + i * 3 + j) << " ";
      cout << endl;
    }
    cout << endl;
  }*/
  
  
  for( std::size_t i = 0; i < MAXPATCHES * MAXPATCHES * PATCHSIZE; i++ )
  {
    std::cout << ptsV[i] << " ";
    if( (i+1)%(PATCHPOINTS*3) == 0 )
      std::cout << std::endl;
  }
}

//----------------------------------------------------------------

// Given u,v in parametric space, return x,y,z
Point3d Bezier::EvalCoord( double u, double v )
{
  // Make sure u and v between 0.0 and 1.0
  u = InBounds(u, 0.0, 1.0);
  v = InBounds(v, 0.0, 1.0);

  unsigned int iu = 0;
  unsigned int iv = 0;

  // Scale so entire surface fits between 0 and 1
  if(count > 1)
  {
    u *= count;
    v *= count;

    iu = (unsigned int)u;
    iv = (unsigned int)v;

    if(iu == count)
    {
      u = 1.0;
      iu = count - 1;
    }
    else
      u = u - iu;
    
    if(iv == count)
    {
      v = 1.0;
      iv = count - 1;
    }
    else
      v = v - iv;
  }

  // Find correct patch
  float *pts = ptsV + (iv * MAXPATCHES + iu) * PATCHSIZE;

  //std::cout << "u: " << u << " v: " << v << " iu: " << iu << " iv: " << iv << std::endl;
  
  // Evaluate
  double x = 0, y = 0, z = 0;
  for(unsigned int j = 0; j < PATCHPOINTS; j++)
    for(unsigned int i = 0; i < PATCHPOINTS; i++)
    {
      double s = (double)Choose(PATCHPOINTS - 1, i) * std::pow(u, (int)i) * 
                          pow(1.0 - u, (int)(PATCHPOINTS - 1 - i)) *
                 (double)Choose(PATCHPOINTS - 1, j) * std::pow(v, (int)j) * 
                    pow(1.0 - v, (int)(PATCHPOINTS - 1 - j));
      x += s * *pts++;
      y += s * *pts++;
      z += s * *pts++;
    }
  Point3d p(x, y, z);
  //std::cout << p << std::endl;
  return(p);
}

//----------------------------------------------------------------

// Find binomial coefficient
int Bezier::Choose( unsigned int n, unsigned int k )
{
  static bool first = true;
  static int choose[PATCHPOINTS][PATCHPOINTS];

  // Only calculate what we need
  if(n > PATCHPOINTS || k > n) {
    cerr << "Bezier::Choose:Error Can't compute " << n << " choose " << k 
                                                               << endl;
    return(0);
  }

  // If first time in, fill the table
  if(first) {
    first = false;
    for(unsigned int j = 0; j < PATCHPOINTS; j++)
      for(unsigned int i = j; i < PATCHPOINTS; i++)
        if(j == 0 || j == i)
          choose[i][j] = 1;
        else 
          choose[i][j] = choose[i-1][j] + choose[i-1][j-1];
  }

  // Just grab value from the table
  return(choose[n][k]);
}

#else
//----------------------------------------------------------------

void Bezier::Load( const std::string bezFile )
{
  std::ifstream bIn( bezFile.c_str() );
  if(!bIn)
  {
    cerr << "Bezier::Bezier:Error opening " << bezFile << endl;
    exit(1);
  }

  std::size_t numPatches = 0;
  std::string line;
  // determine the number of patches
  // choices are for example 1, 4 or 16
  while(bIn)
  {
    getline( bIn, line );
    std::transform( line.begin(), line.end(), line.begin(), ::tolower );
    if(line.substr(0, 5) == "patch")
      numPatches++;
  }
  
  bIn.close();
  bIn.open( bezFile.c_str() );
  if(!bIn)
  {
    cerr << "Bezier::Bezier:Error opening " << bezFile << endl;
    exit(1);
  }
  
  const std::size_t elements = 3*sqrt(numPatches) + 1;
  
  // initialize control points matrix
  _cpMatrix.resize( elements );
  for( std::size_t r=0;r<elements;r++ )
    _cpMatrix.at(r).resize( elements );
  
  // begin at first patch
  for( std::size_t l=0; l<8; l++ )
    getline( bIn, line );
  std::transform( line.begin(), line.end(), line.begin(), ::tolower );
  
  unsigned int patchcount = 0;
  while( bIn && patchcount < numPatches )
  {
    std::transform( line.begin(), line.end(), line.begin(), ::tolower );
    if(!bIn || line.substr(0, 5) != "patch")
      break;
    
    // skip color info
    getline( bIn, line );
    // Skip nhbd info
    getline( bIn, line );
    getline( bIn, line );
    getline( bIn, line );

    // read first line of patch
    getline( bIn, line );
    std::stringstream ss( line );
    
    // read points: TODO: at the moment this is applied for having only
    // one or 4 bezier patches
    for(std::size_t i = 0; i < PATCHPOINTS; i++)
    { 
      for(std::size_t j = 0; j < PATCHPOINTS; j++)
      {
        std::size_t xIndex, yIndex;
        
        switch( patchcount )
        {
          case 0: default: xIndex = j; yIndex = i; break;
          case 1: xIndex = PATCHPOINTS - 1 + j; yIndex = i; break;
          case 2: xIndex = j; yIndex = PATCHPOINTS - 1 + i; break;
          case 3: xIndex = PATCHPOINTS - 1 + j; yIndex = PATCHPOINTS - 1 + i; break;
        }
        
        double x,y,z;
        ss >> x >> y >> z;
        Point3d pos( x, y, z );
        _cpMatrix.at(yIndex).at(xIndex) = pos;
      }
      
      getline( bIn, line );
      ss.str( line );
    }
    
    patchcount++;
  }

  count = (unsigned int)pow((double)patchcount, 0.5); 
  //this->Print();
  
  bIn.close();
}

//----------------------------------------------------------------

void Bezier::LoadGrowthBezier( const std::string bezFile,
                               const std::size_t numPatches,
                               const std::size_t numControlPoints )
{
  std::ifstream bIn( bezFile.c_str() );
  if(!bIn)
  {
    cerr << "Bezier::Bezier:Error opening " << bezFile << endl;
    exit(1);
  }

  std::string line;
  
  // initialize control points matrix
  _cpMatrix.resize( numControlPoints );
  for( std::size_t r=0;r<numControlPoints;r++ )
    _cpMatrix.at(r).resize( numControlPoints );
  
  unsigned int patchcount = 0;
  while( bIn && patchcount < numPatches )
  {
    // read first line of patch
    getline( bIn, line );
    std::stringstream ss( line );
    
    // read points
    for(std::size_t i = 0; i < numControlPoints; i++)
    { 
      for(std::size_t j = 0; j < numControlPoints; j++)
      { 
        double x,y,z;
        ss >> x >> y >> z;
        Point3d pos( x, y, z );
        _cpMatrix.at(i).at(j) = pos;
      }
      
      getline( bIn, line );
      ss.str( line );
    }
    
    patchcount++;
  }

  count = (unsigned int)pow((double)patchcount, 0.5); 
  //this->Print();
  
  bIn.close();
}

//----------------------------------------------------------------

// Copy and scale
void Bezier::Scale( Bezier &b, double scale )
{
  conpoi mat = b._cpMatrix;
  count = b.count;
  for(std::size_t i = 0; i < _cpMatrix.size(); i++) 
    for(std::size_t j = 0; j < _cpMatrix.at(i).size(); j++)
      _cpMatrix.at(i).at(j) = mat.at(i).at(j) * scale;
}

//----------------------------------------------------------------

// Interpolate (plus scale) between two patches
void Bezier::Interpolate( Bezier &src1, Bezier &src2, 
                          double scale1, double scale2, double s )
{
  if(src1.count != src2.count)
    cerr << "Bezier::Interpolate Error must have matching patch count" << endl;
  count = src1.count;
  
  s = InBounds(s, 0.0, 1.0);
  conpoi mat1 = src1._cpMatrix;
  conpoi mat2 = src2._cpMatrix;
  
  // initialize control points matrix
  _cpMatrix.resize( mat1.size() );
  for( std::size_t r=0;r<_cpMatrix.size();r++ )
    _cpMatrix.at(r).resize( mat1.size() );
  
  for(std::size_t i = 0; i < _cpMatrix.size(); i++) 
    for(std::size_t j = 0; j < _cpMatrix.at(i).size(); j++)
      _cpMatrix.at(i).at(j) = ((1. - s) * scale1 * mat1.at(i).at(j) +
                              s * scale2 * mat2.at(i).at(j));
}

//----------------------------------------------------------------

// Print out points (debugging)
void Bezier::Print()
{ 
  cout << "Patch count " << count << endl;
  
  for(std::size_t i = 0; i < _cpMatrix.size(); i++)
  {
    for(std::size_t j = 0; j < _cpMatrix.at(i).size(); j++)
      std::cout << _cpMatrix.at(_cpMatrix.size()-1-i).at(j) << " ";
    
    std::cout << std::endl;
  }
}

//----------------------------------------------------------------

void Bezier::applyGrowthOnlyInHeight( const Bezier &surfaceE )
{
  conpoi sur = surfaceE.getControlPoints();
  for(std::size_t i = 0; i < _cpMatrix.size(); i++)
  {
    for(std::size_t j = 0; j < _cpMatrix.at(i).size(); j++)
    {
      // copy all x positions from the last time step
      _cpMatrix.at(i).at(j).i() = sur.at(i).at(j).i();
      _cpMatrix.at(i).at(j).j() = sur.at(0).at(j).j() + (double)i*5.;
    }
  }
}

//----------------------------------------------------------------

void Bezier::increaseDomeTipHeight( const double topPoint )
{
  double sideTopPoints = topPoint + 10.; // 180.; -> default
  _cpMatrix.at(6).at(3).j() = topPoint;
  _cpMatrix.at(6).at(2).j() = sideTopPoints;
  _cpMatrix.at(6).at(4).j() = sideTopPoints;
  for(std::size_t i = 1; i < _cpMatrix.size()-1; i++)
  {
    for(std::size_t j = 2; j < _cpMatrix.at(i).size()-2; j++)
    {
      // update y positions
      double t = (double)i/(double)(_cpMatrix.at(i).size()-1);
      _cpMatrix.at(i).at(j).j() =
      (1-t) * _cpMatrix.at(0).at(j).j() + t * _cpMatrix.at(6).at(j).j();
    }
  }
}

//----------------------------------------------------------------

void Bezier::setRadialSurface( const int numControlPoints,
                               const bool start )
{
  // initialize control points matrix
  _cpMatrix.resize( numControlPoints );
  for( std::size_t r=0;r<numControlPoints;r++ )
    _cpMatrix.at(r).resize( numControlPoints );
  
  // first dimension is the row while second one is the column
  
  double r = 300.;
  double rBottom = 7.*r/10.;
  double rTop = 8.*r;
  double angle = 15.;
  double angleStep = 25.;
  int variant = 1;
  
  // set bottom boundary points
  for(std::size_t j = 0; j < _cpMatrix.front().size(); j++)
  {
    _cpMatrix.front().at(j).i() = -rBottom * cos( angle * M_PI/180. );
    _cpMatrix.front().at(j).j() = rBottom * sin( angle * M_PI/180. );
    _cpMatrix.front().at(j).k() = 0.;
    angle += angleStep;
  }
  
  angle = 15.;
  
  // set top boundary points
  if( start )
  {
    for(std::size_t j = 0; j < _cpMatrix.back().size(); j++)
    {
      _cpMatrix.back().at(j).i() = -r * cos( angle * M_PI/180. );
      _cpMatrix.back().at(j).j() = r * sin( angle * M_PI/180. );
      _cpMatrix.back().at(j).k() = 0.;
      angle += angleStep;
    }
  }
  else
  {
    if( variant == 0 )
    {
      for(std::size_t j = 0; j < _cpMatrix.back().size(); j++)
      {
        _cpMatrix.back().at(j).i() = 1.1 * -r * cos( angle * M_PI/180. );
        _cpMatrix.back().at(j).j() = 4. * r * sin( angle * M_PI/180. );
        
        if( j == 1 || j == 5 )
          _cpMatrix.back().at(j).j() *= 1.5;
        
        if( j == 3 )
          _cpMatrix.back().at(j).j() /= 1.1;
          
        _cpMatrix.back().at(j).k() = 0.;
        angle += angleStep;
      }
    }
    else
    {
      // manually set all interior and top points of the surface in order
      // to get a more parabolic growth behavior
      
      for(std::size_t i = 1; i < _cpMatrix.size(); i++)
      {
        double s = (double)i/(double)(_cpMatrix.size()-1);
        
        // mid is just growing in height and interpolated linearly
        Point3d top = Point3d( -r * cos( 90. * M_PI/180. ),
                               rTop/1.5 * sin( 90. * M_PI/180. ), 0 );
        Point3d bottom = Point3d( -rBottom * cos( 90. * M_PI/180. ),
                               rBottom * sin( 90. * M_PI/180. ), 0 );
        _cpMatrix.at(i).at(3) = (1.-s) * bottom + s * top;
        
        // set the x positions for the inner columns
        double t = -r * cos( 40. * M_PI/180. )*1.5;
        double b = -rBottom * cos( 40. * M_PI/180. );
        _cpMatrix.at(i).at(1).i() = (1.-s) * b + s * t;
        t = -r * cos( 65. * M_PI/180. )*3.5;
        b = -rBottom * cos( 65. * M_PI/180. );
        _cpMatrix.at(i).at(2).i() = (1.-s) * b + s * t;
        t = -r * cos( 115. * M_PI/180. )*3.5;
        b = -rBottom * cos( 115. * M_PI/180. );
        _cpMatrix.at(i).at(4).i() = (1.-s) * b + s * t;
        t = -r * cos( 140. * M_PI/180. )*1.5;
        b = -rBottom * cos( 140. * M_PI/180. );
        _cpMatrix.at(i).at(5).i() = (1.-s) * b + s * t;
        
        // most left
        top = Point3d( -r * cos( 15. * M_PI/180. ),
                       r * sin( 15. * M_PI/180. ), 0 );
        bottom = Point3d( -rBottom * cos( 15. * M_PI/180. ),
                          rBottom * sin( 15. * M_PI/180. ), 0 );
        _cpMatrix.at(i).front() = (1.-s) * bottom + s * top;
        
        // most right
        top = Point3d( -r * cos( 165. * M_PI/180. ),
                       r * sin( 165. * M_PI/180. ), 0 );
        bottom = Point3d( -rBottom * cos( 165. * M_PI/180. ),
                          rBottom * sin( 165. * M_PI/180. ), 0 );
        _cpMatrix.at(i).back() = (1.-s) * bottom + s * top;
        
      }
      
      std::size_t j = 1;
      double height = fabs( rTop - rBottom );
      double para[] = { rBottom+0.1*height,
                        rBottom+0.6*height,
                        rBottom+0.8*height,
                        rBottom+0.9*height,
                        rBottom+0.95*height,
                        rTop };
      double prevY = _cpMatrix.at(0).at(j).j();
      for( std::size_t c = 1; c<=6; c++ )
        _cpMatrix.at(c).at(j).j() = prevY + para[c-1]/8.;
      
      j = 2;
      prevY = _cpMatrix.at(0).at(j).j();
      for( std::size_t c = 1; c<=6; c++ )
        _cpMatrix.at(c).at(j).j() = prevY + para[c-1]/2.;
      
      j = 4;
      prevY = _cpMatrix.at(0).at(j).j();
      for( std::size_t c = 1; c<=6; c++ )
        _cpMatrix.at(c).at(j).j() = prevY + para[c-1]/2.;
      
      j = 5;
      prevY = _cpMatrix.at(0).at(j).j();
      for( std::size_t c = 1; c<=6; c++ )
        _cpMatrix.at(c).at(j).j() = prevY + para[c-1]/8.;
    }
  }
  
  if( variant == 0 || start )
  {
    // interpolate in between
    for(std::size_t i = 1; i < _cpMatrix.size()-1; i++)
      for(std::size_t j = 0; j < _cpMatrix.at(i).size(); j++)
      {
        Point3d top = _cpMatrix.back().at(j);
        Point3d bottom = _cpMatrix.front().at(j);
        double dist = r;
        double s = (double)i/(double)(_cpMatrix.size()-1);
        _cpMatrix.at(i).at(j) = (1.-s) * bottom + s * top;
      }
  }
}

//----------------------------------------------------------------

void Bezier::setRadialSurfaces( const std::size_t surf )
{
  std::size_t numControlPoints = 7;
  // initialize control points matrix
  _cpMatrix.resize( numControlPoints );
  for( std::size_t r=0;r<numControlPoints;r++ )
    _cpMatrix.at(r).resize( numControlPoints );
  
  double r = 300.;
  double rBottom = 6.*r/10.;
  double rTop = 6.*r;
  double rMid = 3.*r;
  double angle = 15.;
  double angleStep = 25.;
  
  // start
  if( surf == 0 )
  {
    // set bottom and top boundary points
    for(std::size_t j = 0; j < numControlPoints; j++)
    {
      double a = angle * M_PI/180.;
      // bottom
      if( j >= 2 && j <= 4 )
        this->setControlPoint( 0, j,
                               -rBottom*cos( a ),
                               rBottom*sin( a ) - 30. );
      else
        this->setControlPoint( 0, j,
                               -rBottom*cos( a ),
                               rBottom*sin( a ) );
      // top
      if( j >= 2 && j <= 4 )
        this->setControlPoint( 6, j,
                               -r*cos( a ),
                               r*sin( a ) - 30. );
      else
        this->setControlPoint( 6, j,
                               -r*cos( a ),
                               r*sin( a ) );
      angle += angleStep;
    }
    
    // interpolate in between
    for(std::size_t i = 1; i < _cpMatrix.size()-1; i++)
      for(std::size_t j = 0; j < _cpMatrix.at(i).size(); j++)
      {
        Point3d top = _cpMatrix.back().at(j);
        Point3d bottom = _cpMatrix.front().at(j);
        double s = (double)i/(double)(_cpMatrix.size()-1);
        _cpMatrix.at(i).at(j) = (1.-s) * bottom + s * top;
      }
  }
  // mid and end
  else if( surf == 1 || surf == 2 )
  {
    // set bottom boundary points
    for(std::size_t j = 0; j < numControlPoints; j++)
    {
      double a = angle * M_PI/180.;
      // bottom
      if( j >= 2 && j <= 4 )
        this->setControlPoint( 0, j,
                               -rBottom*cos( a ),
                               rBottom*sin( a ) - 30. );
      else if( j == 0 || j == 6 )
        this->setControlPoint( 0, j,
                               -rBottom*cos( a ),
                               rBottom*sin( a ) + 20. );
      else
        this->setControlPoint( 0, j,
                               -rBottom*cos( a ),
                               rBottom*sin( a ) );
      angle += angleStep;
    }
    
    for(std::size_t i = 1; i < _cpMatrix.size(); i++)
    {
      double s = (double)i/(double)(_cpMatrix.size());
      
      // mid is just growing in height and interpolated linearly
      double a = 90.* M_PI/180.;
      Point3d top = Point3d( -rMid * cos( a ), rMid * sin( a ), 0 );
      Point3d bottom = Point3d( -rBottom * cos( a ), rBottom * sin( a ), 0 );
      _cpMatrix.at(i).at(3) = (1.-s) * bottom + s * top;
      
      // set the x positions for the inner columns
      a = 15.* M_PI/180.;
      double t = -r * cos( a ) * 1.7;
      double b = -rBottom * cos( a );
      _cpMatrix.at(i).at(0).i() = (1.-s) * b + s * t;
      
      a = 40.* M_PI/180.;
      t = -r * cos( a ) * 1.8;
      b = -rBottom * cos( a );
      _cpMatrix.at(i).at(1).i() = (1.-s) * b + s * t;
      
      a = 65.* M_PI/180.;
      t = -r * cos( a ) * 4.2;
      b = -rBottom * cos( a );
      _cpMatrix.at(i).at(2).i() = (1.-s) * b + s * t;
      
      a = 115.* M_PI/180.;
      t = -r * cos( a ) * 4.2;
      b = -rBottom * cos( a );
      _cpMatrix.at(i).at(4).i() = (1.-s) * b + s * t;
      
      a = 140.* M_PI/180.;
      t = -r * cos( a ) * 1.8;
      b = -rBottom * cos( a );
      _cpMatrix.at(i).at(5).i() = (1.-s) * b + s * t;

      a = 165.* M_PI/180.;
      t = -r * cos( a ) * 1.7;
      b = -rBottom * cos( a );
      _cpMatrix.at(i).at(6).i() = (1.-s) * b + s * t;
    }
    
    std::size_t j = 0;
    double height = fabs( rMid - r );
    double para[] = { r+0.1*height,
                      r+0.6*height,
                      r+0.8*height,
                      r+0.9*height,
                      r+0.95*height,
                      rMid };
    double para2[] = { 0.01*height,
                      0.1*height,
                      0.2*height,
                      0.4*height,
                      0.8*height,
                      2.*height };
    double prevY = _cpMatrix.at(0).at(j).j();
    for( std::size_t c = 1; c<=6; c++ )
      _cpMatrix.at(c).at(j).j() = prevY + para2[c-1]/16.;
    
    j = 1;
    prevY = _cpMatrix.at(0).at(j).j();
    for( std::size_t c = 1; c<=6; c++ )
      _cpMatrix.at(c).at(j).j() = prevY + para[c-1]/8.;
    
    j = 2;
    prevY = _cpMatrix.at(0).at(j).j();
    for( std::size_t c = 1; c<=6; c++ )
      _cpMatrix.at(c).at(j).j() = prevY + para[c-1]/2.;
    
    j = 4;
    prevY = _cpMatrix.at(0).at(j).j();
    for( std::size_t c = 1; c<=6; c++ )
      _cpMatrix.at(c).at(j).j() = prevY + para[c-1]/2.;
    
    j = 5;
    prevY = _cpMatrix.at(0).at(j).j();
    for( std::size_t c = 1; c<=6; c++ )
      _cpMatrix.at(c).at(j).j() = prevY + para[c-1]/8.;
    
    j = 6;
    prevY = _cpMatrix.at(0).at(j).j();
    for( std::size_t c = 1; c<=6; c++ )
      _cpMatrix.at(c).at(j).j() = prevY + para2[c-1]/16.;
   
    if( surf == 2 )
    {
      double height = fabs( rTop - r );
      double para[] = { r+0.1*height,
                      r+0.1*height,
                      r+0.2*height,
                      r+0.6*height,
                      r+0.2*height,
                      r+0.1*height,
                      r+0.1*height };
      // set top boundary points
      for( std::size_t i = 1; i < numControlPoints; i++ )
        for( std::size_t j = 1; j < numControlPoints-1; j++ )
        { 
          if( j == 1 || j == 5 )
            _cpMatrix.at(i).at(j).j() += para[j] * i/16.;
          else
            _cpMatrix.at(i).at(j).j() += para[j];
        }
        
      _cpMatrix.at(6).at(0).j() += 60.;
      _cpMatrix.at(5).at(0).j() += 30.;
      _cpMatrix.at(6).at(6).j() += 60.;
      _cpMatrix.at(5).at(6).j() += 30.;
      
      // increase width
      for(std::size_t i = 1; i < _cpMatrix.size(); i++)
      {
        double s = (double)i/(double)(_cpMatrix.size()-1);
        double t = _cpMatrix.at(6).at(1).i() - 150.;
        double b = _cpMatrix.at(0).at(1).i();
        _cpMatrix.at(i).at(1).i() = (1.-s) * b + s * t;
        
        t = _cpMatrix.at(6).at(5).i() + 150.;
        b = _cpMatrix.at(0).at(5).i();
        _cpMatrix.at(i).at(5).i() = (1.-s) * b + s * t;
      }
    }
  }
}

//----------------------------------------------------------------

void Bezier::setControlPoint( const std::size_t row,
                              const std::size_t column,
                              const double x,
                              const double y )
{
  _cpMatrix.at(row).at(column).i() = x;
  _cpMatrix.at(row).at(column).j() = y;
  _cpMatrix.at(row).at(column).k() = 0.;
}

//----------------------------------------------------------------

void Bezier::focusCPs( const bool focusTop )
{
  double factor = 2./3.;
  std::size_t numInner = _cpMatrix.size()-2;
  
  // only process the inner control points
  for(std::size_t i = 1; i < _cpMatrix.size()-1; i++)
  {
    for(std::size_t j = 1; j < _cpMatrix.at(i).size()-1; j++)
    {
      Point3d top = _cpMatrix.back().at(j);
      Point3d bottom = _cpMatrix.front().at(j);
      double dist = fabs( top.j() - bottom.j() );
  
      // if focus on top then change the cps in such a way that
      // the upper cps are located in a band with bandwidth dist*(1.-factor)
      if( focusTop )
        bottom.j() += dist*factor;
      // if focus on bottom then change the cps in such a way that
      // the lower cps are located in a band with bandwidth dist*(1.-factor)
      else
        top.j() -= dist*factor;
      
      double s = (double)i/(double)(numInner);
      _cpMatrix.at(i).at(j).j() = (1.-s) * bottom.j() + s * top.j();
    }
  }  
}

//----------------------------------------------------------------

// Given u,v in parametric space, return x,y,z
Point3d Bezier::EvalCoord( double u, double v )
{
  // Make sure u and v between 0.0 and 1.0
  u = InBounds(u, 0.0, 1.0);
  v = InBounds(v, 0.0, 1.0);

  // Evaluate
  Point3d pos( 0., 0., 0. );
  for(std::size_t i = 0; i < _cpMatrix.size(); i++) 
    for(std::size_t j = 0; j < _cpMatrix.at(i).size(); j++)
    {
      double s = (double)Choose(_cpMatrix.size() - 1, i) * pow(v, i) * 
                 pow(1.0 - v, (int)(_cpMatrix.size() - 1 - i)) *
                 (double)Choose(_cpMatrix.at(i).size() - 1, j) * pow(u, j) * 
                 pow(1.0 - u, (int)(_cpMatrix.at(i).size() - 1 - j));
                 
      pos += s * _cpMatrix.at(i).at(j);
    }

  return pos;
}

//----------------------------------------------------------------

// Find binomial coefficient
int Bezier::Choose( unsigned int n, unsigned int k )
{
  static bool first = true;
  const std::size_t size = 7;
  static int choose[size][size];

  // Only calculate what we need
  if(n > size || k > n) {
    cerr << "Bezier::Choose:Error Can't compute " << n << " choose " << k 
                                                               << endl;
    return(0);
  }

  // If first time in, fill the table
  if(first) {
    first = false;
    for(unsigned int j = 0; j < size; j++)
      for(unsigned int i = j; i < size; i++)
        if(j == 0 || j == i)
          choose[i][j] = 1;
        else 
          choose[i][j] = choose[i-1][j] + choose[i-1][j-1];
  }

  // Just grab value from the table
  return(choose[n][k]);
}

#endif
