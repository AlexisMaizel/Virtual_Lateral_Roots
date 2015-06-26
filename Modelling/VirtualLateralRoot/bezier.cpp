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
      
      double s = (double)j/(double)(numInner);
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
