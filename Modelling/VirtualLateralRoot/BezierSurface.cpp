#include "BezierSurface.h"

#include <fstream>

/**
  @file   BezierSurface.cpp
  @brief  Class for handling a bezier surface
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

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

void BezierSurface::Load( const std::string &bezFile )
{
  std::ifstream bIn( bezFile.c_str() );
  if(!bIn)
  {
    std::cerr << "Bezier::Bezier:Error opening " << bezFile << std::endl;
    exit(1);
  }
  
  std::string line;
  while(bIn)
  {
    getline( bIn, line );
    if(line.substr(0, 9) == "BezPoints")
      break;
  }
    
  std::stringstream ss( line );
  std::size_t elements;
  std::string dummy;
  ss >> dummy >> elements;
  
  //std::cout << "elements: " << elements << std::endl;
  
  // skip next line
  getline( bIn, line );
  
  // initialize control points matrix
  _cpMatrix.resize( elements );
  for( std::size_t r=0;r<elements;r++ )
    _cpMatrix.at(r).resize( elements );
  
  double factor = 1000000.;
  
  for(std::size_t i = 0; i < elements; i++)
    for(std::size_t j = 0; j < elements; j++)
    {
      getline( bIn, line );
      std::stringstream ss( line );
      
      double x,y,z;
      ss >> dummy >> x >> y >> z;
      Point3d pos( x/factor, y/factor, 0. );
      _cpMatrix.at(i).at(j) = pos;
      //std::cout << "i:" << i << " j:" << j << " pos: " << pos << std::endl;
    }
    
  bIn.close();
}

//----------------------------------------------------------------

// Copy and scale
void BezierSurface::Scale( BezierSurface &b, double scale )
{
  std::vector< std::vector<Point3d> > mat = b._cpMatrix;
  for(std::size_t i = 0; i < _cpMatrix.size(); i++) 
    for(std::size_t j = 0; j < _cpMatrix.at(i).size(); j++)
      _cpMatrix.at(i).at(j) = mat.at(i).at(j) * scale;
}

//----------------------------------------------------------------

// Interpolate (plus scale) between two patches
void BezierSurface::Interpolate( BezierSurface &src1, BezierSurface &src2, 
                                 double scale1, double scale2, double s )
{
  s = InBounds(s, 0.0, 1.0);
  std::vector< std::vector<Point3d> > mat1 = src1._cpMatrix;
  std::vector< std::vector<Point3d> > mat2 = src2._cpMatrix;
  
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
void BezierSurface::Print()
{
  for(std::size_t i = 0; i < _cpMatrix.size(); i++)
  {
    for(std::size_t j = 0; j < _cpMatrix.at(i).size(); j++)
      std::cout << _cpMatrix.at(_cpMatrix.size()-1-i).at(j) << " ";
    
    std::cout << std::endl;
  }
}

//----------------------------------------------------------------

// Given u,v in parametric space, return x,y,z
Point3d BezierSurface::EvalCoord( double u, double v )
{
  // Make sure u and v between 0.0 and 1.0
  u = InBounds(u, 0.0, 1.0);
  v = InBounds(v, 0.0, 1.0);

  // Evaluate
  Point3d pos( 0., 0., 0. );
  for(std::size_t i = 0; i < _cpMatrix.size(); i++) 
    for(std::size_t j = 0; j < _cpMatrix.at(i).size(); j++)
    {
      double s = (double)binom(_cpMatrix.size() - 1, i) * pow(u, i) * 
                 pow(1.0 - u, (int)(_cpMatrix.size() - 1 - i)) *
                 (double)binom(_cpMatrix.at(i).size() - 1, j) * pow(v, j) * 
                 pow(1.0 - v, (int)(_cpMatrix.at(i).size() - 1 - j));
                 
      pos += s * _cpMatrix.at(i).at(j);
    }

  return pos;
}

//----------------------------------------------------------------

// Find binomial coefficient
int BezierSurface::binom( unsigned int n, unsigned int k )
{
  const std::size_t size = 8;
  if( n > size || k > n )
  {
    std::cerr << "Bezier::binom:Error Can't compute "
              << n << " choose " << k << std::endl;
    return(0);
  }

  static int binom[size][size];
  static bool first = true;
  
  if( first )
  {
    first = false;
    
    for(unsigned int j = 0; j < size; j++)
      for(unsigned int i = j; i < size; i++)
      {
        if( j == 0 || j == i )
          binom[i][j] = 1;
        else 
          binom[i][j] = binom[i-1][j] + binom[i-1][j-1];
      }
  }
  
  return binom[n][k];
}

//----------------------------------------------------------------
