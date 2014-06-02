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

// Load from file written by patch editor
void Bezier::Load(string bezFile)
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

// Copy and scale
void Bezier::Scale(Bezier &b, double scale)
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
void Bezier::Interpolate(Bezier &src1, Bezier &src2, 
                                     double scale1, double scale2, double s)
{
  if(src1.count != src2.count)
    cerr << "Bezier::Interpolate Error must have matching patch count" << endl;
  count = src1.count;

  s = InBounds(s, 0.0, 1.0);
  for(unsigned int i = 0; i < src1.count; i++)
    for(unsigned int j = 0; j < src1.count; j++) {
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
  
  for(unsigned int k = 0; k < count * count; k++) {
    for(int i = 0; i < 16; i++) {
      for(int j = 0; j < 3; j++)
        cout << *(ptsV + k*PATCHSIZE + i * 3 + j) << " ";
      cout << endl;
    }
    cout << endl;
  }
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

//----------------------------------------------------------------

// Given u,v in parametric space, return x,y,z
Point3d Bezier::EvalCoord(double u, double v)
{
  // Make sure u and v between 0.0 and 1.0
  u = InBounds(u, 0.0, 1.0);
  v = InBounds(v, 0.0, 1.0);

  unsigned int iu = 0;
  unsigned int iv = 0;

  // Scale so entire surface fits between 0 and 1
  if(count > 1) {
    u *= count;
    v *= count;

    iu = (unsigned int)u;
    iv = (unsigned int)v;

    if(iu == count) {
      u = 1.0;
      iu = count - 1;
    } else
      u = u - iu;
    if(iv == count) {
      v = 1.0;
      iv = count - 1;
    } else
      v = v - iv;
  }

  // Find correct patch
  float *pts = ptsV + (iv * MAXPATCHES + iu) * PATCHSIZE;

  // Evaluate
  double x = 0, y = 0, z = 0;
  for(unsigned int j = 0; j < PATCHPOINTS; j++)
    for(unsigned int i = 0; i < PATCHPOINTS; i++) {
      double s = (double)Choose(PATCHPOINTS - 1, i) * pow(u, i) * 
	                        pow(1.0 - u, (int)(PATCHPOINTS - 1 - i)) *
                 (double)Choose(PATCHPOINTS - 1, j) * pow(v, j) * 
		                pow(1.0 - v, (int)(PATCHPOINTS - 1 - j));
      x += s * *pts++;
      y += s * *pts++;
      z += s * *pts++;
    }
  Point3d p(x, y, z);
  return(p);
}
