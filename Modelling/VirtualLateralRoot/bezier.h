// bezier.h - Load Bezier from files output by interactive editor.
//
// Handles multiple patch files (advanced editor) and combines into 
// one surface parameterized by u, v, both between 0 and 1.
//
// Only patch info is used, not heading, etc. Order for multiple patch files
// is defined below, order info NOT read from file.
//

#ifndef __BEZIER_H__
#define __BEZIER_H__

#include <string>
#include <vector>

#include <util/parms.h>
#include <util/vector.h>

typedef util::Vector<3, double> Point3d;
typedef std::vector< std::vector<Point3d> > conpoi;

// Size of individual patches 4 x 4
const unsigned int PATCHPOINTS = 4;
const unsigned int PATCHSIZE = PATCHPOINTS * PATCHPOINTS * 3;

// MaxPatches is square of this 
const unsigned int MAXPATCHES = 4;
// Patch editor writes them this way.
const unsigned int PATCHMAP[MAXPATCHES * MAXPATCHES] = 
                      {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};

using namespace std;

class Bezier {
  public:
    Bezier() {};
    ~Bezier() {};

    // Load from file written by interactice editor.
    void Load( const std::string bezFile );
    
    // Load bezier surface with included growth tensor information from file
    void LoadGrowthBezier( const std::string bezFile,
                           const std::size_t numPatches,
                           const std::size_t numControlPoints );
    
    // Copy and/or scale
    void Scale( Bezier &b, double Scale );

    // Create new interpolated (and scaled) surface
    void Interpolate( Bezier &src1, Bezier &src2, 
                      double scale1, double scale2, double s );

    // Return x,y,z point from u,v parameters
    Point3d EvalCoord( double u, double v );

    // Print points (debugging)
    void Print();
    // Get points (debugging)
    float *Points();
    
    conpoi getControlPoints() const
    { return _cpMatrix; }
    
    void setCount( const unsigned int cnt )
    { count = cnt; }
    
    void focusCPs( const bool focusTop );
    
    void applyGrowthOnlyInHeight( const Bezier &surfaceE );

  private:
    int Choose( unsigned int n, unsigned int k );
    float ptsV[MAXPATCHES * MAXPATCHES * PATCHSIZE];

    // matrix storing unique control points
    conpoi _cpMatrix;
    
    // Note this is square root of count of # of patches
    unsigned int count;
};

#endif
