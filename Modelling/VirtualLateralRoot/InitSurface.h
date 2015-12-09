#ifndef InitSurface_H
#define InitSurface_H

#include "ModelHeader.h"
#include "ModelUtils.h"
#include "SurfaceBaseClass.h"

/**
  @file   InitSurface.h
  @brief  Contains class for setting initial surfaces of model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

class InitSurface
{
public:
  InitSurface(){}
  
  void init( const std::size_t lod,
             const std::string &lineageFileName,
             const bool forceInitialSituation,
             const bool useAutomaticContourPoints,
             const double surfaceScale,
             const std::string &dataset );
  
  void initIdealizedCells( MyTissue &T,
                           const std::size_t founderCells,
                           SurfaceBaseClass &surface );
  
  void initRealDataCells( MyTissue &T, SurfaceBaseClass &surface );
  
  void initRadialCells( MyTissue &T, SurfaceBaseClass &surface );
  
  void initLateralRootBasedOnRealData( MyTissue &T, SurfaceBaseClass &surface );
  
  void generateCell( MyTissue &T,
                     const std::pair<double, double> &start,
                     const std::pair<double, double> &length,
                     const std::size_t treeId,
                     SurfaceBaseClass &surface,
                     const std::size_t addJunctionToWall = 4,
                     const int cellFile = 0 );
  
  Point3d determinePos( const double u, const double v );
  
  void incrementTime(){ _time++; }
  void incrementCellID(){ _IDCounter++; }
  
  std::size_t getTime() const { return _time; }
  std::size_t getCellID() const { return _IDCounter; }
private:
  
  void findNearestPointToMerge( MyTissue &T, junction &js );
  
  // sub division level for number of junctions of initial cells
  std::size_t _lod;
  // current counter for cell IDs
  std::size_t _IDCounter;
  // junction IDs for the initial junctions
  std::size_t _jIDCounter;
  // current time step
  std::size_t _time;
  
  std::string _lineageFileName;
  
  bool _forceInitialSituation;
  
  bool _useAutomaticContourPoints;
  
  // only valid for surface based on real data sets
  double _surfaceScale;
  
  std::string _dataset;
  
  bool _radialModel;
};

#endif // InitSurface_H