#ifndef SurfaceClass_H
#define SurfaceClass_H

#include "ModelHeader.h"
#include "ModelUtils.h"

/**
  @file   SurfaceClass.h
  @brief  Contains class for setting initial surfaces of model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

const double EPS = 0.0001;

class SurfaceClass
{
public:
  SurfaceClass(){}
  
  void init( const std::size_t lod,
             const std::string &lineageFileName,
             const bool exportLineage,
             const std::size_t timeSteps,
             const SurfaceType::type sType );
  
  void initModelBasedOnBezier( MyTissue &T,
                               const std::size_t cellNumber,
                               Surface &lateralRoot );
  
  void initLateralRootBasedOnBezier( MyTissue &T,
                                     const std::string &dataset,
                                     Surface &lateralRoot );
  
  void initLateralRootBasedOnRealData( MyTissue &T,
                                       RealSurface &lateralRoot,
                                       const std::string &dataset,
                                       const double surfaceScale,
                                       const bool useAutomaticContourPoints );
  
  void generateCell( MyTissue &T,
                     const std::pair<double, double> &start,
                     const std::pair<double, double> &length,
                     const std::size_t treeId,
                     Surface &lateralRoot );
  
  void generateCell( MyTissue &T,
                     const std::pair<double, double> &start,
                     const std::pair<double, double> &length,
                     const std::size_t treeId,
                     RealSurface &lateralRoot,
                     const std::vector<Point3d> &conPoints,
                     const bool oneCell = false );
  
  Point3d determinePos( const std::pair<double, double> &coord,
                        const std::vector<Point3d> &conPoints );
  
  void incrementTime(){ _time++; }
  void incrementCellID(){ _IDCounter++; }
  
  std::size_t getTime() const { return _time; }
  std::size_t getMaxTime() const { return _maxTime; }
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
  // max time step of data set
  std::size_t _maxTime;
  
  std::string _lineageFileName;
  
  bool _exportLineage;
  
  SurfaceType::type _sType;
};

#endif // SurfaceClass_H