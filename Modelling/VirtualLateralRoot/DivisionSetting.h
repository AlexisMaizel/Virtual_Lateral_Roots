#ifndef DivisionSetting_HH
#define DivisionSetting_HH

#include "ModelHeader.h"

/**
  @file   DivisionSetting.h
  @brief  Class for defining division steps in model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

class DivisionSetting
{
public:
  
  DivisionSetting( MyTissue &T,
                     const std::size_t initialSituationType,
                     const std::string &divisionType,
                     const std::size_t timeDelay,
                     const double firstDivisionsAreaRatio,
                     const double secondDivisionsAreaRatio,
                     const bool useAlternativeDT,
                     const bool accurateCenterOfMass,
                     const double probabilityOfDecussationDivision,
                     const bool useAreaRatio,
                     const double equalAreaRatio,
                     const bool useCombinedAreaRatio,
                     const bool useWallRatio,
                     const double divisionArea,
                     const double divisionAreaRatio,
                     const double divisionWallRatio,
                     const double LODThreshold,
                     const double avoidTrianglesThreshold,
                     const bool loadLastModel,
                     const bool onlyGrowthInHeight,
                     const std::size_t surfaceType,
                     Surface &VLRBezierSurface,
                     RealSurface &VLRDataPointSurface );
  
  void clear();

  void step_cellWalls( const std::string &fileName );

  void step_tracking( const std::string &fileName );

  void step_growth( const double dt );
  
  void setAreaRatios( const bool loadLastModel,
                      const std::size_t initialSituationType );
  
  const std::vector<std::size_t> &getRandomChoices() const
  { return _randomChoices; }

  void setRandomChoices( const std::vector<std::size_t> &randomChoices )
  { _randomChoices = randomChoices; _choiceCounter = 0; }
  
  void setCellDivisionSettings( const cell &c,
                                const std::size_t curTime,
                                const bool fixedCenterPos = false );
  
  MyTissue::division_data setDivisionPoints( const cell& c );
  
  MyTissue::division_data getEnergyDivisionData( const cell& c, bool &empty );
  
  MyTissue::division_data getBessonDumaisDivisionData( const cell& c, bool &empty );
  
  MyTissue::division_data getRandomDivisionData( const cell& c, bool &empty,
                                                 const bool equalArea );
  
  double generateNoiseInRange( const double start, const double end,
                               const unsigned int steps );

private:
  
  MyTissue *_T;
  std::size_t _initialSituationType;
  std::string _divisionType;
  std::size_t _timeDelay;
  double _firstDivisionsAreaRatio;
  double _secondDivisionsAreaRatio;
  bool _useAlternativeDT;
  bool _accurateCenterOfMass;
  double _probabilityOfDecussationDivision;
  bool _useAreaRatio;
  double _equalAreaRatio;
  bool _useCombinedAreaRatio;
  bool _useWallRatio;
  double _divisionArea;
  double _divisionAreaRatio;
  double _divisionWallRatio;
  double _LODThreshold;
  double _avoidTrianglesThreshold;
  bool _loadLastModel;
  bool _onlyGrowthInHeight;
  std::size_t _surfaceType;
  
  std::vector<std::size_t> _randomChoices;
  std::size_t _choiceCounter;
  
  Surface *_VLRBezierSurface;
  RealSurface *_VLRDataPointSurface;
  
  std::vector<std::pair<Point3d, Point3d> > _pcLines;
};

#endif // DivisionSetting_HH