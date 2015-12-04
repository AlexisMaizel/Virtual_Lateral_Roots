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
                     const bool useCombinedAreaRatio,
                     const bool useWallRatio,
                     const double divisionArea,
                     const double divisionAreaRatio,
                     const double divisionWallRatio,
                     const double LODThreshold,
                     const double avoidTrianglesThreshold,
                     const bool loadLastModel,
                     const double cellPinch,
                     const double cellMaxPinch,
                     const bool onlyGrowthInHeight,
                     Surface &VLRBezierSurface,
                     RealSurface &VLRDataPointSurface );
  
  void step_divisions( const std::size_t curTime );

  void step_cellWalls( const std::string &fileName );

  void step_tracking( const std::string &fileName );

  void step_growth( const double dt );
  
  void setAreaRatios( const bool loadLastModel,
                      const std::size_t initialSituationType );
  
private:

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
  bool _useCombinedAreaRatio;
  bool _useWallRatio;
  double _divisionArea;
  double _divisionAreaRatio;
  double _divisionWallRatio;
  double _LODThreshold;
  double _avoidTrianglesThreshold;
  bool _loadLastModel;
  double _cellPinch;
  double _cellMaxPinch;
  bool _onlyGrowthInHeight;
  Surface *_VLRBezierSurface;
  RealSurface *_VLRDataPointSurface;
  
  std::size_t _timeFourCellStage;
  std::size_t _timeSixCellStage;
  
  std::vector<std::pair<Point3d, Point3d> > _pcLines;
};

#endif // DivisionSetting_HH