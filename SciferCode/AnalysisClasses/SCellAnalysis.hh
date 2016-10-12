#ifndef SCELLANALYSIS_HH
#define SCELLANALYSIS_HH

/**
  @file   SCellAnalysis.hh
  @brief  Contains class for analyzing cell growth and movement directions
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SVoronoiDiagram.hh"

struct SLengthAnalysis
{
  std::size_t cellId;
  std::size_t timeStep;
  std::size_t cellCycle;
  std::pair<double,double> daughterDivisionAngles;
  std::pair<osg::Vec3,osg::Vec3> daughterDivisionDirections;
  std::vector<osg::Vec3> neighborDisplacementVectors;
  osg::Vec3 principalGrowthVector;
  osg::Vec3 mostDisplacementNeighborVector;
  std::pair<double,double> correlationPrincipalGrowthDivisionAngles;
  std::pair<double,double> correlationMostDisplacementGrowthDivisionAngles;
  std::size_t numberOfNeighbors;
};

struct SSingleNeighborDisplacement
{
  std::size_t cellId;
  std::size_t startTime;
  std::size_t cellCycle;
  std::map<std::size_t,osg::Vec3> singleDisplacements;
  osg::Vec3 principalDisplacement;
};

class SCellAnalysis
{
public:
  SCellAnalysis( boost::shared_ptr<SVoronoiDiagram> voronoi,
                 boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                 osg::ref_ptr<SSliderCallback> sliderCallback,
                 osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                 const bool firstPosFixed,
                 const bool useCycleSlider );

  void lengthAnalysis();

  //void checkCompleteCellCycle( const SLineageTree* node );

  void checkCellCycle( const SLineageTree* node );

  osg::Vec3 determinePrincipalDirectionVector( const std::vector<osg::Vec3> &neighborDisplacementVectors );

  double computeDivisionAngle( const SLineageTree *cell );

  void drawArrows( const osg::Vec3 &start,
                   const osg::Vec3 &end,
                   const std::size_t renderTimeStep,
                   const std::size_t arrowType );

//  void drawArrows( const std::vector<osg::Vec3> &nodePositions,
//                   const std::vector<osg::Vec3> &nNodePositions,
//                   const int startTimeStep,
//                   osg::Vec3 &startMeanPos,
//                   osg::Vec3 &endMeanPos,
//                   const std::size_t renderStep,
//                   const bool drawOnlyMean );

//  void drawDiffArrows( const osg::Vec3 &startPos,
//                       const osg::Vec3 &endPos,
//                       const osg::Vec3 &nStartPos,
//                       const osg::Vec3 &nEndPos,
//                       const osg::Vec3 &startMeanPos,
//                       const osg::Vec3 &endMeanPos,
//                       const std::size_t renderStep );

  void writeAnalysisResults( const std::string &fileName );

  osg::Vec3 getDivisionDirectionFromDaughterCells( const SLineageTree *cell );

  double computeAngleBetweenThisAndPriorDivision( const SLineageTree *cell );

  double collinearConversion( const double angle );

  void writeLengthParameters( const std::string &fileName );

  void writeSingleDisplacements( const std::string &fileName );

private:
  boost::shared_ptr<SVoronoiDiagram> _voronoi;
  boost::shared_ptr<SAdvancedLineageTreeCollection> _lineages;

  std::vector<osg::ref_ptr<osg::Group> > _analysisGroupPerTimeStep;

  std::vector<osg::ref_ptr<osg::Group> > _analysisGroupPerCellCycle;

  std::vector<SSingleNeighborDisplacement> _allSingleDisplacements;

  osg::Vec4 _neighborColor;

  osg::Vec4 _divisionColor;

  osg::Vec4 _principalColor;

  osg::Vec4 _realPrincipalColor;

  std::vector<SLengthAnalysis> _lengthAnalysisStats;

  /// compute the angles and draw the cyclinders for the
  /// first fixed position of the current node
  bool _firstPosFixed;

  /// use cycle slider callback or else time slider callback
  bool _useCycleSlider;

  /// radius of rendered arrows
  double _arrowRadius;
};

#endif // SCELLANALYSIS_HH
