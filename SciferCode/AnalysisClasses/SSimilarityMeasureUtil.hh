#ifndef SSimilarityMeasureUtil_HH
#define SSimilarityMeasureUtil_HH

/**
  @file   SSimilarityMeasureUtil.hh
  @brief  Contains namespace for similarity measure tools
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SSimilarityMeasureHeader.hh"

namespace SSimilarityMeasureUtil
{

unsigned int determineLayerValue( const SLineageTree *cell,
                                  boost::shared_ptr<cellLayerVector> cellLayersPerTimeStep );

std::vector<int> determineAveragedCellCycleTimeSteps( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                                                      const std::size_t numCellCycles );

void determineCellLifeDuration( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                                const std::vector<osg::Vec3> &centresOfMass,
                                const std::map<int,int> &cellFiles,
                                const osg::Matrix &rotMat,
                                std::map<std::size_t, std::size_t> &lifeDuration,
                                std::size_t &min, std::size_t &max );

SphereC computeSphericalCoordinates( const osg::Vec3 &vec );

osg::Vec3 determineDivDirection( const SLineageTree *cell );

osg::Vec3 determineDirectionVector( const SLineageTree *cell1,
                                    const SLineageTree *cell2 );

double computeDivisionAngle( const SLineageTree *cell );

int getDivisionNumber( const SLineageTree* const cell );

double binomCoeff( const std::size_t n, std::size_t k );

std::vector<double> computeIndexes( const std::size_t numObjects,
                                    const std::vector< std::set<std::size_t> > &cluster1,
                                    const std::vector< std::set<std::size_t> > &cluster2 );

double computeSumOfSquaredErrors( const std::vector<std::vector<double> > &data,
                                  const std::vector< std::set<std::size_t> > &clusters );

}

#endif // SSimilarityMeasureUtil_HH
