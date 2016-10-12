#include "SSimilarityMeasureUtil.hh"

/**
  @file   SSimilarityMeasureUtil.cc
  @brief  Contains namespace for similarity measure tools
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include <boost/foreach.hpp>
//#include <boost/numeric/conversion/bounds.hpp>
//#include <boost/limits.hpp>

#include "../similarityMeasure/STreeDivisionSequence.hh"

namespace SSimilarityMeasureUtil
{

// ---------------------------------------------------------------------

unsigned int determineLayerValue( const SLineageTree *cell,
                                  boost::shared_ptr<cellLayerVector> cellLayersPerTimeStep )
{
  unsigned int layerValue = 0;

  // for each layer
  BOOST_FOREACH( cellSet cSIL, cellLayersPerTimeStep->at(cell->timeStep-1) )
  {
    cellSet::const_iterator iter = cSIL.find( cell );
    if( iter == cSIL.end() )
      layerValue++;
    else
      break;
  }

  return layerValue;
}

// ---------------------------------------------------------------------

std::vector<int> determineAveragedCellCycleTimeSteps( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                                                      const std::size_t numCellCycles )
{
  std::vector<int> averagedCycleShapeTimes;
  std::vector<int> divisionPerCellCycle;
  averagedCycleShapeTimes.resize( numCellCycles, 0 );
  divisionPerCellCycle.resize( numCellCycles, 0 );

  for( SAdvancedLineageTreeCollection::const_iterator l = lineages->begin();
       l != lineages->end(); ++l )
  {
    for( SLineageTree::const_iterator iter = l->second->begin();
         iter != l->second->end(); ++iter )
    {
      if( iter->children.size() == 2 )
      {
        osg::Vec3 divPos( iter->getX(), iter->getY(), iter->getZ() );

        osg::Vec3 c1( iter->children[0]->getX(),
            iter->children[0]->getY(),
            iter->children[0]->getZ() );

        osg::Vec3 c2( iter->children[1]->getX(),
            iter->children[1]->getY(),
            iter->children[1]->getZ() );

        std::size_t cellCycle = iter->getCellCycleId();
        averagedCycleShapeTimes.at(cellCycle) += iter->timeStep-1;
        divisionPerCellCycle.at(cellCycle)++;
      }
    }
  }

  for( std::size_t i=0; i < numCellCycles; i++ )
    averagedCycleShapeTimes.at(i) /= divisionPerCellCycle.at(i);

  return averagedCycleShapeTimes;
}

//----------------------------------------------------------------------------

void determineCellLifeDuration( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                                const std::vector<osg::Vec3> &centresOfMass,
                                const std::map<int,int> &cellFiles,
                                const osg::Matrix &rotMat,
                                std::map<std::size_t, std::size_t> &lifeDuration,
                                std::size_t &min, std::size_t &max )
{
  min = 300;
  max = 0;

  // determine division sequence in order to focus life duration
  // analysis on these cells rising in the outer layer
  STreeDivisionSequence ds( lineages, false, rotMat, false, AXISTYPE::CELLS );
  std::map<std::size_t, std::string> seq = ds.getDivisionSequences();

  for( SAdvancedLineageTreeCollection::const_iterator l = lineages->begin();
       l != lineages->end(); ++l )
  {
    for( SLineageTree::const_iterator iter = l->second->begin();
         iter != l->second->end(); ++iter )
    {
      if( !iter->parent || iter->parent->children.size() == 2 )
      {
        const SLineageTree *tree = *iter;
        std::size_t life = 1;

        while( tree->children.size() == 1 )
        {
          life++;
          tree = tree->children.at(0);
        }

        // ignore life durations that end in a leaf node
        if( tree->children.size() == 0 )
          continue;

        osg::Vec3 pos( tree->getX(), tree->getY(), tree->getZ() );
        pos = pos - centresOfMass.at( tree->timeStep-1 );
        pos = pos * rotMat;

        // skip outliers seen from the side view to focus on cells in QC (quiescent center)
//        if( pos.x() < -50. || pos.x() > 50. )
//          continue;

        int cellFile = cellFiles.at( tree->treeId );

//        if( cellFile < -1 || cellFile > 1 )
//          continue;

        auto seqIter = seq.find( tree->cellId );
        if( seqIter != seq.end() )
        {
//          std::cout << "seq: " << seqIter->second << std::endl;

//          std::size_t len = seqIter->second.size();
//          if( len < 3 )
//            continue;

          std::size_t found = seqIter->second.find( "012" );

          if( found == std::string::npos )
            continue;

//          std::cout << "pass" << std::endl;
        }
        else
        {
          std::cout << "cellid: " << tree->cellId << " t: " << tree->timeStep-1 << std::endl;
          std::cout << "NOT FOUND" << std::endl;
        }

//        if( tree->getCellCycleId() > 4 )
//          continue;

//        if( life > 150 )
//          continue;

//        if( life < 30 )
//          continue;

        lifeDuration.insert( std::make_pair( iter->cellId, life ) );

        if( life < min )
          min = life;

        if( life >= max )
          max = life;
      }
    }
  }
}

//----------------------------------------------------------------------------

SphereC computeSphericalCoordinates( const osg::Vec3 &vec )
{
  // theta (inclination or latitude) is in [0, M_PI]
  // phi (azimuth or longitude) is in (-M_PI, M_PI]
  SphereC s;
  double x,y,z;
  x = vec.x();
  y = vec.y();
  z = vec.z();
  s.r = std::sqrt( x*x + y*y + z*z );
  s.theta = acosf( z/s.r );
  s.phi = atan2( y, x );

  return s;
}

// ---------------------------------------------------------------------

osg::Vec3 determineDivDirection( const SLineageTree *cell )
{
  osg::Vec3 parent( cell->parent->getX(),
                    cell->parent->getY(),
                    cell->parent->getZ() );
  osg::Vec3 child( cell->getX(),
                   cell->getY(),
                   cell->getZ() );

  osg::Vec3 divDir = child - parent;
  divDir.normalize();
  return divDir;
}

// ---------------------------------------------------------------------

osg::Vec3 determineDirectionVector( const SLineageTree *cell1,
                                    const SLineageTree *cell2 )
{
  osg::Vec3 c1( cell1->getX(), cell1->getY(), cell1->getZ() );
  osg::Vec3 c2( cell2->getX(), cell2->getY(), cell2->getZ() );

  osg::Vec3 dir = c2 - c1;
  dir.normalize();
  return dir;
}

// ---------------------------------------------------------------------

double computeDivisionAngle( const SLineageTree *cell )
{
  // go past the division of the provided cell
  const SLineageTree *lastDiv = cell->parent;

  // find the last division
  while( lastDiv->parent && lastDiv->parent->children.size() < 2 )
    lastDiv = lastDiv->parent;

  // there is no prior division
  if( lastDiv == cell->parent || !lastDiv->parent ||
      lastDiv->parent->children.size() < 2 )
  {
    return -1;
  }

  // compute angle
  SArray dir1( lastDiv->getX() - lastDiv->parent->getX(),
               lastDiv->getY() - lastDiv->parent->getY(),
               lastDiv->getZ() - lastDiv->parent->getZ() );
  SArray dir2( cell->getX() - cell->parent->getX(),
               cell->getY() - cell->parent->getY(),
               cell->getZ() - cell->parent->getZ() );

  dir1.normalize();
  dir2.normalize();

  return acos( dir1 * dir2 )*180./M_PI;
}

// ---------------------------------------------------------------------

// computation of several indexes to measure the similarity between
// cluster hierarchies: Rand index, adjusted Rand index, Jaccard index,
// Fowlkes-Mallows index (E. B. Fowlkes & C. L. Mallows 1983), and F_1 index
std::vector<double> computeIndexes( const std::size_t numObjects,
                                    const std::vector< std::set<std::size_t> > &cluster1,
                                    const std::vector< std::set<std::size_t> > &cluster2 )
{
  if( cluster1.size() != cluster2.size() )
  {
    std::cerr << "Cluster hierarchies must have the same number of clusters!" << std::endl;
    std::vector<double> temp;
    return temp;
  }

  std::size_t numClusters = cluster1.size();
  VMatrix M( numClusters );
  std::vector<double> rowSums;
  std::vector<double> columnSums;
  rowSums.resize( numClusters, 0. );
  columnSums.resize( numClusters, 0. );

  double T = -(double)numObjects;
  double P = -(double)numObjects;
  double Q = -(double)numObjects;
  double allBinoms = 0.;

  for( std::size_t i=0; i < numClusters; ++i )
    for( std::size_t j=0; j < numClusters; ++j )
    {
      // first determine the number of common objects in cluster1 and cluster2
      // and store the result in M(i,j)
      std::size_t commonObjects = 0;
      for( std::set<std::size_t>::const_iterator setIter = cluster1.at(i).begin();
           setIter != cluster1.at(i).end(); ++setIter )
      {
        if( cluster2.at(j).find( *setIter ) != cluster2.at(j).end() )
          commonObjects++;
      }

      M(i,j) = (double)commonObjects;
      T += M(i,j)*M(i,j);
      rowSums.at(i) += M(i,j);
      columnSums.at(j) += M(i,j);
      allBinoms += SSimilarityMeasureUtil::binomCoeff( commonObjects, 2 );
    }

  double rowBinoms = 0.;
  double columnBinoms = 0.;
  for( std::size_t i=0; i < numClusters; ++i )
  {
    P += rowSums.at(i)*rowSums.at(i);
    Q += columnSums.at(i)*columnSums.at(i);

    // temporal results for adjusted RI
    rowBinoms += SSimilarityMeasureUtil::binomCoeff( static_cast<std::size_t>(rowSums.at(i)), 2 );
    columnBinoms += SSimilarityMeasureUtil::binomCoeff( static_cast<std::size_t>(columnSums.at(i)), 2 );
  }

  // compute the adjusted Rand index
  double numBinom = SSimilarityMeasureUtil::binomCoeff( numObjects, 2 );
  double t1 = allBinoms - rowBinoms*columnBinoms/numBinom;
  double t2 = (rowBinoms+columnBinoms)/2. - rowBinoms*columnBinoms/numBinom;
  double aRI = t1/t2;

  double tp = T;
  double fp = P-T;
  double fn = Q-T;
  // number of objects that are in different clusters can be computed depending on the other results
  double tn = (double)numObjects*((double)numObjects-1.)/2. - tp - fp - fn;

  double precision = tp/(tp+fp);
  double recall = tp/(tp+fn);
  std::vector<double> indexSet;

  // Rand index
  //indexSet.push_back( (tp+tn)/(tp+fp+fn+tn) );
  // adjusted Rand index
  //indexSet.push_back( aRI );
  // Jaccard index
  indexSet.push_back( tp/(tp+fp+fn) );
  // Fowlkes-Mallows index
  indexSet.push_back( std::sqrt(precision*recall) );
  // F_1 index
  indexSet.push_back( 2.*precision*recall/(precision+recall) );
  return indexSet;
}

// ---------------------------------------------------------------------

double binomCoeff( const std::size_t n, std::size_t k )
{
  double res = 1.;
  if( k > n-k )
    k = n-k;

  for( auto i=0; i<k; ++i )
  {
    res *= (n-i);
    res /= (i+1);
  }

  return res;
}

// ---------------------------------------------------------------------

double computeSumOfSquaredErrors( const std::vector<std::vector<double> > &data,
                                  const std::vector< std::set<std::size_t> > &clusters )
{
  double sse = 0.;
  std::vector< std::vector<double> > mean;
  mean.resize( clusters.size() );

  // determine the mean vector of clusters
  for( std::size_t c = 0; c < clusters.size(); ++c )
  {
    // resize to number of elements
    mean.at(c).resize( data.at(0).size(), 0. );

    for( std::size_t e = 0; e <  mean.at(c).size(); ++e )
    {
      for( std::set<std::size_t>::iterator setIter = clusters.at(c).begin();
           setIter != clusters.at(c).end(); ++setIter )
      {
        std::vector<double> vec = data.at(*setIter);
        mean.at(c).at(e) += vec.at(e);
      }
      mean.at(c).at(e) /= clusters.at(c).size();
    }
  }

  for( std::size_t c = 0; c < clusters.size(); ++c )
  {
    double temp = 0.;
    for( std::set<std::size_t>::iterator setIter = clusters.at(c).begin();
         setIter != clusters.at(c).end(); ++setIter )
    {
      std::vector<double> vec = data.at(*setIter);
      double inner = 0.;
      // compute the Euclidean norm of vec and mean
      for( std::size_t e = 0; e < vec.size(); ++e )
        inner += (vec.at(e) - mean.at(c).at(e))*(vec.at(e) - mean.at(c).at(e));

      inner = std::sqrt(inner);
      temp += inner;
    }
    sse += temp;
  }

  return sse;
}

// ---------------------------------------------------------------------

int getDivisionNumber( const SLineageTree* const cell )
{
  int nDivs = 0;

  const SLineageTree *currCell = cell;

  while( currCell->parent != 0 )
  {
    if( currCell->children.size() > 1 )
      nDivs++;

    currCell = currCell->parent;
  }

  return nDivs;
}

// ---------------------------------------------------------------------

}
