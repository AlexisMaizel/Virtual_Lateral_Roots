#include "SCellAnalysis.hh"

/**
  @file   SCellAnalysis.cc
  @brief  Contains class for analyzing cell growth and movement directions
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SSimilarityMeasureGraphics.hh"
#include "SInteractiveLengthAnalysis.hh"

#include <fstream>

#include <osg/io_utils>

// ---------------------------------------------------------------------

SCellAnalysis::SCellAnalysis( boost::shared_ptr<SVoronoiDiagram> voronoi,
                              boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                              osg::ref_ptr<SSliderCallback> sliderCallback,
                              osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                              const bool firstPosFixed,
                              const bool useCycleSlider )
  : _voronoi( voronoi ),
    _lineages( lineages ),
    _neighborColor( osg::Vec4( 0., 1., 1., 1. ) ),
    _divisionColor( osg::Vec4( 0., 0. ,1., 1. ) ),
    _principalColor( osg::Vec4( 1., 0. ,0., 1. ) ),
    _realPrincipalColor( osg::Vec4( 1., 1. ,0., 1. ) ),
    _firstPosFixed( firstPosFixed ),
    _useCycleSlider( useCycleSlider ),
    _arrowRadius( 1.5 )
{
  // initialize callback for arrows
  osg::ref_ptr<SInteractiveLengthAnalysis> ila = new SInteractiveLengthAnalysis();
  ila->addUpdateCallback( sliderCallback );
  ila->setName( "Cell Analysis" );
  addToThisGroup->addChild( ila );

  std::vector<std::string> groupNames;
  groupNames.push_back( "Neighbor Arrows" );
  groupNames.push_back( "Division Arrows" );
  groupNames.push_back( "Principal Arrows" );
  groupNames.push_back( "TimeStep-based Principal Arrows" );

  if( _useCycleSlider )
  {
    std::size_t maxCellCycle = 0;

    // first get maximum cell cycle amount
    for( SAdvancedLineageTreeCollection::const_iterator l = _lineages->begin();
         l != _lineages->end(); ++l )
    {
      SLineageTree *tree = (*_lineages)[l->first];

      for( SLineageTree::iterator iter = tree->begin();
           iter != tree->end(); ++iter )
      {
        // only check the number of cell cycle for the leaves
        if( iter->children.size() == 0 )
        {
          std::size_t cc = iter->getCellCycleId();

          if( cc > maxCellCycle )
            maxCellCycle = cc;
        }
      }
    }

    _analysisGroupPerCellCycle.resize( maxCellCycle );
    for( std::size_t c = 0; c < maxCellCycle; c++ )
    {
      _analysisGroupPerCellCycle.at(c) = new osg::Group;

      for( std::size_t i = 0; i < groupNames.size(); i++ )
      {
        osg::ref_ptr<osg::Group> group = new osg::Group;
        group->setName( groupNames[i] );
        _analysisGroupPerCellCycle.at(c)->addChild( group );
      }
    }
  }
  else
  {
    _analysisGroupPerTimeStep.resize( _lineages->getOrComputeNumberOfTimesteps() );
    for( std::size_t t = 0; t < _analysisGroupPerTimeStep.size(); t++ )
    {
      _analysisGroupPerTimeStep.at(t) = new osg::Group;

      for( std::size_t i = 0; i < groupNames.size(); i++ )
      {
        osg::ref_ptr<osg::Group> group = new osg::Group;
        group->setName( groupNames[i] );
        _analysisGroupPerTimeStep.at(t)->addChild( group );
      }
    }
  }

  // add the initial groups for the interactive rendering
  for( std::size_t i = 0; i < groupNames.size(); i++ )
  {
    osg::ref_ptr<osg::Group> group = new osg::Group;
    group->setName( groupNames[i] );
    ila->addChild( group );
  }

  // perform length analysis
  this->lengthAnalysis();

  // afterwards set all osg groups for the callback
  if( _useCycleSlider )
  {
    for( std::size_t c = 0; c < _analysisGroupPerCellCycle.size(); c++ )
    {
      for( std::size_t a = 0; a < groupNames.size(); a++ )
        ila->addCycleArrows( c+1, a, _analysisGroupPerCellCycle.at(c)->getChild(a)->asGroup() );
    }
  }
  else
  {
    for( std::size_t t = 0; t < _analysisGroupPerTimeStep.size(); t++ )
    {
      for( std::size_t a = 0; a < groupNames.size(); a++ )
        ila->addArrows( t+1, a, _analysisGroupPerTimeStep.at(t)->getChild(a)->asGroup() );
    }
  }
}

// ---------------------------------------------------------------------

void SCellAnalysis::lengthAnalysis()
{
  // loop over all lineages
  for( SAdvancedLineageTreeCollection::const_iterator l = _lineages->begin();
       l != _lineages->end(); ++l )
  {
    SLineageTree *tree = (*_lineages)[l->first];

    // loop over all nodes of the tree
    // only start the length analysis for daughter cells
    // and the root node
    for( SLineageTree::iterator treeIter = tree->begin();
         treeIter != tree->end(); ++treeIter )
    {
      // root node
      if( !treeIter->parent )
        this->checkCellCycle( *treeIter );
      // daughter cells
      else if( treeIter->parent->children.size() == 2 )
        this->checkCellCycle( *treeIter );
    }
  }
}

// ---------------------------------------------------------------------

//void SCellAnalysis::checkCompleteCellCycle( const SLineageTree* node )
//{
//  // start time step
//  int startTime = node->timeStep;

//  // start positions of current node during the cell cycle
//  std::vector<osg::Vec3> nodePositions;

//  // store also the positions of the children of the last node
//  // which is a division node or in the case of a leaf do not
//  // draw the division arrows
//  osg::Vec3 c1Pos, c2Pos;

//  // store the cell cycle value of the last division node
//  std::size_t cellCycle;

//  // get end time step at division node
//  int endTime = startTime;
//  const SLineageTree* temp = node;
//  while( temp->children.size() == 1 )
//  {
//    nodePositions.push_back( osg::Vec3( temp->getX(), temp->getY(), temp->getZ() ) );
//    temp = temp->children[0];
//    endTime++;
//  }

//  // if the last node is a leaf node, then we ignore the analysis for now
//  if( temp->children.size() == 0 )
//    return;
//  // else if it is a division node then add the final position of it and also store the
//  // positions of its daughter cells for computing the division directions
//  else if( temp->children.size() == 2 )
//  {
//    nodePositions.push_back( osg::Vec3( temp->getX(), temp->getY(), temp->getZ() ) );
//    cellCycle = temp->getCellCycleId();
//    c1Pos = osg::Vec3( temp->children[0]->getX(), temp->children[0]->getY(), temp->children[0]->getZ() );
//    c2Pos = osg::Vec3( temp->children[1]->getX(), temp->children[1]->getY(), temp->children[1]->getZ() );
//  }

//  // get the direction between the daughter cells in order
//  // to generate only one angle describing the division
//  osg::Vec3 divisionDirection = c2Pos - c1Pos;
//  divisionDirection.normalize();

//  //std::cout << "node: " << node->cellId << " " << node->timeStep << std::endl;

//  // get neighbors of current node
//  std::set<const SLineageTree*> neighbors =
//      _voronoi->getNeighbors( node->cellId, node->timeStep );

//  // exclude these neighbors for which a division occurs while the current
//  // node has not divided yet
//  for( std::set<const SLineageTree*>::const_iterator setIter = neighbors.begin();
//       setIter != neighbors.end(); )
//  {
//    int nEndTime = startTime;
//    SLineageTree* temp = new SLineageTree( **setIter );
//    while( temp->children.size() == 1 )
//    {
//      nEndTime++;
//      temp = temp->children[0];
//    }

//    delete temp;

//    // if this neighbor node divides or ends earlier than
//    // the considered node then we ignore it at the moment
//    if( nEndTime < endTime )
//      neighbors.erase( setIter++ );
//    else
//      ++setIter;
//  }

//  std::size_t renderStep;

//  if( _useCycleSlider )
//    renderStep = cellCycle;
//  else
//    renderStep = endTime-1;

//  // store the growth angles for the greatest displacement
//  // between cell c and neighbors n_i
//  double greatestDisplacement = -1.;
//  double growthAngle;

//  // also store the average of all growth directions of each
//  // neighbor
//  osg::Vec3 averageGrowthDir( 0., 0., 0. );

//  // loop over all remaining neighbor nodes
//  for( std::set<const SLineageTree*>::const_iterator setIter = neighbors.begin();
//       setIter != neighbors.end(); ++setIter )
//  {
//    //std::cout << "node: " << (*setIter)->cellId << " " << (*setIter)->timeStep << std::endl;

//    std::size_t index = 1;
//    std::vector<osg::Vec3> nNodePositions;

//    // get neighbor positions and draw arrows
//    for( SLineageTree::const_iterator treeIter = (*setIter)->begin();
//         treeIter != (*setIter)->end(); ++treeIter, index++ )
//    {
//      // only draw the arrows until the division of the
//      // current considered node
//      if( index > nodePositions.size() )
//        break;

//      osg::Vec3 nPos( treeIter->getX(), treeIter->getY(), treeIter->getZ() );
//      nNodePositions.push_back( nPos );
//    }

//    osg::Vec3 startMeanPos;
//    osg::Vec3 endMeanPos;

//    // draw distance/length arrows
//    this->drawArrows( nodePositions, nNodePositions, startTime-1,
//                      startMeanPos, endMeanPos, renderStep, true );

//    // and afterwards the difference arrows along the same direction
//    // vector from the previous computed mean start and end positions
////    this->drawDiffArrows( nodePositions.front(), nodePositions.back(),
////                          nNodePositions.front(), nNodePositions.back(),
////                          startMeanPos, endMeanPos, renderStep );

//    // store the stats of the current division/length analysis
//    // but only for division nodes
//    SLengthAnalysis lAnalysis;
//    lAnalysis.cellId = node->cellId;
//    lAnalysis.timeStep = endTime;

//    // get direction fo growth and add value to the average
//    osg::Vec3 dirGrowth = endMeanPos - startMeanPos;
//    averageGrowthDir += dirGrowth;
//    dirGrowth.normalize();

//    // compute the angle between the direction vector of the two daughter
//    // cells and the growth direction
//    double angle = acos( dirGrowth * divisionDirection );
//    angle = angle * 180./M_PI;
//    lAnalysis.growthDivisionAngle = angle;

//    // compute the angle between the previous and the current division direction
//    lAnalysis.consecutiveDivisionAngle = this->computeAngleBetweenThisAndPriorDivision( temp );

//    // store the cell cycle of the last division node
//    lAnalysis.cellCycle = cellCycle+1;

//    lAnalysis.startLength = (nNodePositions.front() - nodePositions.front() ).length();

//    // compute the length development for the first
//    // fixed cell position of the current cell
//    if( _firstPosFixed )
//      lAnalysis.endLength = (nNodePositions.back() - nodePositions.front() ).length();
//    else
//      lAnalysis.endLength = (nNodePositions.back() - nodePositions.back() ).length();

//    double length = std::fabs( lAnalysis.endLength - lAnalysis.startLength );

//    if( length > greatestDisplacement )
//    {
//      greatestDisplacement = length;
//      growthAngle = lAnalysis.growthDivisionAngle;
//    }

//    _lengthAnalysisStats.push_back( lAnalysis );
//  }

//  if( neighbors.size() > 0 )
//  {
//    // store only one unique analysis structure for each occuring division node
//    SLengthAnalysis uniqueLengthAnalysis;
//    uniqueLengthAnalysis.cellId = node->cellId;
//    uniqueLengthAnalysis.timeStep = endTime;
//    uniqueLengthAnalysis.cellCycle = cellCycle+1;
//    uniqueLengthAnalysis.consecutiveDivisionAngle = this->computeAngleBetweenThisAndPriorDivision( temp );
//    // only store the angle for which the greatest displacement occurs
//    uniqueLengthAnalysis.growthDivisionAngle = growthAngle;

//    averageGrowthDir /= neighbors.size();
//    averageGrowthDir.normalize();
//    uniqueLengthAnalysis.averageGrowthDir = averageGrowthDir;
//    double angle = acos( averageGrowthDir * divisionDirection );
//    angle = angle * 180./M_PI;
//    uniqueLengthAnalysis.averageGrowthDivisionAngle = angle;

//    _uniqueCellLengthAnalysisStats.push_back( uniqueLengthAnalysis );
//  }
//}

// ---------------------------------------------------------------------

void SCellAnalysis::checkCellCycle( const SLineageTree* node )
{
  // start time step
  std::size_t startTime = node->timeStep;

  // only used for considering displacements between
  // consecutive time steps
  std::vector<osg::Vec3> nodePositions;

  // start position of current node
  osg::Vec3 startPosCell( node->getX(), node->getY(), node->getZ() );
  osg::Vec3 endPosCell;

  // store the cell cycle value of the last division node
  std::size_t cellCycle;

  // get end time step at division node
  std::size_t endTime = startTime;
  const SLineageTree* temp = node;

  //std::cout << "id: " << node->cellId << " t: " << node->timeStep << std::endl;

  while( temp->children.size() == 1 )
  {
    // store each position
    nodePositions.push_back( osg::Vec3( temp->getX(), temp->getY(), temp->getZ() ) );

    //std::cout << "pos: " << nodePositions.back() << std::endl;

    temp = temp->children[0];
    endTime++;
  }

  // if the last node is a leaf node, then we ignore the analysis for now
  if( temp->children.size() == 0 )
    return;
  // else if it is a division node then store the final position of it
  else if( temp->children.size() == 2 )
  {
    cellCycle = temp->getCellCycleId();
    endPosCell = osg::Vec3( temp->getX(), temp->getY(), temp->getZ() );
    nodePositions.push_back( endPosCell );
  }

  // get neighbors of current node
  std::set<const SLineageTree*> neighbors =
      _voronoi->getNeighbors( node->cellId, node->timeStep );

  // exclude these neighbors for which a division occurs while the current
  // node has not divided yet
  for( std::set<const SLineageTree*>::const_iterator setIter = neighbors.begin();
       setIter != neighbors.end(); )
  {
    std::size_t nEndTime = startTime;
    SLineageTree* temp = new SLineageTree( **setIter );
    while( temp->children.size() == 1 )
    {
      nEndTime++;
      temp = temp->children[0];
    }

    delete temp;

    // if this neighbor node divides or ends earlier than
    // the considered node then we ignore it at the moment
    if( nEndTime < endTime )
      neighbors.erase( setIter++ );
    else
      ++setIter;
  }

  // if no neighbors exists anymore then return
  if( neighbors.size() == 0 )
    return;

  std::size_t renderTimeStep;

  if( _useCycleSlider )
    renderTimeStep = cellCycle;
  else
    renderTimeStep = endTime-1;

  std::vector<osg::Vec3> startPosNeighbors;
  std::vector<osg::Vec3> endPosNeighbors;

  // store the stats of the current division/length analysis
  // but only for division nodes
  SLengthAnalysis lAnalysis;
  lAnalysis.cellId = node->cellId;
  lAnalysis.timeStep = endTime;
  // store the cell cycle of the last division node
  lAnalysis.cellCycle = cellCycle+1;
  // store the number of relevant neighbors
  lAnalysis.numberOfNeighbors = neighbors.size();
  std::vector<osg::Vec3> neighborDisplacementVectors;
  std::vector<osg::Vec3> neighborCombinedDisplacementVectors;

  // loop over all remaining neighbor nodes
  for( std::set<const SLineageTree*>::const_iterator setIter = neighbors.begin();
       setIter != neighbors.end(); ++setIter )
  {
    // only used for considering displacements between
    // consecutive time steps
    std::vector<osg::Vec3> neighborPositions;

    // get neighbor positions and draw arrows
    for( SLineageTree::const_iterator treeIter = (*setIter)->begin();
         treeIter != (*setIter)->end(); ++treeIter )
    {
      osg::Vec3 currentPos( treeIter->getX(), treeIter->getY(), treeIter->getZ() );

      // store each neighbor position
      neighborPositions.push_back( currentPos );

      // store first position of neighbor node
      if( treeIter->timeStep == startTime )
        startPosNeighbors.push_back( currentPos );

      // get the last position of the current neighbor
      if( treeIter->timeStep == endTime )
      {
        endPosNeighbors.push_back( currentPos );
        break;
      }
    }

    // determine the displacement vectors between the current division cell
    // and the the current neighbor cell
    osg::Vec3 startDiff = startPosCell - startPosNeighbors.back();
    osg::Vec3 endDiff = endPosCell - endPosNeighbors.back();

    osg::Vec3 displacementVector = endDiff - startDiff;
    neighborDisplacementVectors.push_back( displacementVector );

    // neighbor arrows
    //this->drawArrows( startPosNeighbors.back(), startPosCell, startTime, 0 );
    //this->drawArrows( endPosNeighbors.back(), endPosCell, endTime, 0 );

    // each displacement arrow of a neighbor
    //this->drawArrows( endPosCell, endPosCell + displacementVector, renderTimeStep, 0 );


    // instance of single displacement type
    SSingleNeighborDisplacement snd;
    snd.cellCycle = cellCycle + 1;
    snd.startTime = startTime;
    snd.cellId = node->cellId;
    snd.principalDisplacement = displacementVector;
    std::map<std::size_t,osg::Vec3> singleDisplacementMap;

    // after traversing all time steps of the current considered neighbor cell
    // compute the single displacement vectors between two consecutive time steps
    osg::Vec3 combinedDisplacementVector( 0., 0., 0. );
    // and combine them to one principal displacement vector
    for( std::size_t p = 1; p < neighborPositions.size(); p++ )
    {
      // determine the displacement vectors between the current cell position
      // and the the current neighbor cell position
      osg::Vec3 start = nodePositions.at(p-1) - neighborPositions.at(p-1);
      osg::Vec3 end = nodePositions.at(p) - neighborPositions.at(p);

      osg::Vec3 singleDisplacement = end - start;
      // testing of normalization of vectors
      //singleDisplacement.normalize();
      //singleDisplacement *= 10.;

      singleDisplacementMap.insert( std::make_pair( startTime + p - 1, singleDisplacement ) );

//      if( !_useCycleSlider )
//      {
//        this->drawArrows( nodePositions.at(p), nodePositions.at(p) + singleDisplacement, startTime + p - 1, 0 );
//        // each displacement arrow of a neighbor
//        this->drawArrows( nodePositions.at(p), nodePositions.at(p) + displacementVector, startTime + p - 1, 2 );
//      }
//      else
//      {
//        this->drawArrows( endPosCell, endPosCell + singleDisplacement, renderTimeStep, 0 );
//      }

      combinedDisplacementVector += singleDisplacement;

      osg::Vec3 a,b;
      a = singleDisplacement;
      b = displacementVector;
      a.normalize();
      b.normalize();
      double angle = std::acos( a * b );
      angle = angle * 180./M_PI;

//      if( angle > 90. )
//        std::cout << "diff" << std::endl;
    }

    snd.singleDisplacements = singleDisplacementMap;
    _allSingleDisplacements.push_back( snd );

    // TESTING by selecting the mean of all positions of the current
    // considered node and its current neighbor node
//    osg::Vec3 startPos( 0., 0., 0. );
//    osg::Vec3 endPos( 0., 0., 0. );
//    for( std::size_t p = 0; p < neighborPositions.size(); p++ )
//    {
//      startPos += neighborPositions.at(p);
//      endPos += nodePositions.at(p);
//    }

//    startPos /= neighborPositions.size();
//    endPos /= neighborPositions.size();

//    osg::Vec3 displacement = end - start;

    //combinedDisplacementVector /= neighborPositions.size();
    neighborCombinedDisplacementVectors.push_back( combinedDisplacementVector );
  }

  lAnalysis.neighborDisplacementVectors = neighborDisplacementVectors;

  // get the direction vector with the most displacement
  double mostDisplacementLength = -1.;
  std::size_t index;
  for( std::size_t n = 0; n < neighborDisplacementVectors.size(); n++ )
  {
    if( neighborDisplacementVectors.at(n).length() > mostDisplacementLength )
    {
      mostDisplacementLength = neighborDisplacementVectors.at(n).length();
      index = n;
    }
  }
  lAnalysis.mostDisplacementNeighborVector = neighborDisplacementVectors.at(index);

  // draw arrow with the most displacement
  osg::Vec3 startDiff = startPosCell - startPosNeighbors.at(index);
  osg::Vec3 endDiff = endPosCell - endPosNeighbors.at(index);
  osg::Vec3 displacementVector = endDiff - startDiff;
  //this->drawArrows( endPosCell, endPosCell + displacementVector, renderTimeStep, 2 );

  // get the principal direction vector of all displacement vectors
  lAnalysis.principalGrowthVector = this->determinePrincipalDirectionVector( neighborDisplacementVectors );

  //std::cout << "vector1: " << lAnalysis.principalGrowthVector << std::endl;

  // principal arrow
  this->drawArrows( endPosCell - lAnalysis.principalGrowthVector, endPosCell, renderTimeStep, 2 );

  // real principal arrow which considers not only the first and prelast position before
  // the next division but all positions in between
  osg::Vec3 realPrincipalGrowthVector = this->determinePrincipalDirectionVector( neighborCombinedDisplacementVectors );
  this->drawArrows( endPosCell - realPrincipalGrowthVector, endPosCell, renderTimeStep, 3 );

  //std::cout << "vector2: " << realPrincipalGrowthVector << std::endl;

  // at last determine the angles related to the division
  std::pair<double,double> daughterDivisionAngles;
  std::pair<osg::Vec3,osg::Vec3> daughterDivisionDirections;
  osg::Vec3 divisionDirection( temp->children[0]->getX() - temp->getX(),
                               temp->children[0]->getY() - temp->getY(),
                               temp->children[0]->getZ() - temp->getZ() );
  daughterDivisionAngles.first = this->computeDivisionAngle( temp->children[0] );
  daughterDivisionDirections.first = divisionDirection;
  daughterDivisionAngles.second = this->computeDivisionAngle( temp->children[1] );
  divisionDirection = osg::Vec3( temp->children[1]->getX() - temp->getX(),
      temp->children[1]->getY() - temp->getY(),
      temp->children[1]->getZ() - temp->getZ() );
  daughterDivisionDirections.second = divisionDirection;
  lAnalysis.daughterDivisionAngles = daughterDivisionAngles;
  lAnalysis.daughterDivisionDirections = daughterDivisionDirections;

  // division arrows
  this->drawArrows( endPosCell, endPosCell + daughterDivisionDirections.first, renderTimeStep, 1 );
  this->drawArrows( endPosCell, endPosCell + daughterDivisionDirections.second, renderTimeStep, 1 );

  // and the correlation of the principal growth direction with the division
  // as well as the most displacement direction with the division direction
  std::pair<double,double> correlationMostDisplacementGrowthDivisionAngles;
  osg::Vec3 dir1,dir2;
  dir1 = lAnalysis.mostDisplacementNeighborVector;
  dir2 = lAnalysis.daughterDivisionDirections.first;
  dir1.normalize();
  dir2.normalize();
  correlationMostDisplacementGrowthDivisionAngles.first = std::acos( dir1*dir2 )*180./M_PI;
  dir2 = lAnalysis.daughterDivisionDirections.second;
  dir2.normalize();
  correlationMostDisplacementGrowthDivisionAngles.second = std::acos( dir1*dir2 )*180./M_PI;
  lAnalysis.correlationMostDisplacementGrowthDivisionAngles = correlationMostDisplacementGrowthDivisionAngles;

  std::pair<double,double> correlationPrincipalGrowthDivisionAngles;
  dir1 = lAnalysis.principalGrowthVector;
  dir2 = lAnalysis.daughterDivisionDirections.first;
  dir1.normalize();
  dir2.normalize();
  correlationPrincipalGrowthDivisionAngles.first = std::acos( dir1*dir2 )*180./M_PI;
  dir2 = lAnalysis.daughterDivisionDirections.second;
  dir2.normalize();
  correlationPrincipalGrowthDivisionAngles.second = std::acos( dir1*dir2 )*180./M_PI;
  lAnalysis.correlationPrincipalGrowthDivisionAngles = correlationPrincipalGrowthDivisionAngles;

  // at last push the length analysis for this division cell
  _lengthAnalysisStats.push_back( lAnalysis );
}

// ---------------------------------------------------------------------

osg::Vec3 SCellAnalysis::determinePrincipalDirectionVector( const std::vector<osg::Vec3> &neighborDisplacementVectors )
{
  osg::Vec3 principalVector( 0., 0., 0. );

  for( std::size_t n = 0; n < neighborDisplacementVectors.size(); n++ )
    principalVector += neighborDisplacementVectors.at(n);

  // TODO: perhaps divide by size of neighborDisplacementVectors

  return principalVector;
}

//----------------------------------------------------------------------------

double SCellAnalysis::computeDivisionAngle( const SLineageTree *cell )
{
  // go past the division of the provided cell
  const SLineageTree *lastDiv = cell->parent;

  // find the last division
  while( lastDiv->parent && lastDiv->parent->children.size() < 2 )
    lastDiv = lastDiv->parent;

  // there is no prior division
  if( lastDiv == cell->parent || !lastDiv->parent ||
      lastDiv->parent->children.size() < 2 )
    return -1;

  // compute angle
  SArray dir1( lastDiv->getX() - lastDiv->parent->getX(),
               lastDiv->getY() - lastDiv->parent->getY(),
               lastDiv->getZ() - lastDiv->parent->getZ() );
  SArray dir2( cell->getX() - cell->parent->getX(),
               cell->getY() - cell->parent->getY(),
               cell->getZ() - cell->parent->getZ() );

  dir1.normalize();
  dir2.normalize();

  return std::acos( dir1 * dir2 )*180./M_PI;
}

// ---------------------------------------------------------------------

//void SCellAnalysis::drawDiffArrows( const osg::Vec3 &startPos,
//                                    const osg::Vec3 &endPos,
//                                    const osg::Vec3 &nStartPos,
//                                    const osg::Vec3 &nEndPos,
//                                    const osg::Vec3 &startMeanPos,
//                                    const osg::Vec3 &endMeanPos,
//                                    const std::size_t renderStep )
//{
//  // compute length of direction vector at the beginning
//  // and at the start
//  double startLength = (nStartPos - startPos).length();
//  double endLength;

//  // compute the length development for the first
//  // fixed cell position of the current cell
//  if( _firstPosFixed )
//    endLength = (nEndPos - startPos).length();
//  else
//    endLength = (nEndPos - endPos).length();

//  double deltaLength = fabs( endLength - startLength );

//  // has the length increased or decreased
//  int positive;
//  if( endLength > startLength )
//    positive = 1;
//  else
//    positive = 2;

//  osg::Vec3 end = endMeanPos - startMeanPos;
//  end.normalize();

//  end = end * deltaLength;
//  end = end + startMeanPos;

//  osg::Vec4 color;
//  if( positive == 1 )
//    color = _positiveLengthColor;
//  else
//    color = _negativeLengthColor;

//  osg::ref_ptr<osg::Group> renderGroup;

//  if( _useCycleSlider )
//    renderGroup = _analysisGroupPerCellCycle.at(renderStep)->getChild(positive)->asGroup();
//  else
//    renderGroup = _analysisGroupPerTimeStep.at(renderStep)->getChild(positive)->asGroup();

//  SSimilarityMeasureGraphics::addCylinderBetweenPoints( end, startMeanPos, 3.0,
//                                                    color, renderGroup );
//}

// ---------------------------------------------------------------------

void SCellAnalysis::drawArrows( const osg::Vec3 &start,
                                const osg::Vec3 &end,
                                const std::size_t renderTimeStep,
                                const std::size_t arrowType )
{
  if( arrowType > 3 )
    return;

  osg::Vec4 color;

  switch( arrowType )
  {
  case 0:
  default: color = _neighborColor; break;
  case 1: color = _divisionColor; break;
  case 2: color = _principalColor; break;
  case 3: color = _realPrincipalColor; break;
  }

  // draw arrow of neighbor length at the given render time step
  osg::ref_ptr<osg::Group> renderGroup;

  if( _useCycleSlider )
    renderGroup = _analysisGroupPerCellCycle.at(renderTimeStep)->getChild(arrowType)->asGroup();
  else
    renderGroup = _analysisGroupPerTimeStep.at(renderTimeStep)->getChild(arrowType)->asGroup();

  SSimilarityMeasureGraphics::addSimpleArrowBetweenPoints( end, start, color, color, renderGroup, _arrowRadius );
}

// ---------------------------------------------------------------------

//void SCellAnalysis::drawArrows( const std::vector<osg::Vec3> &nodePositions,
//                                const std::vector<osg::Vec3> &nNodePositions,
//                                const int startTimeStep,
//                                osg::Vec3 &startMeanPos,
//                                osg::Vec3 &endMeanPos,
//                                const std::size_t renderStep,
//                                const bool drawOnlyMean )
//{
//  osg::Vec3 m1( 0., 0., 0. );
//  osg::Vec3 m2( 0., 0., 0. );

//  int start = startTimeStep;

//  for( std::size_t i = 0; i < nodePositions.size(); i++, start++ )
//  {
//    m1 = m1 + nodePositions.at(i);
//    m2 = m2 + nNodePositions.at(i);

//    // of not the mean arrow should be drawn then render an arrow
//    // for each pair of node positions between the cell neighbors
//    if( !drawOnlyMean )
//    {
//      SSimilarityMeasureGraphics::addSimpleArrowBetweenPoints(
//            nNodePositions.at(i), nodePositions.at(i),
//            osg::Vec4( 0., 0., 1., 1. ), _neighborColor,
//            _analysisGroupPerTimeStep.at(start)->getChild(0)->asGroup(), 2.0 );
//    }
//  }

//  // decrease start parameter since it was increased by one in the
//  // last loop step
//  start--;

//  m1 = m1 / nodePositions.size();
//  m2 = m2 / nodePositions.size();

//  // if the first position should be fixed
//  // then only consider the first element
//  if( _firstPosFixed )
//    startMeanPos = nodePositions.front();
//  else
//    startMeanPos = m1;

//  endMeanPos = m2;

//  // draw mean arrow of neighbor length at the last time step
//  if( drawOnlyMean )
//  {
//    osg::ref_ptr<osg::Group> renderGroup;

//    if( _useCycleSlider )
//      renderGroup = _analysisGroupPerCellCycle.at(renderStep)->getChild(0)->asGroup();
//    else
//      renderGroup = _analysisGroupPerTimeStep.at(start)->getChild(0)->asGroup();

//    SSimilarityMeasureGraphics::addSimpleArrowBetweenPoints( endMeanPos, startMeanPos, _neighborColor, _neighborColor,
//                                                         renderGroup, 2.0 );
//  }
//}

// ---------------------------------------------------------------------

void SCellAnalysis::writeAnalysisResults( const std::string &fileName )
{
  // open file at the end of the stream and allow writing
  std::ofstream out( fileName.c_str(), std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to write into data "
                     + fileName );
    return;
  }

  out << "CellId TimeStep CellCycle NumberOfNeighbors "
      << "FirstAngleBetweenPrincipalGrowthAndDivision "
      << "SecondAngleBetweenPrincipalGrowthAndDivision "
      << "FirstAngleBetweenLongestDisplacementAndDivision "
      << "SecondAngleBetweenLongestDisplacementAndDivision "
      << "FirstConsecutiveDivisionAngle "
      << "SecondConsecutiveDivisionAngle\n";

  for( std::vector<SLengthAnalysis>::iterator iter = _lengthAnalysisStats.begin();
       iter != _lengthAnalysisStats.end(); ++iter )
  {
    out << iter->cellId << " "
        << iter->timeStep << " "
        << iter->cellCycle << " "
        << iter->numberOfNeighbors << " "
        << iter->correlationPrincipalGrowthDivisionAngles.first << " "
        << iter->correlationPrincipalGrowthDivisionAngles.second << " "
        << iter->correlationMostDisplacementGrowthDivisionAngles.first << " "
        << iter->correlationMostDisplacementGrowthDivisionAngles.second << " "
        << iter->daughterDivisionAngles.first << " "
        << iter->daughterDivisionAngles.second << "\n";
  }

  // at the end put an endline
  out << "\n";

  out.close();
}

//----------------------------------------------------------------------------

void SCellAnalysis::writeLengthParameters( const std::string &fileName )
{
  std::fstream out( fileName.c_str(), std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data "
                     + fileName );
    return;
  }

  out << "CellId TimeStep "
      << "FirstAngleBetweenPrincipalGrowthAndDivision "
      << "SecondAngleBetweenPrincipalGrowthAndDivision "
      << "FirstAngleBetweenLongestDisplacementAndDivision "
      << "SecondAngleBetweenLongestDisplacementAndDivision "
      << "MagnitudeOfPrincipalGrowthDirection "
      << "MagnitudeOfLongestDisplacementDirection\n";

  for( std::vector<SLengthAnalysis>::iterator iter = _lengthAnalysisStats.begin();
       iter != _lengthAnalysisStats.end(); ++iter )
  {
    out << iter->cellId << " "
        << iter->timeStep << " "
        << iter->correlationPrincipalGrowthDivisionAngles.first << " "
        << iter->correlationPrincipalGrowthDivisionAngles.second << " "
        << iter->correlationMostDisplacementGrowthDivisionAngles.first << " "
        << iter->correlationMostDisplacementGrowthDivisionAngles.second << " "
        << iter->principalGrowthVector.length() << " "
        << iter->mostDisplacementNeighborVector.length() << "\n";
  }

  out.close();
}

//----------------------------------------------------------------------------

void SCellAnalysis::writeSingleDisplacements( const std::string &fileName )
{
  std::fstream out( fileName.c_str(), std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data "
                     + fileName );
    return;
  }

  out << "CellId CellCycle TimeStep Angle Magnitude\n";
      //<< "SingleDisplacement PrincipalDisplacement\n";

  for( std::vector<SSingleNeighborDisplacement>::const_iterator iter = _allSingleDisplacements.begin();
       iter != _allSingleDisplacements.end(); ++iter )
  {
    for( std::map<std::size_t,osg::Vec3>::const_iterator mapIter = iter->singleDisplacements.begin();
         mapIter != iter->singleDisplacements.end(); mapIter++ )
    {
      osg::Vec3 a,b;
      a = iter->principalDisplacement;
      b = mapIter->second;
      a.normalize();
      b.normalize();

      double angle = std::acos( a * b );
      angle = angle * 180./M_PI;

      out << iter->cellId << " "
          << iter->cellCycle << " "
          << mapIter->first << " "
          << angle << " "
          << mapIter->second.length() << "\n";
    }
  }

  out.close();
}

//----------------------------------------------------------------------------

double SCellAnalysis::collinearConversion( const double angle )
{
  if( angle <= 90. )
    return angle;
  else
    return 180. - angle;
}

//----------------------------------------------------------------------------

double SCellAnalysis::computeAngleBetweenThisAndPriorDivision( const SLineageTree *cell )
{
  const SLineageTree *node = cell;

  // store current division direction by determining
  // the direction vector between both daughter cells
  osg::Vec3 currentDivDir = this->getDivisionDirectionFromDaughterCells( node );

  node = node->parent;

  while( node )
  {
    if( node->children.size() == 2 )
    {
      osg::Vec3  previousDivDir = this->getDivisionDirectionFromDaughterCells( node );
      double angle = acos( currentDivDir * previousDivDir );
      angle = angle * 180./M_PI;

      return angle;
    }
    node = node->parent;
  }

  return -1.;
}

//----------------------------------------------------------------------------

osg::Vec3 SCellAnalysis::getDivisionDirectionFromDaughterCells( const SLineageTree *cell )
{
  if( cell->children.size() == 2 )
  {
    osg::Vec3 child1( cell->children[0]->getX(), cell->children[0]->getY(), cell->children[0]->getZ() );
    osg::Vec3 child2( cell->children[1]->getX(), cell->children[1]->getY(), cell->children[1]->getZ() );
    osg::Vec3 divisionDir( child2 - child1 );
    divisionDir.normalize();

    return divisionDir;
  }
  else
  {
    std::cout << "No dividing cell.\n";
    return osg::Vec3( 0., 0., 0. );
  }
}

// ---------------------------------------------------------------------

