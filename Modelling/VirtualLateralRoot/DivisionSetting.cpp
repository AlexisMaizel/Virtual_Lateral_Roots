#include "DivisionSetting.h"

#include "ModelUtils.h"
#include "ModelExporter.h"

/**
  @file   DivisionSetting.cpp
  @brief  Class for defining division steps in model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

//----------------------------------------------------------------

DivisionSetting::DivisionSetting( MyTissue &T,
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
                     RealSurface &VLRDataPointSurface )
  : _T( &T ),
    _initialSituationType( initialSituationType ),
    _divisionType( divisionType ),
    _timeDelay( timeDelay ),
    _firstDivisionsAreaRatio( firstDivisionsAreaRatio ),
    _secondDivisionsAreaRatio( secondDivisionsAreaRatio ),
    _useAlternativeDT( useAlternativeDT ),
    _accurateCenterOfMass( accurateCenterOfMass ),
    _probabilityOfDecussationDivision( probabilityOfDecussationDivision ),
    _useAreaRatio( useAreaRatio ),
    _useCombinedAreaRatio( useCombinedAreaRatio ),
    _useWallRatio( useWallRatio ),
    _divisionArea( divisionArea ),
    _divisionAreaRatio( divisionAreaRatio ),
    _divisionWallRatio( divisionWallRatio ),
    _LODThreshold( LODThreshold ),
    _avoidTrianglesThreshold( avoidTrianglesThreshold ),
    _loadLastModel( loadLastModel ),
    _cellPinch( cellPinch ),
    _cellMaxPinch( cellMaxPinch ),
    _onlyGrowthInHeight( onlyGrowthInHeight ),
    _VLRBezierSurface( &VLRBezierSurface ),
    _VLRDataPointSurface( &VLRDataPointSurface )
{
  
}

//----------------------------------------------------------------
  
void DivisionSetting::step_divisions( const std::size_t curTime )
{
  // for which number of cells should the division area ratio check apply
    // this is required since at the beginning only few divisions occur due
    // to the small increasing of area based on the initial area size
    std::size_t areaRatioStart;
    if( _T->C.size() > 1 )
      areaRatioStart = 0;
    else
      areaRatioStart = 1;
    
    // Find cells to be divided
    std::list<cell> to_divide;
    
    // wait after the six cell stage until the predefined time has passed
    bool wait = false;
    if( _T->C.size() == 4 && _initialSituationType == 1 )
    {
      if( curTime - _timeFourCellStage > _timeDelay )
        wait = false;
      else
        wait = true;
    }
    
    if( _T->C.size() == 6 &&
        ( _initialSituationType == 2 || _initialSituationType == 3 ) )
    {
      if( curTime - _timeSixCellStage > _timeDelay )
        wait = false;
      else
        wait = true;
    }
    
    // at first determine the min and max values for each cell
    forall(const cell& c, _T->C)
      ModelUtils::determineXMinMax( c, *_T );
    
    // wait for the four-cell stage some time steps such that the
    // the future divisions are not occurring too fast
    // and only update the center of cells
    if( _T->C.size() == 4 && _initialSituationType == 1 && wait )
    {
      forall(const cell& c, _T->C)
        this->setCellDivisionSettings( c, curTime );
    }
    // force the initial start of the VLR
    else if( _T->C.size() < 6 && _initialSituationType != 0 )
    {
      forall(const cell& c, _T->C)
      {
        this->setCellDivisionSettings( c, curTime, true );
        
        double a = c->area;
        double l = c->longestWallLength;
        
        // divide cells if their area size has increased by a certain percentage amount
        double initialArea = c->initialArea;
        if( _T->C.size() < 4 )
          initialArea += initialArea*_firstDivisionsAreaRatio;
        else
          initialArea += initialArea*_secondDivisionsAreaRatio;
        
        if( _T->C.size() == 3 )
        {
          // when one of the two founder cells has already divided
          // then assure that the second founder cell will divide
          // before one of the two daughter cells of the first division
          // divide again
          if( c->id == 2 || c->id == 1 )
          {
            if( a > initialArea )
              to_divide.push_back(c);
          }
          
          _timeFourCellStage = curTime;
        }
        // else perform the "normal" division routine
        else
        {
          if( a > initialArea )
          {
            // if there are 5 cells and if we are in this loop
            // a next division will occur; so we store the time
            // of the six cell stage
            if( _T->C.size() == 5 )
              _timeSixCellStage = curTime;
              
            to_divide.push_back(c);
          }
        }
      }
    }
    // wait for the six-cell stage some time steps such that the
    // the future divisions are not occurring too fast
    // and only update the center of cells
    else if( _T->C.size() == 6 && _initialSituationType != 0 && wait )
    {
      forall(const cell& c, _T->C)
        this->setCellDivisionSettings( c, curTime );
    }
    else
    {      
      forall(const cell& c, _T->C)
      {
        this->setCellDivisionSettings( c, curTime );
        
        double a = c->area;
        double l = c->longestWallLength;
        if( _useAreaRatio && _useWallRatio && c->id > _areaRatioStart )
        {
          // divide cells if their area size has increased by a certain percentage amount
          double initialArea = c->initialArea;
          initialArea += initialArea*_divisionAreaRatio;
          // divide cell if its wall length has increased by a certain percentage amount
          double initialLongestLength = c->initialLongestWallLength;
          initialLongestLength += initialLongestLength*_divisionWallRatio;
          
          if( a > initialArea || l > initialLongestLength )
            to_divide.push_back(c);
        } 
        // apply a division if the cells area exceeds a certain threshold or area ratio
        else if( _useCombinedAreaRatio && c->id > _areaRatioStart )
        {
          // divide cells if their area size has increased by a certain percentage amount
          double initialArea = c->initialArea;
          initialArea += initialArea*_divisionAreaRatio;
          if( a > initialArea && a > _divisionArea )
            to_divide.push_back(c);
        }
        // only apply the division based on ratio with at least areaRatioStart cells
        else if( _useAreaRatio && c->id > _areaRatioStart )
        {
          // divide cells if their area size has increased by a certain percentage amount
          double initialArea = c->initialArea;
          initialArea += initialArea*_divisionAreaRatio;
          if( a > initialArea )
            to_divide.push_back(c);
        }
        // only apply the division based on ratio with at least areaRatioStart cells
        else if( _useWallRatio && c->id > _areaRatioStart )
        {
          // divide cell if its wall length has increased by a certain percentage amount
          double initialLongestLength = c->longestWallLength;
          initialLongestLength += initialLongestLength*_divisionWallRatio;
          if( l > initialLongestLength )
            to_divide.push_back(c);
        }
        else
        {
          if( a > _divisionArea )
            to_divide.push_back(c);
        }
      }
    }
    
    // Divide the cells
    forall(const cell& c, to_divide)
    {
      // perform the "normal" division routine for the initial two founder cells
      if( _T->C.size() < 4 && _initialSituationType != 0 )
        _T->divideCell(c);
      // set the division properties for the second division
      // periclinal division resulting in six cells
      // **********************************
      // *         *     *     *          *
      // *         *************          *
      // *         *     *     *          *
      // **********************************
      else if( _T->C.size() > 3 && _T->C.size() < 6 && _initialSituationType == 2 )
      {
        if( c->innerCell )
        {
          c->angle = 180.;
          MyTissue::division_data ddata = this->setDivisionPoints( c );
          _T->divideCell( c, ddata );
        }
      }
      // set the division properties for the second division
      // anticlinal division resulting in six cells
      // **********************************
      // *         *  *  *  *  *          *
      // *         *  *  *  *  *          *
      // *         *  *  *  *  *          *
      // **********************************
      else if( _T->C.size() > 3 && _T->C.size() < 6 && _initialSituationType == 3 )
      {
        if( c->innerCell )
        {
          c->angle = 90.;
          MyTissue::division_data ddata = this->setDivisionPoints( c );
          _T->divideCell( c, ddata );
        }
      }
      else if( _divisionType == "Decussation" && c->id > _areaRatioStart &&
               _useAlternativeDT )
      {
        // if true then the next division is perpendicular to the last one
        // else the division is collinear to the previous division direction
        if( ModelUtils::setNextDecussationDivision( _probabilityOfDecussationDivision ) )
        {          
          double noise = this->generateNoiseInRange( -10., 10., 1000 );
          c->angle = fmod( c->angle + 90., 360. ) + noise;
        }
        
        MyTissue::division_data ddata = this->setDivisionPoints( c );
        _T->divideCell( c, ddata );
      }
      else if( _divisionType == "PerToGrowth" && c->id > _areaRatioStart &&
               _useAlternativeDT )
      {
        // if true then the next division is perpendicular to the principal
        // component growth of the center deformation since the cell was born
        // determine principal growth direction
        if( c->centerPos.size() > 2 )
          c->principalGrowthDir = ModelUtils::determineLongestPCGrowth( c->centerPos );
        else
          std::cout << "Too few center positions!" << std::endl;
        
        Point3d xaxisDir = Point3d( 1., 0., 0. );
        double noise = this->generateNoiseInRange( -10., 10., 1000 );
        c->angle = 180./M_PI * acos( c->principalGrowthDir*xaxisDir );
        if( ModelUtils::determineSlope( c->principalGrowthDir ) > 0 )
        {
          // only add the property that it is perpendicular to the PCA dir
          c->angle += 90. + noise;
        }
        else
        {
          // add the property AND handle the missing degrees for having a negative slope
          if( c->principalGrowthDir.i() > 0 )
            c->angle = 180. - c->angle;
          
          c->angle += 90. + noise;
        }
        // update the division angle
        MyTissue::division_data ddata = this->setDivisionPoints( c );
        
        // store the line positions for later drawing
        Point3d midDiv = (ddata.pu+ddata.pv)/2.;
        _pcLines.push_back( std::make_pair( midDiv - 5.*c->principalGrowthDir,
                                            midDiv + 5.*c->principalGrowthDir ) );
        
        _T->divideCell( c, ddata );
      }
      else if( _divisionType == "Energy" && c->id > _areaRatioStart &&
               _useAlternativeDT )
      {
        // if true then all division planes are determined preserving
        // almost non-triangle cells with same area sizes and the shortest cell wall
        // has the lowest energy while the longest one has the highest energy;
        // probability values are then assigned like this:
        // lowest energy -> high probability
        // highest energy -> low probability
        bool empty = false;
        MyTissue::division_data ddata = this->getEnergyDivisionData( c, empty );
        // if the division data is empty then just use the default division
        // rule (e.g. ShortestWall)
        if( empty )
          _T->divideCell( c );
        else
          _T->divideCell( c, ddata );
      }
      else if( _divisionType == "Besson-Dumais" && c->id > _areaRatioStart &&
               _useAlternativeDT )
      {
        // the division plane is chosen based on Gibbs measure as described
        // in the paper of Besson and Dumais, 2011 (Universal rule for the 
        // symmetric division of plant cells)
        bool empty = false;
        MyTissue::division_data ddata = this->getBessonDumaisDivisionData( c, empty );
        // if the division data is empty then just use the default division
        // rule (e.g. ShortestWall)
        if( empty )
          _T->divideCell( c );
        else
          _T->divideCell( c, ddata );
      }
      else if( _divisionType == "RandomEqualAreas" && c->id > _areaRatioStart &&
               _useAlternativeDT )
      {
        // if true then all division planes are determined randomly preserving
        // cells with almost same area sizes
        bool empty = false;
        MyTissue::division_data ddata = this->getRandomDivisionData( c, empty, true );
        // if the division data is empty then just use the default division
        // rule (e.g. ShortestWall)
        if( empty )
          _T->divideCell( c );
        else
          _T->divideCell( c, ddata );
      }
      else if( _divisionType == "Random" && c->id > _areaRatioStart && _useAlternativeDT )
      {
        // if true then all division planes are determined randomly
        bool empty = false;
        MyTissue::division_data ddata = this->getRandomDivisionData( c, empty, false );
        // if the division data is empty then just use the default division
        // rule (e.g. ShortestWall)
        if( empty )
          _T->divideCell( c );
        else
          _T->divideCell( c, ddata );
      }
      else
        _T->divideCell(c);
    }

    return !to_divide.empty();  
}

//----------------------------------------------------------------

void DivisionSetting::setCellDivisionSettings( const cell &c,
                                               const std::size_t curTime,
                                               const bool fixedCenterPos )
{
  double area;
  Point3d center = ModelUtils::computeCellCenter( *_T, c, area, _accurateCenterOfMass );
  
  if( fixedCenterPos && _T->C.size() < 4 )
  {
    double xLength = std::fabs( c->xMax - c->xMin );
    if( c->id == 1 )
      center.i() = c->xMin + 2.*xLength/3.; 
    else if( c->id == 2 )
      center.i() = c->xMin + 1.*xLength/3.;
  }
  
  if( SURFACETYPE == 0 )
    _VLRBezierSurface->setPos( c->sp, center );
  else
    _VLRDataPointSurface->setPos( c->sp, center );
  c->area = area;
  c->center = center;
  c->centerPos.push_back( center );
  c->timeStep = curTime;
  c->longestWallLength = ModelUtils::determineLongestWallLength( c, &_T ); 
}

//----------------------------------------------------------------
  
MyTissue::division_data DivisionSetting::setDivisionPoints( const cell& c )
{
  MyTissue::division_data ddata;
  const Point3d& center = c->center;
  double a = M_PI/180. * c->angle;
  Point3d direction = Point3d(-sin(a), cos(a), 0);
  forall( const junction& j,_T->S.neighbors(c) )
  {
    Point3d jpos, jnpos;
    const junction& jn = _T->S.nextTo(c, j);
    jpos = j->getPos();
    jnpos = jn->getPos();
    
    Point3d u;
    double s;
    if(geometry::planeLineIntersection(u, s, center, direction, jpos, jnpos) and s >= 0 and s <= 1)
    {
      if((jpos - center)*direction > 0)
      {
        ddata.v1 = j;
        ddata.pv = u;
      }
      else if((jpos - center)*direction < 0)
      {
        ddata.u1 = j;
        ddata.pu = u;
      }
      if(ddata.u1 and ddata.v1)
        break;
    }
  }
  
  vvcomplex::testDivisionOnVertices(c, ddata, *_T, 0.01);
  
  // apply cell pinching
  tissue::CellPinchingParams params;
  params.cellPinch = _cellPinch;
  params.cellMaxPinch = _cellMaxPinch;
  tissue::cellPinching( c, &_T, ddata, params );
  
  return ddata;
}

//----------------------------------------------------------------
  
MyTissue::division_data DivisionSetting::getEnergyDivisionData( const cell& c,
                                                                bool &empty )
{
  std::vector<MyTissue::division_data> divData;
  divData = ModelUtils::determinePossibleDivisionData(
    c, _avoidTrianglesThreshold, _LODThreshold, *_T );
  
  // compute the lengths of all division lines and sort them
  // in ascending order
  std::vector<double> lengths;
  lengths.resize( divData.size() );
  double minLength = 5000.;
  double maxLength = 0.;
  double sumLength = 0.;
  std::size_t minIndex = 0;
  
  for( std::size_t l=0; l<lengths.size(); l++ )
  {
    Point3d pu = divData.at(l).pu;
    Point3d pv = divData.at(l).pv;
    //std::cout << "pu: " << pu << std::endl;
    //std::cout << "pv: " << pv << std::endl;
    double length = norm( pu - pv );
    lengths.at(l) = length;
    
    if( length < minLength )
    {
      minLength = length;
      minIndex = l;
    }
    
    if( length >= maxLength )
      maxLength = length;
    
    sumLength += length;
  }
  
  //std::cout << "possible Divisions: " << divData.size() << std::endl;
  //std::cout << "minLength: " << minLength << std::endl;
  //std::cout << "maxLength: " << maxLength << std::endl;
  //std::cout << "sumLength: " << sumLength << std::endl;
  
  std::size_t probDist = 1;
  
  MyTissue::division_data ddata;
  if( divData.size() != 0 )
  {
    std::vector<double> probValues;
    probValues.resize( divData.size() );
    
    if( probDist == 2 )
    {
      // get the probability distribution function for all lengths
      for( std::size_t l=0; l<probValues.size(); l++ )
      {
        double prob = (lengths.at(l)*100.)/sumLength;
        // the current probability value assigns a higher prob to longer lengths;
        // however, this should be inverse. Consequently, we substract the prob
        // from 100 percent and divide it by lengths.size()-1 to get the inverse
        // probability value
        if( lengths.size() > 1 )
          probValues.at(l) = (100. - prob);// /(lengths.size()-1);
        else
          probValues.at(l) = 100.;
        
        //std::cout << "prob value: " << probValues.at(l) << std::endl;
      }
    }
    else
    {
      double mu = 0.;
      double sd = std::sqrt(0.5);
      std::vector<double> normLengths;
      normLengths.resize( divData.size() );
      for( std::size_t l=0; l<normLengths.size(); l++ )
      {
        normLengths.at(l) = (lengths.at(l)-minLength)/(maxLength-minLength);
        double prob = ModelUtils::getNormalDistribution( normLengths.at(l), mu, sd );
    
        if( lengths.size() > 1 )
          probValues.at(l) = prob*100.;
        else
          probValues.at(l) = 100.;
        
        //std::cout << l << " prob value: " << probValues.at(l) << std::endl;
      }
    }
    
    std::size_t choice = ModelUtils::getRandomResultOfDistribution( probValues );
    
    // store choice for later usage
    if( !_loadLastModel )
      _randomChoices.push_back( choice );
    else
    {
      if( _choiceCounter < _randomChoices.size() )
      {
        choice = _randomChoices.at(_choiceCounter);
        ddata = divData.at( choice );
        _choiceCounter++;
      }
    }
    
    //std::cout << "choice: " << choice << std::endl;
    ddata = divData.at( choice );
    vvcomplex::testDivisionOnVertices(c, ddata, *_T, 0.01);
    
    // apply cell pinching
    tissue::CellPinchingParams params;
    params.cellPinch = _cellPinch;
    params.cellMaxPinch = _cellMaxPinch;
    tissue::cellPinching( c, *_T, ddata, params );
    empty = false;
  }
  else
    empty = true;
  
  return ddata;
}

//----------------------------------------------------------------

MyTissue::division_data DivisionSetting::getBessonDumaisDivisionData( const cell& c, bool &empty )
{
  std::vector<MyTissue::division_data> divData;
  divData = ModelUtils::determinePossibleDivisionData(
    c, _avoidTrianglesThreshold, _LODThreshold, *_T );
  
  std::size_t numLines = divData.size();
  //std::cout << "possible Divisions: " << numLines << std::endl;
  
  // fixed parameter from paper determined by experimental results
  const double beta = 20.6;
  
  // compute the lengths of all division lines and sort them
  // in ascending order
  std::vector<double> lengths;
  lengths.resize( numLines );
  
  // compute area size and mean cell diameter
  std::vector<Point3d> polygon;
  forall(const junction& j, _T->S.neighbors(c))
    polygon.push_back( j->getPos() );
  
  double area = geometry::polygonArea( polygon );
  // mean cell diameter
  double rho = std::sqrt( area );
  
  for( std::size_t l=0; l<numLines; l++ )
  {
    Point3d pu = divData.at(l).pu;
    Point3d pv = divData.at(l).pv;
    //std::cout << "pu: " << pu << std::endl;
    //std::cout << "pv: " << pv << std::endl;
    lengths.at(l) = norm( pu - pv );
    //std::cout << l << " length: " << lengths.at(l) << std::endl;
  }
  
  MyTissue::division_data ddata;
  if( numLines != 0 )
  {
    double sum = 0.;
    for( std::size_t l=0; l<numLines; l++ )
      sum += std::exp( (-1. * beta * lengths.at(l)) / rho );
    
    std::vector<double> probValues;
    probValues.resize( numLines );
    for( std::size_t l=0; l<numLines; l++ )
    {
      probValues.at(l) = 100. * std::exp( (-1. * beta * lengths.at(l)) / rho ) / sum;
      //std::cout << l << " prob value: " << probValues.at(l) << std::endl;
      //ddata = divData.at( l );
      //ModelUtils::checkDivisionArea( c, ddata, T, _equalAreaRatio );
    }
    
    //unsigned int attempts = 0;
    std::size_t choice = ModelUtils::getRandomResultOfDistribution( probValues );
    //std::cout << "choice: " << choice << std::endl;
    //std::cout << " selected prob value: " << probValues.at(choice) << std::endl;
    //ddata = divData.at( choice );
    //ModelUtils::checkDivisionArea( c, ddata, T, _equalAreaRatio );
    
//       if( true )
//       {
//         // check the two possible area of the daughter cells
//         do
//         {
//           choice = ModelUtils::getRandomResultOfDistribution( probValues );
//   //         std::cout << "c: " << choice << std::endl;
//           ddata = divData.at( choice );
//           attempts++;
//           
//           if( attempts > 50 )
//             break;
//         }
//         while( !ModelUtils::checkDivisionArea( c, ddata, T, _equalAreaRatio ) );
//       }
//       else
//       {
//         choice = ModelUtils::getRandomResultOfDistribution( probValues );
//         ddata = divData.at( choice );
//       }
    
    // store choice for later usage
    if( !_loadLastModel )
      _randomChoices.push_back( choice );
    else
    {
      if( _choiceCounter < _randomChoices.size() )
      {
        choice = _randomChoices.at(_choiceCounter);
        _choiceCounter++;
      }
    }
    
    ddata = divData.at( choice );
    
    vvcomplex::testDivisionOnVertices(c, ddata, *_T, 0.01);
    
    _probValues.push_back( probValues );
    _lengths.push_back( lengths );
    _choices.push_back( choice );
    
    // apply cell pinching
    tissue::CellPinchingParams params;
    params.cellPinch = _cellPinch;
    params.cellMaxPinch = _cellMaxPinch;
    tissue::cellPinching( c, *_T, ddata, params );
    empty = false;
  }
  else
    empty = true;
  
  return ddata;
}

//----------------------------------------------------------------
  
MyTissue::division_data DivisionSetting::getRandomDivisionData( const cell& c,
                                                                bool &empty,
                                                                const bool equalArea )
{
  std::vector<MyTissue::division_data> divData;
  divData = ModelUtils::determinePossibleDivisionData(
    c, _avoidTrianglesThreshold, _LODThreshold, *_T );
  
//     std::cout << "possible Divisions: " << divData.size() << std::endl;
  
  MyTissue::division_data ddata;
  unsigned int attempts = 0;
  if( divData.size() != 0 )
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> d( 0, divData.size()-1 );
    
    // get a random choice of all possible division data
    std::size_t choice;
    double areaDiff;
    
    if( equalArea )
    {
      // testing
//         choice = ModelUtils::determineDivisionDataBasedOnAlmostEqualArea( c, divData, T, 50 );
//         ddata = divData.at( choice );
      
      // check the two possible area of the daughter cells
      do
      {
        choice = d(gen);
        ddata = divData.at( choice );
        attempts++;
        
        if( attempts > 50 )
          break;
      }
      while( !ModelUtils::checkDivisionArea( c, ddata, T, _equalAreaRatio, areaDiff ) );
    }
    else
    {
      choice = d(gen);
      ddata = divData.at( choice );
    }
    
    // store choice for later usage
    if( !_loadLastModel )
      _randomChoices.push_back( choice );
    else
    {
      if( _choiceCounter < _randomChoices.size() )
      {
        choice = _randomChoices.at(_choiceCounter);
        ddata = divData.at( choice );
        _choiceCounter++;
      }
    }
    
    vvcomplex::testDivisionOnVertices(c, ddata, *_T, 0.01);
    
    // apply cell pinching
    tissue::CellPinchingParams params;
    params.cellPinch = _cellPinch;
    params.cellMaxPinch = _cellMaxPinch;
    tissue::cellPinching( c, *_T, ddata, params );
    empty = false;
  }
  else
    empty = true;
  
  return ddata;
}

//----------------------------------------------------------------

double DivisionSetting::generateNoiseInRange( const double start,
                                              const double end,
                                              const unsigned int steps )
{
  std::vector<double> values;
  double stepSize = fabs( end - start )/steps;
  
  for( double v = start; v < end+stepSize; v+=stepSize )
    values.push_back( v );
  
  // store choice for later usage
  std::size_t choice;
  if( !_loadLastModel )
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> d( 0, values.size()-1 );
    choice = d(gen);
    _randomChoices.push_back( choice );
  }
  else
  {
    if( _choiceCounter < _randomChoices.size() )
    {
      choice = _randomChoices.at(_choiceCounter);
      _choiceCounter++;
    }
  }
  
  return values.at( choice );
}

//----------------------------------------------------------------

void DivisionSetting::step_cellWalls( const std::string &fileName )
{
  forall(const cell& c, _T->C)
    ModelExporter::exportCellWalls( fileName, c, *_T, false );
}

//----------------------------------------------------------------

void DivisionSetting::step_tracking( const std::string &fileName )
{
  forall(const cell& c, _T->C)
    ModelExporter::exportLineageInformation( fileName, c, *_T, false );
}

//----------------------------------------------------------------

void DivisionSetting::step_growth( const double dt )
{
  if( SURFACETYPE == 0 )
  {
    _VLRBezierSurface->growStep( dt );
    forall(const junction& j, _T->W)
      _VLRBezierSurface->getPos( j->sp );
  }
  else
  {
    // generate new triangulation surface based on real data points
    std::vector<SurfacePoint> sps;
    forall(const junction& j, _T->W)
      sps.push_back( j->sp );

    _VLRDataPointSurface->growStep( dt, sps );
    
    int i = 0;
    forall(const junction& j, _T->W)
    {
      j->sp = sps[i++];
      _VLRDataPointSurface->getPos( j->sp );
    }
  }
}

//----------------------------------------------------------------
 
void DivisionSetting::setAreaRatios( const bool loadLastModel,
                                     const std::size_t initialSituationType )
{
  _loadLastModel = loadLastModel;
  _initialSituationType = initialSituationType;
  
  switch( _initialSituationType )
  {
    case 0:
      if( _divisionType == "Besson-Dumais" )
      {
        if( _onlyGrowthInHeight )
        {
          // 90 cells at end: 0.25
          // 64 cells at end: 0.32
          _divisionAreaRatio = 0.3;
        }
        else
        {
          // asymmetric case: 0.455
          // growth in height/width: 0.468
          _divisionAreaRatio = 0.468; 
        }
      }
      else if( _divisionType == "Random" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.308;
        else
          _divisionAreaRatio = 0.4674;
      }
      else if( _divisionType == "RandomEqualAreas" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.308;
        else
          _divisionAreaRatio = 0.4674;
      }
      else if( _divisionType == "Decussation" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.3;
        else
          _divisionAreaRatio = 0.462;
      }
      else if( _divisionType == "PerToGrowth" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.3;
        else
          _divisionAreaRatio = 0.462;
      }
    break;
    case 1:
      if( _divisionType == "Besson-Dumais" )
      {
        if( _onlyGrowthInHeight )
        {
          // 90 cells at end: 0.305
          // 64 cells at end: 0.42
          _divisionAreaRatio = 0.42;
        }
        else
        {
          // asymmetric case: ?
          // growth in height/width: 0.7
          _divisionAreaRatio = 0.7;
        }
      }
      else if( _divisionType == "Random" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.45;
        else
          _divisionAreaRatio = 0.4674;
      }
      else if( _divisionType == "RandomEqualAreas" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.45;
        else
          _divisionAreaRatio = 0.4674;
      }
      else if( _divisionType == "Decussation" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.445;
        else
          _divisionAreaRatio = 0.462;
      }
      else if( _divisionType == "PerToGrowth" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.445;
        else
          _divisionAreaRatio = 0.462;
      }
    break;
    case 2:
    case 3:
      if( _divisionType == "Besson-Dumais" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.45;
        else
        {
          // asymmetric case: 0.6
          // growth in height/width: 0.71
          _divisionAreaRatio = 0.71;
        }
      }
      else if( _divisionType == "Random" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.45;
        else
          _divisionAreaRatio = 0.735;
      }
      else if( _divisionType == "RandomEqualAreas" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.45;
        else
          _divisionAreaRatio = 0.735;
      }
      else if( _divisionType == "Decussation" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.45;
        else
          _divisionAreaRatio = 0.73;
      }
      else if( _divisionType == "PerToGrowth" )
      {
        if( _onlyGrowthInHeight )
          _divisionAreaRatio = 0.45;
        else
          _divisionAreaRatio = 0.71;
      }
    break;
  }
  
  if( _divisionType == "Random" || _divisionType == "RandomEqualAreas" )
    _avoidTrianglesThreshold = 15.;
  else
    _avoidTrianglesThreshold = 0.;
}

//----------------------------------------------------------------
  
}