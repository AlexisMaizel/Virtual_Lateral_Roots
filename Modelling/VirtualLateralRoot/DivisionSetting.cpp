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
    _equalAreaRatio( equalAreaRatio ),
    _useCombinedAreaRatio( useCombinedAreaRatio ),
    _useWallRatio( useWallRatio ),
    _divisionArea( divisionArea ),
    _divisionAreaRatio( divisionAreaRatio ),
    _divisionWallRatio( divisionWallRatio ),
    _LODThreshold( LODThreshold ),
    _avoidTrianglesThreshold( avoidTrianglesThreshold ),
    _loadLastModel( loadLastModel ),
    _onlyGrowthInHeight( onlyGrowthInHeight ),
    _surfaceType( surfaceType ),
    _VLRBezierSurface( &VLRBezierSurface ),
    _VLRDataPointSurface( &VLRDataPointSurface ),
    _choiceCounter( 0 )
{
  
}

//----------------------------------------------------------------

void DivisionSetting::clear()
{
  _randomChoices.clear();
  _choiceCounter = 0;
}

//----------------------------------------------------------------

void DivisionSetting::setCellDivisionSettings( const cell &c,
                                               const std::size_t curTime,
                                               const bool fixedCenterPos )
{
  double area;
  Point3d center = ModelUtils::computeCellCenter( *_T, c, area, _accurateCenterOfMass );
  
  if( fixedCenterPos && c->cellCycle == 0 && _surfaceType == 0 )
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
  c->longestWallLength = ModelUtils::determineLongestWallLength( c, *_T ); 
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
   
    //_probValues.push_back( probValues );
    //_lengths.push_back( lengths );
    //_choices.push_back( choice );
    
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
      while( !ModelUtils::checkDivisionArea( c, ddata, *_T, _equalAreaRatio, areaDiff ) );
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
  
  // TODO: at the moment set constant ratio for radial model
  if( _surfaceType == 1 )
  {
    _divisionAreaRatio = 0.45;
    _avoidTrianglesThreshold = 15.;
    return;
  }
  
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