#include "SDivisionAnalysis.hh"

/**
  @file   SDivisionAnalysis.cc
  @brief  Class for performing a division analysis
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SInteractiveCell.hh"
#include "SInteractiveDivisionArrows.hh"
#include "SInteractiveNormals.hh"

#include "kernel/SKernel.hh"

#include "common/VException.hh"

#include <string>
#include <fstream>

#include <boost/foreach.hpp>

//----------------------------------------------------------------------------

SDivisionAnalysis::SDivisionAnalysis( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                                      const std::size_t numTotalTimeSteps,
                                      const std::size_t numCellCycles,
                                      const std::size_t numRegSteps,
                                      const bool onlyMasterCellFile,
                                      osg::ref_ptr<SSliderCallback> sliderCallback,
                                      std::vector<osg::ref_ptr<osg::MatrixTransform> > addToThisGroup,
                                      const std::size_t sliderType,
                                      const int angleType )
  : _lineages( lineages ),
    _numTotalTimeSteps( numTotalTimeSteps ),
    _numCellCycles( numCellCycles ),
    _numRegSteps( numRegSteps ),
    _onlyMasterCellFile( onlyMasterCellFile ),
    _rotMatrix( addToThisGroup.front()->getMatrix() ),
    _inverseRotMatrix( addToThisGroup.front()->getInverseMatrix() ),
    _sliderType( sliderType ),
    _angleType( angleType ),
    _currentView( 0 )
{
  osg::ref_ptr<osg::Group> cellSurfaceGroup = new osg::Group;
  _renderGroup = new osg::Group;

  // first determine the number and positions of cells for each time step
  std::map<int,int> cellFiles = _lineages->getCellFiles();
  _cellsPerTimeStep.resize( _numTotalTimeSteps );
  _centresOfMass.resize( _numTotalTimeSteps );
  for( auto l = _lineages->begin(); l != _lineages->end(); ++l )
  {
    if( _onlyMasterCellFile && cellFiles.at(l->first) != 0 )
      continue;

    SLineageTree *tree = l->second;
    for( auto nodeIt = tree->begin(); nodeIt != tree->end(); ++nodeIt )
    {
      _cellsPerTimeStep.at(nodeIt->timeStep-1).insert( *nodeIt );
      osg::Vec3 pos( nodeIt->getX(), nodeIt->getY(), nodeIt->getZ() );
      _centresOfMass.at(nodeIt->timeStep-1) += pos;
    }
  }

  for( std::size_t t=0; t<_numTotalTimeSteps; t++ )
  {
    if( _cellsPerTimeStep.at(t).size() > 0 )
      _centresOfMass.at(t) /= _cellsPerTimeStep.at(t).size();
  }

  this->performAngleAnalysis();
  this->readAngleData();

  this->performCylindricalAngleAnalysis();
  _surfaceGroup.resize( 3 );
  _minMaxBoundary.resize( 3 );
  this->initRegDivisionProp();
  this->registerBasedOnNumberOfCells();
  this->determineDivisionsPerCellCycle();

  for( std::size_t v=0; v < 3; v++ )
  {
    _minMaxBoundary.at(v).resize( 3, std::make_pair( 50000., -50000. ) );

    _currentView = v;

    if( _sliderType == 0 || _sliderType == 2 )
    {
      SInteractiveNormals *is = new SInteractiveNormals();
      _renderGroup->setUpdateCallback( sliderCallback );
      is->setName( "Divisions" );
      is->addUpdateCallback( sliderCallback );
      is->addChild( cellSurfaceGroup );

      std::size_t numTotalSteps;
      if( _sliderType == 0 )
        numTotalSteps = _numTotalTimeSteps;
      else
        numTotalSteps = _numRegSteps;

      for( std::size_t t=0; t<numTotalSteps; t++ )
      {
        osg::ref_ptr<osg::Group> group = new osg::Group;
        _surfaceGroup.at(v).push_back( group );
      }

      this->renderSurfaceContours();
      //this->renderDivisions( _angleType );

      for( std::size_t t=0; t<_surfaceGroup.at(v).size(); t++ )
        is->addNormals( t+1, _surfaceGroup.at(v).at(t) );

      // for the time slider add 50 empty time steps
      // because only one data set has 350 time steps while the
      // remaining ones only have 300 steps
      if( _sliderType == 0 && _surfaceGroup.at(v).size() == 300 )
      {
        for( std::size_t i = 0; i < 50; i++ )
          is->addNormals( 300 + i, new osg::Group );
      }

      addToThisGroup.at(v)->addChild( is );
    }
    else if( _sliderType == 1 )
    {
      SInteractiveDivisionArrows *is = new SInteractiveDivisionArrows();
      _renderGroup->setUpdateCallback( sliderCallback );
      is->setName( "Divisions" );
      is->addUpdateCallback( sliderCallback );
      is->addChild( cellSurfaceGroup );

      for( std::size_t c=0; c<_numCellCycles; c++ )
      {
        osg::ref_ptr<osg::Group> group = new osg::Group;
        _surfaceGroup.at(v).push_back( group );
      }

      this->renderSurfaceContours();
      //this->renderDivisions( _angleType );

      for( std::size_t t=0; t<_surfaceGroup.at(v).size(); t++ )
        is->addDivisionArrows( t+1, _surfaceGroup.at(v).at(t) );

      // for the cell cycle slider add an empty sixth cycle
      // because only one data set has 5 cell cycles
      if( _sliderType == 1 && _surfaceGroup.at(v).size() == 5 )
        is->addDivisionArrows( 6, new osg::Group );

      addToThisGroup.at(v)->addChild( is );
    }
  }
}

//----------------------------------------------------------------------------

std::size_t SDivisionAnalysis::getRegIndex( const std::size_t numCells )
{
  std::size_t startC = 18;
  std::size_t endC = 143;
  double stepSize = double(endC - startC)/double(_numRegSteps);
  endC = static_cast<std::size_t>( (double)endC - stepSize );

  if( numCells < 18 || numCells > 143 )
  {
    std::cout << "Number of cells is out of range!" << std::endl;
    return 0;
  }

  double step = startC + stepSize;
  std::size_t index = 0;
  while( numCells > (std::size_t)step )
  {
    step += stepSize;
    index++;
  }

  return index;
}

//----------------------------------------------------------------------------

void SDivisionAnalysis::initRegDivisionProp()
{
  std::size_t startC = 18;
  std::size_t endC = 143;
  double stepSize = double(endC - startC)/double(_numRegSteps);
  endC = static_cast<std::size_t>( (double)endC - stepSize );
  double step = startC;

  for( std::size_t s = 0; s < _numRegSteps; s++ )
  {
    std::size_t iStep = (std::size_t)step;
    std::set<Angles, setComp> angleSet;
    _registeredDivisionProperties.insert( std::make_pair( iStep, angleSet ) );

    step += stepSize;
    if( step > endC )
      step = endC;
  }
}

//----------------------------------------------------------------------------

void SDivisionAnalysis::registerBasedOnNumberOfCells()
{
  std::size_t startC = 18;
  std::size_t endC = 143;
  double stepSize = double(endC - startC)/double(_numRegSteps);
  endC = static_cast<std::size_t>( (double)endC - stepSize );
  for( std::size_t t = 0; t < _divisionProperties.size(); t++ )
  {
    for( auto setIter = _divisionProperties.at(t).begin();
         setIter != _divisionProperties.at(t).end(); ++setIter )
    {
      std::size_t numCells = setIter->numberOfCells;
      double step = startC + stepSize;

      // ignore these divisions with number of cells out of the ranges
      if( numCells < 18 || numCells > 143 )
        continue;

      while( numCells > (std::size_t)step )
        step += stepSize;

      std::size_t iStep = (std::size_t)(step - stepSize);

      if( _registeredDivisionProperties.count( iStep ) != 0 )
        _registeredDivisionProperties.at( iStep ).insert( *setIter );
    }
  }

//  std::cout << "size: " << _registeredDivisionProperties.size() << std::endl;

  for( auto mapIter = _registeredDivisionProperties.begin();
       mapIter != _registeredDivisionProperties.end(); ++mapIter )
  {
    //std::cout << "regsteps: " << mapIter->first << std::endl;
    int avgTimeStep = 0;
    for( auto setIter = mapIter->second.begin(); setIter != mapIter->second.end(); ++setIter )
      avgTimeStep += setIter->timeStep;

    if( mapIter->second.size() != 0 )
      avgTimeStep = static_cast<std::size_t>( (double)avgTimeStep/(double)mapIter->second.size() );
    else
      avgTimeStep = -1;

    _registeredShapes.push_back( avgTimeStep );
//    std::cout << "avgTimeteps: " << avgTimeStep << std::endl;
  }
}

//----------------------------------------------------------------------------

std::vector< std::set<Angles, setComp> > SDivisionAnalysis::getRegisteredDivisionProperties()
{
  std::vector< std::set<Angles, setComp> > divProp;

  for( auto mapIter = _registeredDivisionProperties.begin();
       mapIter != _registeredDivisionProperties.end(); ++mapIter )
  {
    divProp.push_back( mapIter->second );
//    std::cout << "range: " << mapIter->first << std::endl;
//    for( auto setIter = mapIter->second.begin(); setIter != mapIter->second.end(); ++setIter )
//      std::cout << "numCells: " << setIter->numberOfCells << std::endl;
  }

  return divProp;
}

//----------------------------------------------------------------------------

void SDivisionAnalysis::determineDivisionsPerCellCycle()
{
  std::vector<int> divisionPerCellCycle;
  _averagedCycleShapeTimes.resize( _numCellCycles, 0 );
  divisionPerCellCycle.resize( _numCellCycles, 0 );

  for( SAdvancedLineageTreeCollection::const_iterator l = _lineages->begin();
       l != _lineages->end(); ++l )
  {
    for( SLineageTree::const_iterator iter = l->second->begin();
         iter != l->second->end(); ++iter )
    {
      if( iter->children.size() == 2 )
      {
        std::size_t cellCycle = iter->getCellCycleId();
        _averagedCycleShapeTimes.at(cellCycle) += iter->timeStep-1;
        divisionPerCellCycle.at(cellCycle)++;
      }
    }
  }

  for( std::size_t i=0; i < _numCellCycles; i++ )
  {
    if( divisionPerCellCycle.at(i) > 0 )
      _averagedCycleShapeTimes.at(i) /= divisionPerCellCycle.at(i);
    else
      _averagedCycleShapeTimes.at(i) = -1;
  }
}

//----------------------------------------------------------------------------

void SDivisionAnalysis::renderSurfaceContours()
{
  std::size_t numTotalSteps;
  if( _sliderType == 0 )
    numTotalSteps = _numTotalTimeSteps;
  else if( _sliderType == 1 )
    numTotalSteps = _numCellCycles;
  else
    numTotalSteps = _numRegSteps;

  // generate the alpha shapes for the averaged cell cycles
  for( std::size_t s=0; s < numTotalSteps; s++ )
  {
    SKernel::get()->setProgress( (double)s/(double)(numTotalSteps - 1),
                                 "Contour generation" );

    int t;
    if( _sliderType == 0 )
      t = s;
    else if( _sliderType == 1 )
      t = _averagedCycleShapeTimes.at(s);
    else
      t = _registeredShapes.at(s);

    if( t != -1 )
    {
      std::vector<osg::Vec3> boundary = _alphaShapes.at(t)->getOuterBoundary( _currentView, _rotMatrix );
      this->renderContour( boundary, s );
      this->updateBoundaryValues( boundary );
    }
  }
}

//----------------------------------------------------------------------------

void SDivisionAnalysis::updateBoundaryValues( const std::vector<osg::Vec3> &boundary )
{
  for( std::size_t p = 0; p < boundary.size(); p++ )
  {
    double value = boundary.at(p).x();
    this->checkMinMaxBoundary( value, 0 );
    value = boundary.at(p).y();
    this->checkMinMaxBoundary( value, 1 );
    value = boundary.at(p).z();
    this->checkMinMaxBoundary( value, 2 );
  }
}

//----------------------------------------------------------------------------

void SDivisionAnalysis::checkMinMaxBoundary( const double value,
                                             const std::size_t dim )
{
  if( value > _minMaxBoundary.at(_currentView).at(dim).second )
    _minMaxBoundary.at(_currentView).at(dim).second = value;

  if( value <= _minMaxBoundary.at(_currentView).at(dim).first )
    _minMaxBoundary.at(_currentView).at(dim).first = value;
}

//----------------------------------------------------------------------------

void SDivisionAnalysis::renderContour( const std::vector<osg::Vec3> &boundary,
                                       const std::size_t step )
{
  osg::Vec4 color( 0., 0., 0., 1. );
  for( std::size_t b = 0; b < boundary.size(); b+=2 )
  {
    osg::Vec3 p1, p2;
    p1 = boundary.at(b);
    p1 = p1 * _inverseRotMatrix;
    p2 = boundary.at(b+1);
    p2 = p2 * _inverseRotMatrix;
//    SSimilarityMeasureGraphics::addCylinderBetweenPoints( p1, p2, 1., color,
//                                                      _surfaceGroup.at(_currentView).at(step) );
    osg::ref_ptr<osg::Geode> geode = new osg::Geode;

    SSimilarityMeasureGraphics::drawLine( geode, p1, p2, 2., color );
    _surfaceGroup.at(_currentView).at(step)->addChild( geode );
  }
}

//----------------------------------------------------------------------------

//void SDivisionAnalysis::renderDivisions( const int angleType )
//{
//  for( SAdvancedLineageTreeCollection::const_iterator l = _lineages->begin();
//       l != _lineages->end(); ++l )
//  {
//    for( SLineageTree::const_iterator iter = l->second->begin();
//         iter != l->second->end(); ++iter )
//    {
//      if( iter->children.size() == 2 )
//      {
//        osg::Vec3 c1( iter->children[0]->getX(),
//            iter->children[0]->getY(),
//            iter->children[0]->getZ() );

//        osg::Vec3 c2( iter->children[1]->getX(),
//            iter->children[1]->getY(),
//            iter->children[1]->getZ() );

//        double ang;
//        switch( angleType )
//        {
//          case 0: ang = _angleValues.at( std::make_pair( iter->cellId, iter->timeStep ) ).VLRangle; break;
//          case 1: ang = _angleValues.at( std::make_pair( iter->cellId, iter->timeStep ) ).primAngle; break;
//          case 2: ang = _angleValues.at( std::make_pair( iter->cellId, iter->timeStep ) ).surfaceAngle; break;
//        }

//        osg::Vec4 color = _colorMaps.at(angleType)->mapToColor( ang );
//        if( _sliderType == 0 )
//          SSimilarityMeasureGraphics::addCylinderBetweenPoints( c1, c2, 2, color,
//                                                            _surfaceGroup.at(_currentView).at(iter->timeStep-1) );
//        else if( _sliderType == 1 )
//        {
//          std::size_t cellCycle = iter->getCellCycleId();
//          SSimilarityMeasureGraphics::addCylinderBetweenPoints( c1, c2, 2., color,
//                                                            _surfaceGroup.at(_currentView).at(cellCycle) );
//        }
//        else if( _sliderType == 2 )
//        {
//          std::size_t numCells = _angleValues.at( std::make_pair( iter->cellId, iter->timeStep ) ).numberOfCells;
//          if( numCells >= 18 && numCells <= 143 )
//          {
//            std::size_t regIndex = this->getRegIndex( numCells );
//            SSimilarityMeasureGraphics::addCylinderBetweenPoints( c1, c2, 2., color,
//                                                              _surfaceGroup.at(_currentView).at(regIndex) );
//          }
//        }
//      }
//    }
//  }
//}

//----------------------------------------------------------------------------

void SDivisionAnalysis::performAngleAnalysis()
{
  SArray ce = _lineages->getCenter();
  osg::Vec3 primCen( ce[0], ce[1], ce[2] );
  primCen = primCen * _rotMatrix;

  _minMaxAngleValues.resize( 6, std::make_pair( 360., 0. ) );

  std::string name = _lineages->getName();

  // generate the alpha shapes for each time step
  for( std::size_t t=0; t < _numTotalTimeSteps; t++ )
  {
    SKernel::get()->setProgress( (double)t/(double)(_numTotalTimeSteps - 1),
                                 name + ": AS generation" );

    osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array;
    BOOST_FOREACH( const SLineageTree *tree, _cellsPerTimeStep.at(t) )
    {
      osg::Vec3 p( tree->getX(), tree->getY(), tree->getZ() );
      points->push_back( p );
    }

    boost::shared_ptr<SAlphaShape> tri = boost::shared_ptr<SAlphaShape>(
          new SAlphaShape( points, t, osg::Vec4(1., 0., 1., 1.), false, false, 0. ) );
    _alphaShapes.push_back( tri );
  }

  bool firstDivTime = true;

  std::string fileName = "/tmp/DivisionAnalysis_" + _lineages->getName() + ".csv";
  // open file at the end of the stream and allow writing
  std::ofstream out( fileName.c_str(), std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to write into data " + fileName );
    return;
  }

  // header
  out << "CellID TimeStep CellCycle DivPosX DivPosY DivPosZ";
  out << " DauPos1X DauPos1Y DauPos1Z DauPos2X DauPos2Y DauPos2Z";
  out << " CentreOfMassX CentreOfMassY CentreOfMassZ";
  out << " VLRAngle DistFromDivToCentreOfMass PrimAngle DistFromDivToCentreOfMainRoot";
  out << " SurfaceAngle DistFromDivToSurface\n";

  _divisionProperties.resize( _numTotalTimeSteps );
  std::size_t totalLineages = _lineages->size();
  std::size_t li = 0;
  for( SAdvancedLineageTreeCollection::const_iterator l = _lineages->begin();
       l != _lineages->end(); ++l, li++ )
  {
    SKernel::get()->setProgress( (double)li/(double)(totalLineages - 1),
                                 name + ": Angle computations" );

    for( SLineageTree::const_iterator iter = l->second->begin();
         iter != l->second->end(); ++iter )
    {
      if( iter->children.size() == 2 )
      {
        osg::Vec3 divPos( iter->getX(), iter->getY(), iter->getZ() );
        divPos = divPos * _rotMatrix;

        osg::Vec3 c1, c2;
        if( !firstDivTime )
        {
          // first child
          SLineageTree::const_iterator iter1 = iter->children[0];
          while( iter1->children.size() != 2 )
          {
            if( iter1->children.size() == 0 )
              break;
            else
              ++iter1;
          }

          c1 = osg::Vec3( iter1->getX(), iter1->getY(), iter1->getZ() );

          // second child
          SLineageTree::const_iterator iter2 = iter->children[1];
          while( iter2->children.size() != 2 )
          {
            if( iter2->children.size() == 0 )
              break;
            else
              ++iter2;
          }

          c2 = osg::Vec3( iter2->getX(), iter2->getY(), iter2->getZ() );
        }
        else
        {
          c1 = osg::Vec3( iter->children[0]->getX(),
              iter->children[0]->getY(),
              iter->children[0]->getZ() );

          c2 = osg::Vec3( iter->children[1]->getX(),
              iter->children[1]->getY(),
              iter->children[1]->getZ() );
        }

        c1 = c1 * _rotMatrix;
        c2 = c2 * _rotMatrix;

        osg::Vec3 cen = _centresOfMass.at(iter->timeStep-1);
        cen = cen * _rotMatrix;
        osg::Vec3 dirCenter = divPos - cen;
        osg::Vec3 dirPrimCenter = divPos - primCen;
        osg::Vec3 dir;
        dirCenter.normalize();
        dirPrimCenter.normalize();

        // always consider the angle for which one of the daughter
        // cells is nearer to the center
        if( (cen - c1).length() < (cen - c2).length() )
          dir = c2 - c1;
        else
          dir = c1 - c2;

        dir.normalize();
        double angle1 = acos( dir * dirCenter )*180./M_PI;

        if( (primCen - c1).length() < (primCen - c2).length() )
          dir = c2 - c1;
        else
          dir = c1 - c2;

        dir.normalize();
        double angle2 = acos( dir * dirPrimCenter )*180./M_PI;

        dir = c2 - c1;
        dir.normalize();
        osg::Vec3 closestPoint;
        osg::Vec3 normal = _alphaShapes.at( iter->timeStep-1 )->getSurfaceNormal( divPos * _inverseRotMatrix,
                                                                                  closestPoint );
        closestPoint = closestPoint * _rotMatrix;
//        if( _sliderType == 0 )
//          this->renderClosestPointAndLine( closestPoint, divPos, iter->timeStep-1 );
//        else if( _sliderType == 1 )
//          this->renderClosestPointAndLine( closestPoint, divPos, iter->getCellCycleId() );

        double distSurfaceDiv = (closestPoint - divPos).length();

        if( distSurfaceDiv < 0.0001 )
          distSurfaceDiv = 0.;

        double angle3 = acos( dir * normal )*180./M_PI;
        if( angle3 > 90. )
          angle3 = 180. - angle3;

        out << iter->cellId << " "
            << iter->timeStep << " "
            << iter->getCellCycleId() << " "
            << divPos.x() << " " << divPos.y() << " " << divPos.z() << " "
            << c1.x() << " " << c1.y() << " " << c1.z() << " "
            << c2.x() << " " << c2.y() << " " << c2.z() << " "
            << cen.x() << " " << cen.y() << " " << cen.z() << " "
            << angle1 << " " << (divPos - cen).length() << " "
            << angle2 << " " << (divPos - primCen).length() << " "
            << angle3 << " " << distSurfaceDiv << "\n";

        // store angle values
        Angles ang;
        ang.VLRangle = angle1;
        this->checkMinMax( angle1, 0 );
        ang.primAngle = angle2;
        this->checkMinMax( angle2, 1 );
        ang.surfaceAngle = angle3;
        this->checkMinMax( angle3, 2 );
        ang.divPos = divPos;
        ang.cellCycle = iter->getCellCycleId();
        ang.numberOfCells = _cellsPerTimeStep.at( iter->timeStep-1 ).size();
        ang.id = iter->cellId;
        ang.timeStep = iter->timeStep-1;
        _angleValues.insert( std::make_pair( std::make_pair( iter->cellId, iter->timeStep ), ang ) );
        _divisionProperties.at( iter->timeStep-1 ).insert( ang );
      }
    }
  }
}

//----------------------------------------------------------------------------

void SDivisionAnalysis::readAngleData()
{
  // TODO: has to be tested!!

  std::string fileName = "/home/necrolyte/Uni/LateralRootGrowth/261115/DivisionOrientations_LateTimeStep_29102015_";
  fileName += _lineages->getName();
  fileName += ".csv";
  std::ifstream in( fileName.c_str(), std::ifstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open data " + fileName );
    return;
  }

  std::string line;
  getline( in, line );

  while( in.good() )
  {
    getline( in, line );
    if( line == "" )
      break;

    std::string dummy;
    std::size_t id, t;
    double angle;
    std::stringstream lineStream( line );
    lineStream >> id >> t;

    std::set<Angles, setComp>::iterator iter = _divisionProperties.at( t-1 ).begin();
    for( ; iter != _divisionProperties.at( t-1 ).end(); ++iter )
    {
      if( iter->id == id )
        break;
    }

    for( std::size_t e = 0; e < 16; e++ )
      lineStream >> dummy;

    std::string layer = "";
    std::size_t numberOfCells, lineage, divType;
    int file;
    double theta, phi, tau;

    lineStream >> numberOfCells >> lineage >> file;
    iter->numberOfCells = numberOfCells;
    iter->cellFile = file;
    iter->lineage = lineage;

    lineStream >> dummy;
    layer += dummy;
    if( dummy.at(0) == '\"' )
    {
      while( *dummy.rbegin() != '\"' )
      {
        lineStream >> dummy;
        layer += dummy;
      }

      layer += dummy;
    }

    std::cout << "layer: " << layer << std::endl;
    iter->layerInfo = layer;

    lineStream >> divType >> theta >> phi >> tau;
    this->checkMinMax( theta, 3 );
    this->checkMinMax( phi, 4 );
    this->checkMinMax( tau, 5 );
    // requires all angles to be mutable in the struct definition
    iter->theta = theta;
    iter->phi = phi;
    iter->tau = tau;
  }

  in.close();
}

//----------------------------------------------------------------------------

void SDivisionAnalysis::performCylindricalAngleAnalysis()
{
  std::string fileName = "/tmp/CyclindricalCoords_" + _lineages->getName() + ".csv";
  // open file at the end of the stream and allow writing
  std::ofstream out( fileName.c_str(), std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to write into data " + fileName );
    return;
  }

  osg::Vec3 domeTip;
  // center is the dome tip of the last time step
  if( _lineages->getName() == "120830" )
    domeTip = osg::Vec3( 277.2, 215.1, 72.3 );
  else if( _lineages->getName() == "121204" )
    domeTip = osg::Vec3( 347.0, 300.1, 41.3 );
  else if( _lineages->getName() == "121211" )
    domeTip = osg::Vec3( 278.2, 251.8, 145.6 );
  else if( _lineages->getName() == "130508" )
    domeTip = osg::Vec3( 248.6, 144.5, 337.9 );
  else if( _lineages->getName() == "130607" )
    domeTip = osg::Vec3( 160.1, 388.0, 298.0 );

  // VLR center
  SArray ce = _lineages->getCenter();
  osg::Vec3 VLRCenter( ce[0], ce[1], ce[2] );

  // header
  out << "DomeTip: " << domeTip.x() << " " << domeTip.y() << " " << domeTip.z() << "\n";
  out << "VLRCenter: " << VLRCenter.x() << " " << VLRCenter.y() << " " << VLRCenter.z() << "\n";
  out << "CellID TimeStep CellCycle DivPosX DivPosY DivPosZ";
  out << " CenterOfMassX CenterOfMassY CenterOfMassZ "
      << "m_D Rho_D Theta_D m_M Rho_M Theta_M m_C Rho_C Theta_C\n";

  std::size_t totalLineages = _lineages->size();
  std::size_t li = 0;
  for( SAdvancedLineageTreeCollection::const_iterator l = _lineages->begin();
       l != _lineages->end(); ++l, li++ )
  {
    SKernel::get()->setProgress( (double)li/(double)(totalLineages - 1),
                                 _lineages->getName() + ": Angle computations" );

    for( SLineageTree::const_iterator iter = l->second->begin();
         iter != l->second->end(); ++iter )
    {
      if( iter->children.size() == 2 )
      {
        osg::Vec3 divPos( iter->getX(), iter->getY(), iter->getZ() );
        divPos = divPos * _rotMatrix;

        std::vector<osg::Vec3> daughters;

        // first child
        SLineageTree::const_iterator iter1 = iter->children[0];
        while( iter1->children.size() != 2 )
        {
          if( iter1->children.size() == 0 )
            break;
          else
            ++iter1;
        }

        osg::Vec3 c1 = osg::Vec3( iter1->getX(), iter1->getY(), iter1->getZ() );
        c1 = c1 * _rotMatrix;

        daughters.push_back( c1 );

        // second child
        SLineageTree::const_iterator iter2 = iter->children[1];
        while( iter2->children.size() != 2 )
        {
          if( iter2->children.size() == 0 )
            break;
          else
            ++iter2;
        }

        osg::Vec3 c2 = osg::Vec3( iter2->getX(), iter2->getY(), iter2->getZ() );
        c2 = c2 * _rotMatrix;

        daughters.push_back( c2 );

//        osg::Vec3 c1( iter->children[0]->getX(),
//            iter->children[0]->getY(),
//            iter->children[0]->getZ() );

//        osg::Vec3 c2( iter->children[1]->getX(),
//            iter->children[1]->getY(),
//            iter->children[1]->getZ() );

        osg::Vec3 cen = _centresOfMass.at(iter->timeStep-1);
        cen = cen * _rotMatrix;

        for( std::size_t c = 0; c < 2; c++ )
        {
          std::vector<double> m;
          std::vector<double> rho;
          std::vector<double> theta;

          // set the dome tip to the origin
          osg::Vec3 p = daughters.at(c);
          p = p - domeTip;
          m.push_back( daughters.at(c).z() - domeTip.z() );
          rho.push_back( std::sqrt( p.x()*p.x() + p.y()*p.y() ) );
          theta.push_back( std::atan2( p.y(), p.x() ) * 180./M_PI );

          // set the center of mass to the origin
          p = daughters.at(c);
          p = p - cen;
          m.push_back( daughters.at(c).z() - cen.z() );
          rho.push_back( std::sqrt( p.x()*p.x() + p.y()*p.y() ) );
          theta.push_back( std::atan2( p.y(), p.x() ) * 180./M_PI );

          // set the VLR center to the origin
          p = daughters.at(c);
          p = p - VLRCenter;
          m.push_back( daughters.at(c).z() - VLRCenter.z() );
          rho.push_back( std::sqrt( p.x()*p.x() + p.y()*p.y() ) );
          theta.push_back( std::atan2( p.y(), p.x() ) * 180./M_PI );

          out << iter->cellId << " "
              << iter->timeStep << " "
              << iter->getCellCycleId() << " "
              << divPos.x() << " " << divPos.y() << " " << divPos.z() << " "
              << cen.x() << " " << cen.y() << " " << cen.z() << " "
              << m.at(0) << " " << rho.at(0) << " " << theta.at(0) << " "
              << m.at(1) << " " << rho.at(1) << " " << theta.at(1) << " "
              << m.at(2) << " " << rho.at(2) << " " << theta.at(2) << "\n";
        }
      }
    }
  }
}

//----------------------------------------------------------------------------

void SDivisionAnalysis::checkMinMax( const double value, const int type )
{
  if( value > _minMaxAngleValues.at(type).second )
    _minMaxAngleValues.at(type).second = value;

  if( value <= _minMaxAngleValues.at(type).first )
    _minMaxAngleValues.at(type).first = value;
}

//----------------------------------------------------------------------------

void SDivisionAnalysis::renderClosestPointAndLine( const osg::Vec3 &p1,
                                                   const osg::Vec3 &p2,
                                                   const std::size_t t )
{
  if( _sliderType == 0 )
  {
    // first render closest point
    SInteractiveCell* closestPoint = new SInteractiveCell( 8., 0, t+1, osg::Vec4( 1., 1., 0., 1. ) );
    closestPoint->addPosition( t+1, p1 );
    _renderGroup->addChild( closestPoint );
  }

  // then the connection between closest point and division center
  SSimilarityMeasureGraphics::addCylinderBetweenPoints( p1, p2,
                                                    0.5, osg::Vec4( 0.8, 0.8, 0.8, 1. ),
                                                    _surfaceGroup.at(_currentView).at(t) );
}

//----------------------------------------------------------------------------
