#include "SCellLayers.hh"

/**
  @file   SCellLayers.cc
  @brief  Contains class for generating cell layers of arabidopsis data
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SLayerInformationIO.hh"
#include "SSimilarityMeasureGraphics.hh"

#include "kernel/SKernel.hh"

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>

#include <utility>

// ---------------------------------------------------------------------

SCellLayers::SCellLayers( const boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                          const double changeLayerThreshold,
                          const osg::Matrix &inverseMatrix,
                          const osg::ref_ptr<osg::Vec4Array> layerColors,
                          const float arrowParam1,
                          const float arrowParam2,
                          const bool wireframe,
                          const int shapeType )
  : _currentNumLayers( 1 ),
    _maxNumLayers( std::numeric_limits<unsigned int>::min() ),
    _lineages( lineages ),
    _changeLayerThreshold( changeLayerThreshold ),
    _inverseMatrix( inverseMatrix ),
    _layerColors( layerColors ),
    _arrowParam1( arrowParam1 ),
    _arrowParam2( arrowParam2 ),
    _wireframe( wireframe ),
    _onlyLayerInfo( false ),
    _shapeType( shapeType ),
    _timerForGeometryGeneration( 0. ),
    _timerForLayerGeneration( 0. )
{
  // at first allocate data types to store layer information
  // and get all cells per time step
  this->determineCellsPerTimeStep();
}

// ---------------------------------------------------------------------

SCellLayers::SCellLayers( const boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                          const double changeLayerThreshold )
  : _currentNumLayers( 1 ),
    _maxNumLayers( std::numeric_limits<unsigned int>::min() ),
    _lineages( lineages ),
    _changeLayerThreshold( changeLayerThreshold ),
    _inverseMatrix( osg::Matrix::identity() ),
    _arrowParam1( 1.0 ),
    _arrowParam2( 1.0 ),
    _wireframe( false ),
    _onlyLayerInfo( true ),
    _shapeType( 0 ),
    _timerForGeometryGeneration( 0. ),
    _timerForLayerGeneration( 0. )
{
  // dummy color
  osg::ref_ptr<osg::Vec4Array> layerColors;
  layerColors = new osg::Vec4Array;
  layerColors->push_back( osg::Vec4( 0., 0., 0., 1. ) );
  _layerColors = layerColors;

  // at first allocate data types to store layer information
  // and get all cells per time step
  this->determineCellsPerTimeStep();
}

// ---------------------------------------------------------------------

void SCellLayers::initializeDataTypes()
{
  _cellsPerTimeStep = boost::shared_ptr<cellTimeVector>( new cellTimeVector() );
  _cellLayersPerTimeStep = boost::shared_ptr<cellLayerVector>( new cellLayerVector() );
  _cellHistories = boost::shared_ptr< std::vector<cellHistoryMap> >( new std::vector<cellHistoryMap>() );

  if( _onlyLayerInfo )
  {
    _normalDivisionMap = boost::shared_ptr< std::map<const SLineageTree*, osg::Vec3> >(
          new std::map<const SLineageTree*, osg::Vec3>() );
  }

  // resize cells per time step vector by the number of time steps
  _cellsPerTimeStep->resize( _lineages->getOrComputeNumberOfTimesteps() );

  // resize outer layer per time step vector to the same size
  _cellLayersPerTimeStep->resize( _cellsPerTimeStep->size() );
  // and allocate memory for the first layer
  _cellLayersPerTimeStep->front().resize( 1 );

  // resize cell history vector by one for the first time step
  _cellHistories->resize( 1 );
}

// ---------------------------------------------------------------------

void SCellLayers::determineCellsPerTimeStep()
{
  // first initialize all required data types
  this->initializeDataTypes();

  for( auto l = _lineages->begin(); l != _lineages->end(); ++l )
  {
    SLineageTree *tree = l->second;
    for( SLineageTree::const_iterator nodeIt = tree->begin();
         nodeIt != tree->end(); ++nodeIt )
    {
      _cellsPerTimeStep->at(nodeIt->timeStep - 1).insert( *nodeIt );

      if( nodeIt->timeStep == 1 )
      {
        _cellLayersPerTimeStep->front().back().insert( *nodeIt );
        cellHistory cs;
        cs.insert( 0 );
        _cellHistories->back().insert( std::make_pair( *nodeIt, cs ) );
      }
    }
  }
}

// ---------------------------------------------------------------------

void SCellLayers::determineLayers( osg::ref_ptr<osg::MatrixTransform> addToGroup,
                                   osg::ref_ptr<SSliderCallback> sliderCallback,
                                   const bool layerInfoLoaded,
                                   const bool renderDelaunay )
{
  if( _onlyLayerInfo )
  {
    std::cout << "Wrong constructor is used!" << std::endl;
    return;
  }

  // interactive delaunay triangulation class in order to update
  // all layers for each time step
  SInteractiveSurface *is = new SInteractiveSurface();
  is->setUpdateCallback( sliderCallback );
  is->setName( "Delaunay" );

  // at the beginning we add an initial amount of assumed layers
  // which we will fill and update after traversing all time steps
  // at the end we clear all layer groups which are empty
  // not the best method but required since the scene graph info
  // is not updated when interacting with the slider -> TODO
  int initialLayers = 10;

  // if the layer information is loaded then we already know
  // the exact number of layers
  if( layerInfoLoaded )
    initialLayers = _maxNumLayers;

  for( int i=0;i<initialLayers;i++ )
  {
    osg::ref_ptr<osg::Group> initialLayerGroup = new osg::Group;
    initialLayerGroup->setName( "Layer " + boost::lexical_cast<std::string>(i) );
    is->addChild( initialLayerGroup );
  }

  if( renderDelaunay )
    addToGroup->addChild( is );

  // initialize normal switch for each time step
  SInteractiveNormals *in = new SInteractiveNormals();
  in->setUpdateCallback( sliderCallback );
  in->setName( "Normals" );
  addToGroup->addChild( in );

  // clear volume vector before generating layers
  // and resize to number of time steps
  _layerVolumes.clear();
  _layerVolumes.resize( _cellsPerTimeStep->size() );

  // travere all time steps and generate the cell layers
  // beginning with time step 0
  this->generateCellLayers( 0, is, in, false, layerInfoLoaded );

  // at the end remove all children which are empty
  if( !layerInfoLoaded )
    is->removeChild( _maxNumLayers, is->getNumChildren() - _maxNumLayers );
}

// ---------------------------------------------------------------------

void SCellLayers::determineLayersWithoutGeometry( const bool layerInfoLoaded )
{
  if( !_onlyLayerInfo )
  {
    std::cout << "Wrong constructor is used!" << std::endl;
    return;
  }

  // we create dummies for storing the delaunay triangulation
  // and the normals but in fact we do not assign them to any osg group
  SInteractiveSurface *is = new SInteractiveSurface();
  SInteractiveNormals *in = new SInteractiveNormals();

  // clear volume vector before generating layers
  // and resize to number of time steps
  _layerVolumes.clear();
  _layerVolumes.resize( _cellsPerTimeStep->size() );

  // travere all time steps and generate the cell layers
  // beginning with time step 0
  this->generateCellLayers( 0, is, in, false, layerInfoLoaded );
}

// ---------------------------------------------------------------------

void SCellLayers::updateLayers( const unsigned int startTimeStep,
                                SInteractiveSurface *is,
                                SInteractiveNormals *in )
{
  if( _onlyLayerInfo )
  {
    std::cout << "Wrong constructor is used!" << std::endl;
    return;
  }

  // clear volume vector before generating layers
  // and resize to number of time steps
  _layerVolumes.clear();
  _layerVolumes.resize( _cellsPerTimeStep->size() );

  // travere all time steps and regenerate the cell layers
  // beginning with startTimeStep
  this->generateCellLayers( startTimeStep, is, in, true, false );
}

// ---------------------------------------------------------------------

void SCellLayers::generateCellLayers( const unsigned int startTimeStep,
                                      SInteractiveSurface *is,
                                      SInteractiveNormals *in,
                                      const bool update,
                                      const bool layerInfoLoaded )
{
  // starting time step
  unsigned int timeStep = 0;

  // for each time step
  BOOST_FOREACH( cellSet cS, *_cellsPerTimeStep )
  {
    double progress = (double)timeStep/(double)(_cellsPerTimeStep->size()-1);

    if( !update )
      SKernel::get()->setProgress( progress, "Layer generation" );
    else
      SKernel::get()->setProgress( progress, "Update layers" );

    // skip all previous time steps until we reach the cell
    // that will be assigned to a new layer
    if( timeStep + 1 < startTimeStep )
    {
      timeStep++;
      continue;
    }

    // check the current layer sizes
    // if a layer has no cells in it then erase it
    // -> this can occur if cells were put into a layer
    // in time step t but vanish in time step t+1
    this->checkLayerSizes( timeStep );

    // get the current layer size
    unsigned int curLayerSize = _cellLayersPerTimeStep->at(timeStep).size();

    if( timeStep % 100 == 0 )
    {
      std::cout << "t: " << timeStep
                << " numLayers: " << curLayerSize << std::endl;
    }

    // set update callbacks and initialize groups per time step
    osg::ref_ptr<osg::Group> normalGroupPerTimeStep = new osg::Group;
    normalGroupPerTimeStep->setName( "Normals for t= " + boost::lexical_cast<std::string>( timeStep+1 ) );

    // initialize new cell history for the future cells
    cellHistoryMap currentCellHistoryMap;

    // and at the beginning of each new time step
    // set layer size to the same size as before
    // which could be updated during division check
    // except for the last time step
    if( !layerInfoLoaded )
    {
      if( timeStep != _cellLayersPerTimeStep->size()-1 )
        _cellLayersPerTimeStep->at(timeStep+1).resize( curLayerSize );
    }

    // this vector stores the point size for all
    // layer values; this is important in order to
    // switch between the normal computation
    std::vector<unsigned int> cellsInLayer;
    cellsInLayer.resize( curLayerSize );

    // delaunay triangulation of the current time step and layers
    // convex hull
    std::vector<boost::shared_ptr<SConvexHull> > currentCHVector;
    std::vector<boost::shared_ptr<SAlphaShape> > currentASVector;
    if( _shapeType == 0 )
      currentCHVector.resize( curLayerSize, boost::shared_ptr<SConvexHull>( new SConvexHull() ) );
    else if( _shapeType == 1 )
      currentASVector.resize( curLayerSize, boost::shared_ptr<SAlphaShape>( new SAlphaShape() ) );

    // generate geometries for each layer,
    // e.i. lines, triangles or surface triangulation
    if( _shapeType == 0 )
    {
      this->generateLayerGeometry( timeStep, is, cS, cellsInLayer, currentCHVector,
                                   currentCellHistoryMap, layerInfoLoaded );
    }
    else if( _shapeType == 1 )
    {
      this->generateLayerGeometry( timeStep, is, cS, cellsInLayer, currentASVector,
                                   currentCellHistoryMap, layerInfoLoaded );
    }

    // stop loop if the last time step is reached since we always
    // set the layer for the next time step which is not possible
    // for the last one
    if( !layerInfoLoaded )
    {
      if( timeStep == _cellLayersPerTimeStep->size()-1 )
        break;
    }

    // for each cell in a time step
    BOOST_FOREACH( cellSet::value_type cellInTimeStep, cS )
    {
      boost::timer timer;

      // determine in which layer the current cell is
      unsigned int layerValue = this->determineLayerValue( cellInTimeStep );

      // only handle division cases
      if( cellInTimeStep->children.size() == 2 )
      {
        osg::Vec3 cellPos( cellInTimeStep->getX(),
                           cellInTimeStep->getY(),
                           cellInTimeStep->getZ() );

        osg::Vec3 normal;

        // compute and render vertex normals in current delaunay triangulation
        if( cellsInLayer[layerValue] > 3 )
        {
          if( _shapeType == 0 )
            normal = currentCHVector[layerValue]->getVertexNormal( cellPos );
          else if( _shapeType == 1 )
            normal = currentASVector[layerValue]->getVertexNormal( cellPos );

          // handle inner nodes in the convex hull which can only
          // occur in the CGAL delaunay triangulation
          if( normal == osg::Vec3( 0., 0., 0. ) )
          {
            osg::Vec3 closestPoint;

            if( _shapeType == 0 )
              normal = currentCHVector[layerValue]->getSurfaceNormal( cellPos, closestPoint );
            else if( _shapeType == 1 )
              normal = currentASVector[layerValue]->getSurfaceNormal( cellPos, closestPoint );
          }
        }
        else if( cellsInLayer[layerValue] == 1 || cellsInLayer[layerValue] == 2 )
          normal = this->getCenterNormal( cellPos );
        else if( cellsInLayer[layerValue] == 3 )
        {
          osg::Vec3Array *points = new osg::Vec3Array;

          BOOST_FOREACH( cellSet::value_type cell, _cellLayersPerTimeStep->at(timeStep).at(layerValue) )
          {
            osg::Vec3 p( cell->getX(), cell->getY(), cell->getZ() );
            points->push_back( p );
          }

          // direction vectors from current cell to cells already in layer
          osg::Vec3 layerDir1 = points->at(1) - points->at(0);
          osg::Vec3 layerDir2 = points->at(2) - points->at(0);

          // get normal vector from cross product
          normal = layerDir1 ^ layerDir2;
          normal.normalize();
        }
        // render normal
        SSimilarityMeasureGraphics::renderArrow( cellPos, normal, normalGroupPerTimeStep,
                                             osg::Vec4( 0., 0., 1., 1. ) );

//        osg::Vec3 c1( cellInTimeStep->children.at(0)->getX(),
//                      cellInTimeStep->children.at(0)->getY(),
//                      cellInTimeStep->children.at(0)->getZ() );
//        osg::Vec3 c2( cellInTimeStep->children.at(1)->getX(),
//                      cellInTimeStep->children.at(1)->getY(),
//                      cellInTimeStep->children.at(1)->getZ() );
//        osg::Vec3 cDir = c2 - c1;
//        cDir.normalize();

//        // render division direction
//        SSimilarityMeasureGraphics::renderArrow( cellPos, cDir, normalGroupPerTimeStep,
//                                             osg::Vec4( 0., 1., 0., 1. ) );
//        SSimilarityMeasureGraphics::renderArrow( cellPos, -cDir, normalGroupPerTimeStep,
//                                             osg::Vec4( 0., 1., 0., 1. ) );

        // if only the layer information is required then additionally
        // store the normal vector for each division node which is later
        // used for the identification of periclinal or anticlinal divisions
        if( _onlyLayerInfo )
          (*_normalDivisionMap)[cellInTimeStep] = normal;

        // if the layer info is not loaded then start the testing
        // for each division node
        if( !layerInfoLoaded )
        {
          // check both children for new layer generation
          int childChoice = this->changeLayerCheck( cellInTimeStep, normal );

          // check if the current cell generates a new layer or should
          // be assigned to an already existing one
          if( childChoice != 0 )
          {
            // layer assignment
            this->assignNewLayer( currentCellHistoryMap, cellInTimeStep, childChoice, layerValue+1 );
          }
          // none of the children have moved to a new layer
          // assign both children the layer value from its parent
          else
          {
            _cellLayersPerTimeStep->at(timeStep+1).at(layerValue).insert( cellInTimeStep->children[0] );
            this->updateCellHistory( currentCellHistoryMap, cellInTimeStep->children[0], layerValue );
            _cellLayersPerTimeStep->at(timeStep+1).at(layerValue).insert( cellInTimeStep->children[1] );
            this->updateCellHistory( currentCellHistoryMap, cellInTimeStep->children[1], layerValue );
          }
        }
        // else only update the cell history
        else
        {
          // get assigned layer value from _cellLayersPerTimeStep
          unsigned int layerValue1 = this->determineLayerValue( cellInTimeStep->children[0] );
          unsigned int layerValue2 = this->determineLayerValue( cellInTimeStep->children[1] );

          this->updateCellHistory( currentCellHistoryMap, cellInTimeStep->children[0], layerValue1 );
          this->updateCellHistory( currentCellHistoryMap, cellInTimeStep->children[1], layerValue2 );
        }
      }
      // else cell movements stay in the same layer as before
      else if( cellInTimeStep->children.size() == 1 )
      {
        if( !layerInfoLoaded )
          _cellLayersPerTimeStep->at(timeStep+1).at(layerValue).insert( cellInTimeStep->children[0] );

        // in both cases update the cell history
        this->updateCellHistory( currentCellHistoryMap, cellInTimeStep->children[0], layerValue );
      }

      _timerForLayerGeneration += timer.elapsed();
    }

    // add normal group to switch node
    in->addNormals( timeStep+1, normalGroupPerTimeStep );

    timeStep++;

    // update the cell histories of the future cells
    _cellHistories->push_back( currentCellHistoryMap );

    if( !layerInfoLoaded )
    {
      // update max number of layers
      if( _currentNumLayers > _maxNumLayers )
        _maxNumLayers = _currentNumLayers;
    }
  }

  // at the end output the volume information
  SLayerInformationIO::writeVolumes( "/tmp/layerVolumeInformation.csv", _cellLayersPerTimeStep,
                                     _layerVolumes, _maxNumLayers );
}

// ---------------------------------------------------------------------

void SCellLayers::generateLayerGeometry( const unsigned int timeStep,
                                         SInteractiveSurface *is,
                                         const cellSet cS,
                                         std::vector<unsigned int> &cellsInLayer,
                                         std::vector<boost::shared_ptr<SConvexHull> > &currentDTVector,
                                         cellHistoryMap &currentCellHistoryMap,
                                         const bool layerInfoLoaded )
{
  boost::timer timer;

  unsigned int layerCounter = 0;

  // resize volume vector to the number of layers existing
  // in this time step
  _layerVolumes.at(timeStep).resize( _cellLayersPerTimeStep->at(timeStep).size() );

  // for each layer generate a delaunay triangulation
  BOOST_FOREACH( cellSet &cSL, _cellLayersPerTimeStep->at(timeStep) )
  {
    osg::ref_ptr<osg::Group> layerGroup = new osg::Group;
    osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array;

    // for each cell in a layer
    BOOST_FOREACH( cellSet::value_type cellInLayer, cSL )
    {
      osg::Vec3 coord = osg::Vec3 ( cellInLayer->getX(),
                                    cellInLayer->getY(),
                                    cellInLayer->getZ() );
      points->push_back( coord );

      if( layerInfoLoaded )
      {
        if( cellInLayer->parent )
        {
          // update cell history
          this->updateCellHistory( currentCellHistoryMap,
                                   cellInLayer,
                                   layerCounter );
        }
      }
    }

    if( layerCounter == 0 )
    {
      // check in _cellsPerTimeStep if some cells begin later
      // if so then also add them to points and _cellLayersPerTimeStep
      // for each cell in a time step
      // at the moment we assign those cells only to the initial layer -> TODO
      BOOST_FOREACH( cellSet::value_type cellInTimeStep, cS )
      {
        // if cell is a root and does not begin at time step 0
        if( !(cellInTimeStep->parent) && timeStep != 0 )
        {
          osg::Vec3 coord = osg::Vec3 ( cellInTimeStep->getX(),
                                        cellInTimeStep->getY(),
                                        cellInTimeStep->getZ() );
          points->push_back( coord );

          cSL.insert( cellInTimeStep );

          // insert those later cells also in the cell history
          cellHistory cs;
          cs.insert( 0 );
          _cellHistories->back().insert( std::make_pair( cellInTimeStep, cs ) );
        }
      }
    }

    // initialize color depending on layer counter
    osg::Vec4 color( (*_layerColors)[layerCounter%_layerColors->size()] );

    // there is never a case of having zero cells
    // and there is nothing to do in the case of having only one cell
    // if only two points are inside the layer then draw a cylinder between them
    if( points->size() == 2 )
    {
      SSimilarityMeasureGraphics::addCylinderBetweenPoints( points->at(0), points->at(1), 0.2, color, layerGroup );
      _layerVolumes.at(timeStep).at(layerCounter) = 0.;
    }
    // if there are three points then draw a triangle
    // which we only need when we generate delaunay with CGAL
    else if( points->size() == 3 )
    {
      SSimilarityMeasureGraphics::drawTriangle( points, color, layerGroup );
      _layerVolumes.at(timeStep).at(layerCounter) = 0.;
    }
    else if( points->size() > 3 )
    {
      boost::shared_ptr<SConvexHull> dT = boost::shared_ptr<SConvexHull>(
            new SConvexHull( points, timeStep, color, true ) );
      currentDTVector[layerCounter] = dT;
      dT->wireFrame( _wireframe );
      layerGroup->addChild( dT->getTriangles() );

      //_layerVolumes.at(timeStep).at(layerCounter) = dT->computeVolume();
    }
    // else there is only a single cell
    else
      _layerVolumes.at(timeStep).at(layerCounter) = 0.;

    // store point size
    cellsInLayer[layerCounter] = points->size();

    // set triangulation of layers for each time step
    is->addTriangulation( timeStep+1, layerCounter, layerGroup );

    // afterwards increment layer counter
    layerCounter++;
  }

  _timerForGeometryGeneration += timer.elapsed();
}

// ---------------------------------------------------------------------

void SCellLayers::generateLayerGeometry( const unsigned int timeStep,
                                         SInteractiveSurface *is,
                                         const cellSet cS,
                                         std::vector<unsigned int> &cellsInLayer,
                                         std::vector<boost::shared_ptr<SAlphaShape> > &currentDTVector,
                                         cellHistoryMap &currentCellHistoryMap,
                                         const bool layerInfoLoaded )
{
  boost::timer timer;

  unsigned int layerCounter = 0;

  // resize volume vector to the number of layers existing
  // in this time step
  _layerVolumes.at(timeStep).resize( _cellLayersPerTimeStep->at(timeStep).size() );

  // for each layer generate a delaunay triangulation
  BOOST_FOREACH( cellSet &cSL, _cellLayersPerTimeStep->at(timeStep) )
  {
    osg::ref_ptr<osg::Group> layerGroup = new osg::Group;
    osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array;

    // for each cell in a layer
    BOOST_FOREACH( cellSet::value_type cellInLayer, cSL )
    {
      osg::Vec3 coord = osg::Vec3 ( cellInLayer->getX(),
                                    cellInLayer->getY(),
                                    cellInLayer->getZ() );
      points->push_back( coord );

      if( layerInfoLoaded )
      {
        if( cellInLayer->parent )
        {
          // update cell history
          this->updateCellHistory( currentCellHistoryMap,
                                   cellInLayer,
                                   layerCounter );
        }
      }
    }

    if( layerCounter == 0 )
    {
      // check in _cellsPerTimeStep if some cells begin later
      // if so then also add them to points and _cellLayersPerTimeStep
      // for each cell in a time step
      // at the moment we assign those cells only to the initial layer -> TODO
      BOOST_FOREACH( cellSet::value_type cellInTimeStep, cS )
      {
        // if cell is a root and does not begin at time step 0
        if( !(cellInTimeStep->parent) && timeStep != 0 )
        {
          osg::Vec3 coord = osg::Vec3 ( cellInTimeStep->getX(),
                                        cellInTimeStep->getY(),
                                        cellInTimeStep->getZ() );
          points->push_back( coord );

          cSL.insert( cellInTimeStep );

          // insert those later cells also in the cell history
          cellHistory cs;
          cs.insert( 0 );
          _cellHistories->back().insert( std::make_pair( cellInTimeStep, cs ) );
        }
      }
    }

    // initialize color depending on layer counter
    osg::Vec4 color( (*_layerColors)[layerCounter%_layerColors->size()] );

    // there is never a case of having zero cells
    // and there is nothing to do in the case of having only one cell
    // if only two points are inside the layer then draw a cylinder between them
    if( points->size() == 2 )
    {
      SSimilarityMeasureGraphics::addCylinderBetweenPoints( points->at(0), points->at(1), 0.2, color, layerGroup );
      _layerVolumes.at(timeStep).at(layerCounter) = 0.;
    }
    // if there are three points then draw a triangle
    // which we only need when we generate delaunay with CGAL
    else if( points->size() == 3 )
    {
      SSimilarityMeasureGraphics::drawTriangle( points, color, layerGroup );
      _layerVolumes.at(timeStep).at(layerCounter) = 0.;
    }
    else if( points->size() > 3 )
    {
      boost::shared_ptr<SAlphaShape> dT = boost::shared_ptr<SAlphaShape>(
            new SAlphaShape( points, timeStep, color, true ) );
      currentDTVector[layerCounter] = dT;
      dT->wireFrame( _wireframe );
      layerGroup->addChild( dT->getTriangles() );

      _layerVolumes.at(timeStep).at(layerCounter) = dT->computeVolume();
    }
    // else there is only a single cell
    else
      _layerVolumes.at(timeStep).at(layerCounter) = 0.;

    // store point size
    cellsInLayer[layerCounter] = points->size();

    // set triangulation of layers for each time step
    is->addTriangulation( timeStep+1, layerCounter, layerGroup );

    // afterwards increment layer counter
    layerCounter++;
  }

  _timerForGeometryGeneration += timer.elapsed();
}

// ---------------------------------------------------------------------

bool SCellLayers::assignNewLayerValue( const unsigned int layerValue,
                                       const int cellId,
                                       const std::size_t timeStep,
                                       unsigned int &oldLayerValue )
{  
  oldLayerValue = 0;

  // for each layer
  BOOST_FOREACH( cellSet cSIL, _cellLayersPerTimeStep->at(timeStep-1) )
  {
    // find selected cell in each set of a layer
    BOOST_FOREACH( cellSet::value_type cell, cSIL )
    {
      if( cell->cellId == cellId && cell->timeStep == timeStep )
      {
        // check if the selected node is a child of a division
        // else we do not allow a layer assignment which only
        // makes sense in the case of considering a cell division
        // the same applies for root nodes
        if( !cell->parent )
        {
          std::cout << "The cell is not a child of a division." << std::endl;
          return false;
        }

        if( cell->parent )
        {
          if( cell->parent->children.size() != 2 )
          {
            std::cout << "The cell is not a child of a division." << std::endl;
            return false;
          }
        }

        // check, if new layer assignment is even necessary
        if( layerValue == oldLayerValue )
        {
          std::cout << "The cell is already in the chosen layer." << std::endl;
          return false;
        }
        else
        {
          // before setting the new layer value we check of we have to
          // resize the cell layer vector since we are allowed to assign
          // a layer value for which the layer vector did not allocate
          // any memory at the moment
          if( _cellLayersPerTimeStep->at(timeStep-1).size() <= layerValue )
            _cellLayersPerTimeStep->at(timeStep-1).resize( layerValue+1 );

          _cellLayersPerTimeStep->at(timeStep-1).at(layerValue).insert( cell );
          _cellLayersPerTimeStep->at(timeStep-1).at(oldLayerValue).erase( cell );

          // copy the cell history of the parent cell if it exists
          // and add the new assigned layer value
          if( cell->parent )
          {
            std::map<const SLineageTree*, cellHistory>::iterator cIter =
                _cellHistories->at(timeStep-1).find( cell );

            std::map<const SLineageTree*, cellHistory>::iterator pIter =
                _cellHistories->at(timeStep-2).find( cell->parent );

            if( cIter != _cellHistories->at(timeStep-1).end() &&
                pIter != _cellHistories->at(timeStep-2).end() )
            {
              cIter->second = pIter->second;
              cIter->second.insert( layerValue );
            }
          }
          return true;
        }
      }
    }

    oldLayerValue++;
  }

  return false;
}

// ---------------------------------------------------------------------

void SCellLayers::clearLayerInformation( const unsigned int startTimeStep )
{
  std::vector<cellHistoryMap>::iterator historyIter = _cellHistories->begin();
  std::vector< std::vector<cellSet> >::iterator layerIter = _cellLayersPerTimeStep->begin();

  for( unsigned int count = 0;count < startTimeStep;count++ )
  {
    historyIter++;
    layerIter++;
  }

  // erase cell histories per time step
  _cellHistories->erase( historyIter, _cellHistories->end() );
  // erase cell layer information per time step
  _cellLayersPerTimeStep->erase( layerIter, _cellLayersPerTimeStep->end() );
  _cellLayersPerTimeStep->resize( _cellsPerTimeStep->size() );
  // update the current and maximum number of layers
  // by setting the layer size of the previous time step
  _currentNumLayers = (*(layerIter-1)).size();
  _maxNumLayers = _currentNumLayers;

  // check layers of they are empty which can occur if a node was
  // selected that was the only representant of a layer
  // note that here also the current number of layers can be changed
  if( startTimeStep > 0 )
    this->checkLayerSizes( startTimeStep-1 );
}

// ---------------------------------------------------------------------

int SCellLayers::changeLayerCheck( const SLineageTree *parent,
                                   const osg::Vec3 &pNormal )
{
  // childChoice indicates which cell child changes its layer
  // 0: none
  // 1: first
  // 2: second
  int childChoice = 0;

  osg::Vec3 cPos1( parent->children[0]->getX(),
                   parent->children[0]->getY(),
                   parent->children[0]->getZ() );
  osg::Vec3 cPos2( parent->children[1]->getX(),
                   parent->children[1]->getY(),
                   parent->children[1]->getZ() );

  // get direction vector between both cells
  osg::Vec3 childrenDir = cPos1 - cPos2;
  childrenDir.normalize();

  // check both possible angles with the children
  // Note that if the second child passes the query
  // then this child is chosen independent of the
  // result of the first one
  // This is valid in our case since a division
  // is normally taking place with an angle of nearly 180
  // degrees between both division direction vectors
  for( int i=0;i<2;i++ )
  {
    double angle = acos( childrenDir * pNormal );
    angle = angle * 180./M_PI;

    if( angle < _changeLayerThreshold )
      childChoice = i+1;

    // after first traversal of loop invert the
    // the direction vector to compute the other angle
    childrenDir *= -1.;
  }

  return childChoice;
}

// ---------------------------------------------------------------------

void SCellLayers::updateCellHistory( cellHistoryMap &cellHistories,
                                     const SLineageTree *cell,
                                     const unsigned int newLayer )
{
  std::map<const SLineageTree*, cellHistory>::const_iterator iter;
  iter = _cellHistories->back().find( cell->parent );

  // search the cell in the cell history of the
  // previous time step
  if( iter != _cellHistories->back().end() )
  {
    // if cell found then copy its cell history
    // to the new one
    cellHistory ch = iter->second;
    // and add the new layer value if it is not already
    // stored in the set
    ch.insert( newLayer );

    cellHistories.insert( std::make_pair( cell, ch ) );
  }
}

// ---------------------------------------------------------------------

void SCellLayers::assignNewLayer( cellHistoryMap &cellHistories,
                                  const SLineageTree *cell,
                                  const int childChoice,
                                  const int prevLayer )
{
  if( prevLayer + 1 > _currentNumLayers )
  {
    _currentNumLayers++;
    _cellLayersPerTimeStep->at(cell->timeStep).resize( _currentNumLayers );
  }

  // the other child stays in its previous layer
  int stayInLayerChild;
  // if left child is selected then set index of child to stay to 1
  if( childChoice == 1 ) stayInLayerChild = 1;
  // if right child is selected then set index of child to stay to 0
  else stayInLayerChild = 0;

  // insert children in their corresponding layers
  _cellLayersPerTimeStep->at(cell->timeStep).at(prevLayer).insert( cell->children[childChoice-1] );
  this->updateCellHistory( cellHistories, cell->children[childChoice-1], prevLayer );
  _cellLayersPerTimeStep->at(cell->timeStep).at(prevLayer-1).insert( cell->children[stayInLayerChild] );
  this->updateCellHistory( cellHistories, cell->children[stayInLayerChild], prevLayer-1 );
}

// ---------------------------------------------------------------------

void SCellLayers::checkLayerSizes( const unsigned int timeStep )
{
  for( std::vector<cellSet>::iterator iter = _cellLayersPerTimeStep->at(timeStep).begin();
       iter != _cellLayersPerTimeStep->at(timeStep).end(); )
  {
    // if no cells are in this layer then erase it
    if( iter->size() == 0 )
    {
      iter = _cellLayersPerTimeStep->at(timeStep).erase( iter );
      _currentNumLayers--;
    }
    else
      ++iter;
  }
}

// ---------------------------------------------------------------------

osg::Vec3 SCellLayers::getCenterNormal( const osg::Vec3 &vPos )
{
  SArray center = _lineages->getCenter();
  osg::Vec3 c( center[0], center[1], center[2] );

  c = c * _inverseMatrix.inverse( _inverseMatrix );
  osg::Vec3 centerNormal( vPos[0] - c[0],
                          vPos[1] - c[1],
                          vPos[2] - c[2] );

  centerNormal.normalize();

  return centerNormal;
}

// ---------------------------------------------------------------------

unsigned int SCellLayers::determineLayerValueFromPair( const std::pair<int, int> &IDandT )
{
  unsigned int layerValue = 0;

  cellPairVector layers;
  // first convert SLineageTree pointers to pairs for each node
  BOOST_FOREACH( cellSet cSIL, _cellLayersPerTimeStep->at(IDandT.second-1) )
  {
    cellPairSet pairCellSet;

    BOOST_FOREACH( const SLineageTree* tree, cSIL )
    {
      pairCellSet.insert( std::make_pair( tree->cellId, tree->timeStep ) );
    }

    layers.push_back( pairCellSet );
  }

  // then determine layer value
  // loop for each layer
  BOOST_FOREACH( cellPairSet cSIL, layers )
  {
    cellPairSet::const_iterator iter = cSIL.find( IDandT );
    if( iter == cSIL.end() )
      layerValue++;
    else
      break;
  }

  return layerValue;
}

// ---------------------------------------------------------------------

unsigned int SCellLayers::determineLayerValue( const SLineageTree *cell )
{
  unsigned int layerValue = 0;

  // for each layer
  BOOST_FOREACH( cellSet cSIL, _cellLayersPerTimeStep->at(cell->timeStep-1) )
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

osg::Vec3 SCellLayers::getNearestNeighborNormal( const std::vector<boost::shared_ptr<SConvexHull> > &currentDTVector,
                                                 const unsigned int layerValue,
                                                 const unsigned int timeStep,
                                                 const osg::Vec3 &cellPos )
{
  // this method is called for inner cells within
  // the delaunay triangulation
  // the vertex normal is given by the nearest outer
  // cell in the same layer
  double distance = std::numeric_limits<double>::max();
  // store vertex normal of nearest cell
  osg::Vec3 nearestNormal;

  // first search in the current layer of the cell
  // for each layer
  BOOST_FOREACH( cellSet::value_type cellInLayer, _cellLayersPerTimeStep->at(timeStep).at(layerValue) )
  {
    osg::Vec3 cPos( cellInLayer->getX(), cellInLayer->getY(), cellInLayer->getZ() );
    osg::Vec3 vNormal = currentDTVector[layerValue]->getVertexNormal( cPos );

    // if the current checked cell is not an inner cell
    if( vNormal != osg::Vec3( 0., 0., 0. ) )
    {
      if( (cPos - cellPos).length() < distance )
      {
        distance = (cPos - cellPos).length();
        nearestNormal = vNormal;
      }
    }
  }

  // if a nearest cell in the same layer was found
  if( distance != std::numeric_limits<double>::max() )
    return nearestNormal;
  // if the layer only consists of at most cell then select the
  // nearest cell in another layer
  else
  {
    unsigned int layerCount = 0;

    // for each layer
    BOOST_FOREACH( cellSet layers, _cellLayersPerTimeStep->at(timeStep) )
    {
      // skip the previous layer
      if( layerCount == layerValue )
        continue;

      // for each cell in the current layer
      BOOST_FOREACH( cellSet::value_type cellInLayer, layers )
      {
        osg::Vec3 cPos( cellInLayer->getX(), cellInLayer->getY(), cellInLayer->getZ() );
        osg::Vec3 vNormal = currentDTVector[layerCount]->getVertexNormal( cPos );

        // if the current checked cell is not an inner cell
        if( vNormal != osg::Vec3( 0., 0., 0. ) )
        {
          if( (cPos - cellPos).length() < distance )
          {
            distance = (cPos - cellPos).length();
            nearestNormal = vNormal;
          }
        }
      }

      layerCount++;
    }

    return nearestNormal;
  }
}

// ---------------------------------------------------------------------

osg::Vec3 SCellLayers::getNearestNeighborNormal( const std::vector<boost::shared_ptr<SAlphaShape> > &currentDTVector,
                                                 const unsigned int layerValue,
                                                 const unsigned int timeStep,
                                                 const osg::Vec3 &cellPos )
{
  // this method is called for inner cells within
  // the delaunay triangulation
  // the vertex normal is given by the nearest outer
  // cell in the same layer
  double distance = std::numeric_limits<double>::max();
  // store vertex normal of nearest cell
  osg::Vec3 nearestNormal;

  // first search in the current layer of the cell
  // for each layer
  BOOST_FOREACH( cellSet::value_type cellInLayer, _cellLayersPerTimeStep->at(timeStep).at(layerValue) )
  {
    osg::Vec3 cPos( cellInLayer->getX(), cellInLayer->getY(), cellInLayer->getZ() );
    osg::Vec3 vNormal = currentDTVector[layerValue]->getVertexNormal( cPos );

    // if the current checked cell is not an inner cell
    if( vNormal != osg::Vec3( 0., 0., 0. ) )
    {
      if( (cPos - cellPos).length() < distance )
      {
        distance = (cPos - cellPos).length();
        nearestNormal = vNormal;
      }
    }
  }

  // if a nearest cell in the same layer was found
  if( distance != std::numeric_limits<double>::max() )
    return nearestNormal;
  // if the layer only consists of at most cell then select the
  // nearest cell in another layer
  else
  {
    unsigned int layerCount = 0;

    // for each layer
    BOOST_FOREACH( cellSet layers, _cellLayersPerTimeStep->at(timeStep) )
    {
      // skip the previous layer
      if( layerCount == layerValue )
        continue;

      // for each cell in the current layer
      BOOST_FOREACH( cellSet::value_type cellInLayer, layers )
      {
        osg::Vec3 cPos( cellInLayer->getX(), cellInLayer->getY(), cellInLayer->getZ() );
        osg::Vec3 vNormal = currentDTVector[layerCount]->getVertexNormal( cPos );

        // if the current checked cell is not an inner cell
        if( vNormal != osg::Vec3( 0., 0., 0. ) )
        {
          if( (cPos - cellPos).length() < distance )
          {
            distance = (cPos - cellPos).length();
            nearestNormal = vNormal;
          }
        }
      }

      layerCount++;
    }

    return nearestNormal;
  }
}

// ---------------------------------------------------------------------

void SCellLayers::transformLayerData()
{
  _cellsPerTimeStepInLayer = boost::shared_ptr<NodeFeatureInfo>( new NodeFeatureInfo );
  _cellsPerTimeStepInLayer->resize( _cellsPerTimeStep->size() );

  // for each time step
  for( std::size_t t = 0; t < _cellsPerTimeStepInLayer->size(); t++ )
  {
    // for each layer
    int layer = 0;
    for( std::vector<cellSet>::iterator lIter = _cellLayersPerTimeStep->at(t).begin();
         lIter != _cellLayersPerTimeStep->at(t).end(); lIter++, layer++ )
    {
      for( cellSet::iterator cIter = lIter->begin(); cIter != lIter->end(); cIter++ )
        _cellsPerTimeStepInLayer->at(t).insert( std::make_pair( (*cIter)->cellId, layer ) );
    }
  }
}

// ---------------------------------------------------------------------
