#include "SCellTreeVisualizer.hh"

#include "../similarityMeasure/SLineageTreeLayers.hh"
#include "../similarityMeasure/STreeDivisionSequence.hh"
#include "../similarityMeasure/SAveragedTreeDivisionSequence.hh"
#include "../similarityMeasure/SLineageCallback.hh"
#include "../similarityMeasure/SSimilarityMeasureUtil.hh"

#include <boost/lexical_cast.hpp>

// ---------------------------------------------------------------------

SCellTreeVisualizer::SCellTreeVisualizer( const std::vector< boost::shared_ptr<SAdvancedLineageTreeCollection> > &lineages,
                                          const std::vector<std::size_t> &numTotalTimesteps,
                                          const std::vector<NodeFeatureInfo> &divisionScheme,
                                          const std::vector<NodeFeatureInfo> &layerValues,
                                          const std::vector< std::map<int,int> > &cellFiles,
                                          const std::vector<osg::Matrix> &rotMatrices,
                                          const std::vector< std::vector<osg::Vec3> > &centresOfMass,
                                          const osg::ref_ptr<osg::Vec4Array> layerColors,
                                          const bool onlyMasterFiles,
                                          const int renderLineageLineType,
                                          const int renderLineageNodeType )
  : _lineages( lineages ),
    _numTotalTimesteps( numTotalTimesteps ),
    _divisionScheme( divisionScheme ),
    _layerValues( layerValues ),
    _cellFiles( cellFiles ),
    _rotMatrices( rotMatrices ),
    _centresOfMass( centresOfMass ),
    _layerColors( layerColors ),
    _onlyMasterFiles( onlyMasterFiles ),
    _renderLineageLineType( renderLineageLineType ),
    _renderLineageNodeType( renderLineageNodeType )
{
}

// ---------------------------------------------------------------------

void SCellTreeVisualizer::generateLineageTrees( const bool registerTrees,
                                                const std::size_t divStop,
                                                osg::ref_ptr<osg::Group> addToThisGroup,
                                                const SNodeInfo &info,
                                                const int windowId )
{
  double xSpacing = 4.;
  double offset = 0.;
  double maxTime = 0.;
  double xMax;

  // manual set time steps for registered number of cells in [18,143]
  std::vector< std::pair<std::size_t,std::size_t> > registeredTimeSteps;
  if( registerTrees )
  {
    registeredTimeSteps.push_back( std::make_pair( 36, 269 ) );
    registeredTimeSteps.push_back( std::make_pair( 2, 277 ) );
    registeredTimeSteps.push_back( std::make_pair( 1, 230 ) );
    registeredTimeSteps.push_back( std::make_pair( 73, 344 ) );
    registeredTimeSteps.push_back( std::make_pair( 2, 213 ) );
  }

  addToThisGroup = new osg::Group;
  addToThisGroup->setName( "Lineage Trees" );

  // TODO
//  SInteractiveCellHighlightingCallback *ichc = new SInteractiveCellHighlightingCallback();
//  _connections.push_back( this->onHighlightChanged( boost::bind( &SInteractiveCellHighlightingCallback::updateCellInfo,
//                                                                 ichc, _1 ) ) );

//  SInteractiveCellHighlighting *ich = new SInteractiveCellHighlighting();
//  ich->setName( "Selection" );
//  ich->setUpdateCallback( ichc );

  // loop over all data sets
  for( std::size_t d = 0; d < _lineages.size(); d++ )
  {
    int time = _numTotalTimesteps.at(d);

    if( maxTime < time )
      maxTime = time;

    int treeCounter = 0;

    std::size_t firstRadialDivision = time;

    osg::ref_ptr<osg::Group> lineageDataGroup = new osg::Group;
    lineageDataGroup->setName( "Lineages of data " + boost::lexical_cast<std::string>( d ) );

    std::size_t min, max;
    std::map<std::size_t, std::size_t> lifeDuration;
    SSimilarityMeasureUtil::determineCellLifeDuration( _lineages.at(d), _centresOfMass.at(d),
                                                       _cellFiles.at(d), _rotMatrices.at(d),
                                                       lifeDuration, min, max );
    SColorMap lifeDurationColorMap;
    lifeDurationColorMap.setCoolToWarm();
    lifeDurationColorMap.scaleColors( min, max );
    lifeDurationColorMap.setNumberColorbarTicks( 8 );

    for( SAdvancedLineageTreeCollection::const_iterator l = _lineages.at(d)->begin();
         l != _lineages.at(d)->end(); ++l )
    {
      SKernel::get()->setProgress( (double)treeCounter/(double)(_lineages.at(d)->size() - 1),
                                   _lineages.at(d)->getName() + " : Lin. gen." );

      // get lineage id
      int id = l->first;

      // we create a new copy of this lineage tree in order
      // to change its layout which should not happen
      // to the original data since for them we store the
      // pointers to identify them
      SLineageTree *tree = new SLineageTree( *(*_lineages.at(d))[id] );

      // get cell file value of root node of the current lineage
      // which is the same for the whole tree and all nodes of it
      int cellFile = _cellFiles.at(d).at(id);

      // if only the master cell files should be considered
      // then ignore all lineages with of different cell files
      if( _onlyMasterFiles && cellFile != 0 )
        continue;

      osg::ref_ptr<osg::Group> treeGroup = new osg::Group;
      treeGroup->setName( "Lineage " + boost::lexical_cast<std::string>( id ) );
      treeGroup->setDataVariance( osg::Object::DYNAMIC );

      osg::ref_ptr<osg::Group> sequenceGroup = new osg::Group;

      // cut all sub trees below the divStop-th division
      for( SLineageTree::iterator treeIter = tree->begin();
           treeIter != tree->end(); ++treeIter )
      {
        if( treeIter->getCellCycleId() == divStop + 1 )
        {
          for( std::size_t c = 0; c < treeIter->children.size(); c++ )
            treeIter->children[c]->cropThisSubtree();
        }
      }

      // compute the layout
      tree->computeLayout( xSpacing );

      // properties setting
      osg::ref_ptr<SLineageTreeLayers> lineageTree =
          new SLineageTreeLayers( tree,
                                  lifeDuration,
                                  lifeDurationColorMap,
                                  _layerColors,
                                  _layerValues.at(d),
                                  _divisionScheme.at(d),
                                  _renderLineageLineType,
                                  _renderLineageNodeType );

      // after computing the layout, generate the division sequence and render them
      for( SLineageTree::iterator treeIter = tree->begin();
           treeIter != tree->end(); ++treeIter )
      {
        if( registerTrees )
        {
          if( treeIter->timeStep == registeredTimeSteps.at(d).second )
          {
            SSimilarityMeasureGraphics::renderDivisionSequence( *treeIter, _divisionScheme.at(d),
                                                                registeredTimeSteps.at(d).first, sequenceGroup );
          }
        }
        else
        {
          // for each leaf of the tree
          if( treeIter->children.size() == 0 )
          {
            SSimilarityMeasureGraphics::renderDivisionSequence( *treeIter, _divisionScheme.at(d),
                                                                1, sequenceGroup );
          }
        }

        // check when the first radial division occurred
        if( treeIter->children.size() == 2 )
        {
          // get division type of current node
          auto iter = _divisionScheme.at(d).at(treeIter->timeStep-1).find( treeIter->cellId );
          if( iter != _divisionScheme.at(d).at(treeIter->timeStep-1).end() )
          {
            DivisionType::DivisionType divType = (DivisionType::DivisionType)iter->second;

            if( divType == DivisionType::RADIAL )
            {
              if( firstRadialDivision > treeIter->timeStep-1 )
                firstRadialDivision = treeIter->timeStep-1;
            }
          }
        }
      }

      lineageTree->setDivisionSize( 10. );
      lineageTree->setRootSize( 5. );
      lineageTree->setLeafSize( 5. );
      lineageTree->setLineWidth( 5. );
      lineageTree->renderNoMoveNodes( true );
      lineageTree->setMoveSize( 2. );
      lineageTree->render();

      // TODO
//      SLineageTree::const_iterator cIter = tree->begin();

//      for( SLineageTreeLayers::iterator cell = lineageTree->begin();
//           cell != lineageTree->end(); ++cell, ++cIter )
//      {
//        // since we consider multiple data sets, we also have to include
//        // the data id to get an unique assignment to trees
//        cellDataInfo cDI;
//        cDI.push_back( cIter->cellId );
//        cDI.push_back( cIter->timeStep );
//        cDI.push_back( d );

//        if( !(cIter->parent) || cIter->children.size() != 1 || _renderMoveNodes )
//        {
//          SActions::signal_node_t signal2 = SActions::get( **cell ).getNodeSignal( "Highlight cell" );
//          _connections.push_back( signal2->connect(
//                                    boost::bind( &SVLRBrowser::highlightCell, this, _1, NodeType::NODE2D ) ) );

//          // TODO: set max layers and refactor assignLayer method
//          if( _draw3DCells && _lineageColorType == 0 && false )
//          {
//            for( std::size_t i=0; i < _cellLayers->getMaxNumLayers();i++ )
//            {
//              SActions::signal_node_t signal1 = SActions::get( **cell ).getNodeSignal(
//                    "Assign layer " + boost::lexical_cast<std::string>(i) );
//              _connections.push_back( signal1->connect(
//                                        boost::bind( &SVLRBrowser::assignLayer, this, _1, i, NodeType::NODE2D ) ) );
//            }
//          }

//          // add entry for setting division type info only if this current cell is a dividing cell
//          if( cIter->children.size() == 2 )
//          {
//            // only two allowed division types: anticlinal and radial
//            // that do not change the layering
//            SActions::signal_node_t signal3 = SActions::get( **cell ).getNodeSignal(
//                  "Set to anticlinal division" );
//            _connections.push_back( signal3->connect(
//                                      boost::bind( &SVLRBrowser::setDivisionType, this, _1, 0, NodeType::NODE2D ) ) );

//            SActions::signal_node_t signal4 = SActions::get( **cell ).getNodeSignal(
//                  "Set to radial division" );
//            _connections.push_back( signal4->connect(
//                                      boost::bind( &SVLRBrowser::setDivisionType, this, _1, 2, NodeType::NODE2D ) ) );
//          }

//          _geomNodeMaps->insertCellInfoGeometry( std::make_pair( cDI[0], cDI[1] ), *cell );
//        }

//        // insert node position for highlighting arrows
//        osg::Vec3 pos( cIter->getX() + offset, -cIter->getY(), 5. );
//        ich->addPosition( cDI, pos );
//      }

      /* get lineage tree geometry */

      osg::ref_ptr<osg::MatrixTransform> trafo;
      osg::Vec3d translate( offset, 0, 0 );
      trafo = new osg::MatrixTransform( osg::Matrix::translate( translate ) );
      trafo->setName( "Lineage" );
      trafo->addChild( lineageTree->getRenderedLineage() );
      treeGroup->addChild( trafo );

      /* get lineage tree labels */

      //osg::ref_ptr<osg::MatrixTransform> m = new osg::MatrixTransform( osg::Matrix::translate( offset + tree->getX(), 10, 0 ) );
      osg::ref_ptr<osg::MatrixTransform> m = new osg::MatrixTransform( osg::Matrix::rotate( M_PI/4, osg::Vec3(0,0,1) ) );
      m->postMult( osg::Matrix::translate( offset + tree->getX(), 23, 0 ) );
      //m->postMult( osg::Matrix::translate( offset + tree->getX(), 5, 0 ) );
      //m->postMult( osg::Matrix::translate( offset + tree->getX(), 40, 0 ) );
      m->addChild( lineageTree->getLineageFileLabels( d, cellFile ) );
      m->setName( "Label" );
      treeGroup->addChild( m );

      /* generate division sequence visualization */

      osg::ref_ptr<osg::MatrixTransform> mS = new osg::MatrixTransform( osg::Matrix::translate( offset, -maxTime-10., 0. ) );
      //osg::MatrixTransform *mS = new osg::MatrixTransform( osg::Matrix::translate( offset, 30., 0. ) );
      mS->addChild( sequenceGroup );
      mS->setName( "Division Sequence" );
      treeGroup->addChild( mS );

      /* set callback for color update */

      osg::ref_ptr<SLineageCallback> lineageCallback = new SLineageCallback();
      lineageTree->addUpdateCallback( lineageCallback );
      treeGroup->addChild( lineageTree );
      lineageCallback->updateColors();

      lineageDataGroup->addChild( treeGroup );
      offset += tree->width() + 8. * xSpacing;
      // update current maximum x value
      // for adequate drawing of time line labels
      xMax = offset;

      treeCounter++;
    }

    // add lineage set of the current data to the lineage group
    addToThisGroup->addChild( lineageDataGroup );

    // print the first occurrence of a radial division
    std::cout << "First radial division: " << firstRadialDivision << std::endl;
  }

  osg::ref_ptr<osg::StateSet> stateset = addToThisGroup->getOrCreateStateSet();
  stateset->setMode( GL_LIGHTING, osg::StateAttribute::OFF);

  // at last render time line labels
  osg::ref_ptr<osg::Geode> labels = new osg::Geode;
  labels->setName( "Time Line Labels" );
  SSimilarityMeasureGraphics::drawTimeLineLabels( labels, 25., maxTime-1., -30., xMax );
  addToThisGroup->addChild( labels );

  // at last add selection callback which draws the four arrows
  // highlighting the current selected node
//  addToThisGroup->addChild( ich );

  SKernel::get()->render( addToThisGroup, info, windowId );
}

// ---------------------------------------------------------------------

void SCellTreeVisualizer::generateTreeDivisionSequence( const std::size_t sequenceAxisType,
                                                        const bool loadModel,
                                                        const SNodeInfo &info,
                                                        const int windowId )
{
  osg::ref_ptr<osg::Group> divisionGroup = new osg::Group;
  divisionGroup->setName( "Tree Division Sequences" );
  double curHeight = 0.;
  double curWidth = 0.;
  double spacing = 25.;
  std::size_t dd = _lineages.size();
  //if( _numData == 10 )
  //   dd = 5;

  bool showAverage = false;

  if( showAverage )
  {
    STreeDivisionSequence seq( "/home/necrolyte/Uni/LateralRootGrowth/230615/averagedLayeringRealData.csv" );
    seq.renderTreeDivisionSequence();
    osg::ref_ptr<osg::MatrixTransform> treeMat = seq.getTree();
    treeMat->postMult( osg::Matrix::translate( 0., curHeight, 0. ) );
    divisionGroup->addChild( treeMat );
  }
  else
  {
    // generate sequences for the 5 VLR data sets
    for( std::size_t d = 0; d < dd; d++ )
    {
      osg::MatrixTransform *rotMatrix = new osg::MatrixTransform( osg::Matrix::identity() );
      // then apply the data rotations
      SArray rotation = _lineages.at(d)->getRotation();
      rotMatrix->postMult( osg::Matrix::rotate( rotation[2]/180.*M_PI, osg::Vec3(0,0,1) ) );
      rotMatrix->postMult( osg::Matrix::rotate( rotation[1]/180.*M_PI, osg::Vec3(0,1,0) ) );
      rotMatrix->postMult( osg::Matrix::rotate( rotation[0]/180.*M_PI, osg::Vec3(1,0,0) ) );

      AXISTYPE type;
      if( sequenceAxisType == 0 )
        type = AXISTYPE::TIMESTEP;
      else
        type = AXISTYPE::CELLS;

      STreeDivisionSequence seq( _lineages.at(d), _onlyMasterFiles,
                                 rotMatrix->getMatrix(), loadModel, type );
      seq.renderTreeDivisionSequence();
      seq.setDataName( _lineages.at(d)->getName(), 3. );
      osg::ref_ptr<osg::MatrixTransform> treeMat = seq.getTree();
      treeMat->postMult( osg::Matrix::translate( 0., curHeight, 0. ) );
      treeMat->setName( "Data " + _lineages.at(d)->getName() );
      divisionGroup->addChild( treeMat );
      curHeight += seq.getTreeWidth() + spacing;

      if( curWidth < seq.getTreeHeight() )
        curWidth = seq.getTreeHeight();
    }
  }

  /*
  curHeight = 0.;

  // generate sequences for the 5 model data sets
  for( std::size_t d = dd; d < _numData; d++ )
  {
    STreeDivisionSequence seq( _lineages.at(d), _onlyMasterFiles,
                               _vlrRoot->getMatrix(), _loadModel,
                               true, AXISTYPE::TIMESTEP );
    seq.renderTreeDivisionSequence();
    seq.setDataName( _lineages.at(d)->getName(), 3. );
    osg::ref_ptr<osg::MatrixTransform> treeMat = seq.getTree();
    treeMat->postMult( osg::Matrix::translate( curWidth+spacing, curHeight, 0. ) );
    treeMat->setName( "Data " + _lineages.at(d)->getName() );
    divisionGroup->addChild( treeMat );
    curHeight += seq.getTreeWidth() + spacing;
  }
  */

  SSVGRenderer svg( true );
  svg.startImage( "/tmp/treeDivisionSequence.svg" );
  svg.renderOSGNode( divisionGroup, 1., 1., 0.25 );

  osg::StateSet* stateset = divisionGroup->getOrCreateStateSet();
  stateset->setMode( GL_LIGHTING, osg::StateAttribute::OFF);

  SKernel::get()->render( divisionGroup, info, windowId );
}

// ---------------------------------------------------------------------

void SCellTreeVisualizer::generateAveragedTreeDivisionSequence( const std::string &modelDataDirectory,
                                                                const SNodeInfo &info,
                                                                const int windowId )
{
  osg::ref_ptr<osg::Group> divisionGroup = new osg::Group;
  divisionGroup->setName( "Averaged Tree Division Sequences" );

  SAveragedTreeDivisionSequence atds( modelDataDirectory );
  atds.renderTreeDivisionSequence();
  //atds.setDataName( _lineages.at(d)->getName(), 3. );
  osg::ref_ptr<osg::MatrixTransform> treeMat = atds.getTree();
  divisionGroup->addChild( treeMat );

  SSVGRenderer svg( true );
  svg.startImage( "/tmp/averagedTreeDivisionSequence.svg" );
  svg.renderOSGNode( divisionGroup, 1., 1., 0.25 );

  osg::StateSet* stateset = divisionGroup->getOrCreateStateSet();
  stateset->setMode( GL_LIGHTING, osg::StateAttribute::OFF);

  SKernel::get()->render( divisionGroup, info, windowId );
}

// ---------------------------------------------------------------------
