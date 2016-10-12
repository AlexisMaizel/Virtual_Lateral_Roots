#include "SVLRBrowser.hh"

/**
  @file   SVLRBrowser.cc
  @brief  Contains class for creating a feature browser of VLR data
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "../similarityMeasure/SInteractiveCellHighlighting.hh"
#include "../similarityMeasure/SInteractiveTexture.hh"
#include "../similarityMeasure/SInteractiveDivisionArrows.hh"
#include "../similarityMeasure/SHighlightCellVisitor.hh"
#include "../similarityMeasure/SColorCellVisitor.hh"
#include "../similarityMeasure/SLayerInformationIO.hh"
#include "../similarityMeasure/SGraphDescriptorIO.hh"
#include "../similarityMeasure/SArrow3DCylindrical.hh"
#include "../similarityMeasure/SArrow3DAngled.hh"
#include "../similarityMeasure/SCellAnalysis.hh"
#include "../similarityMeasure/SCurvatureAnalysis.hh"
#include "../similarityMeasure/SDevelopmentAnalysis.hh"
#include "../similarityMeasure/SDivisionAnalysis.hh"
#include "../similarityMeasure/SScalarGrid.hh"
#include "../similarityMeasure/SSimilarityMeasureGraphics.hh"
#include "../similarityMeasure/SVLRDataAnalysisIO.hh"
#include "../similarityMeasure/SCellTreeVisualizer.hh"
#include "../similarityMeasure/SCorrelationAnalysis.hh"

#include "common/VException.hh"

#include "math/SPrincipalComponentAnalysis.hh"

#include "graphics/SLabelGenerator.hh"
#include "graphics/SActions.hh"
#include "graphics/SSelectableSphere.hh"
#include "graphics/SColorMap.hh"
#include "graphics/SScreenshot.hh"
#include "graphics/SColorMapTF.hh"

#include "kernel/SAlgorithmProfile.hh"
#include "kernel/SKernel.hh"
#include "kernel/SWindowData.hh"

#include "dataSet/SDataManager.hh"

#include <string>

#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/conversion/bounds.hpp>

#include <osg/Geode>
#include <osg/Texture2D>
#include <osg/Image>
#include <osg/PositionAttitudeTransform>
#include <osgDB/ReadFile>
#include <osgText/Text>
#include <osg/io_utils>

// ---------------------------------------------------------------------

typedef SAlgorithm* maker_t();
extern std::map<std::string, maker_t*> algorithmFactory;

// ---------------------------------------------------------------------

extern "C" 
{
SAlgorithm* maker()
{
  return( new SVLRBrowser() );
}

class proxy
{
public:
  proxy()
  {
    algorithmFactory[ SVLRBrowser::algoName() ] = maker;
  }
};

proxy p;
}

// ---------------------------------------------------------------------

SVLRBrowser::SVLRBrowser()
  : _lineageWindowId( SKernel::get()->requestNewWindow(
                      SWindowData( algoName() + " - Lineages ", this, SWindowData::SINGLE2D ) ) ),
    _divisionAnalysisWindowId( SKernel::get()->requestNewWindow(
                                SWindowData( algoName() + " - Division Analysis ", this, SWindowData::SINGLE2D,
                                             osg::Vec4( 1., 1., 1., 1. ), true ) ) ),
    _lineageInfo( SNodeInfo( algoName() + " - Lineages", this ) ),
    _divisionWindowId( SKernel::get()->requestNewWindow(
                          SWindowData( algoName() + " - Divisions ", this, SWindowData::SINGLE2D ) ) ),
    _divisionInfo( SNodeInfo( algoName() + " - Divisions", this ) ),
    _averagedDivisionWindowId( SKernel::get()->requestNewWindow(
                              SWindowData( algoName() + " - Averaged Divisions ", this, SWindowData::SINGLE2D ) ) ),
    _averagedDivisionInfo( SNodeInfo( algoName() + " - Averaged Divisions", this ) ),
    _correlationAnalysisWindowId( SKernel::get()->requestNewWindow(
                              SWindowData( algoName() + " - Correlation Analysis ", this, SWindowData::SINGLE2D ) ) ),
    _correlationAnalysisInfo( SNodeInfo( algoName() + " - Correlation Analysis", this ) ),
    _numData( 10 ),
    _tfoc( 1 ),
    _twin( 20 ),
    _firstClick( true ),
    _selectionColor( osg::Vec4( 49./255., 163./255., 84./255., 1. ) ),
    _cellWallWidth( 0.005 )
    //_selectionColor( osg::Vec4( 190./255., 190./255., 190./255., 1. ) )
{
  _vlrWindowId.resize( _numData );
  _vlrInfo.resize( _numData );
  //_vlrRoot.resize( _numData );

  _sliderLabels.push_back( "Time step" );
  _sliderLabels.push_back( "Cell Cycle" );
  _sliderLabels.push_back( "Registered Step" );

  for( std::size_t d = 0; d < _numData; d++ )
  {
    std::string dataName = boost::lexical_cast<std::string>( d );
    _vlrWindowId.at(d) = SKernel::get()->requestNewWindow(
          SWindowData( algoName() + " " + dataName, this,
                       SWindowData::COMPOSITE3D, osg::Vec4( 1., 1., 1., 1. ), true ) );

    _vlrInfo.at(d) = SNodeInfo( algoName() + " " + dataName, this );
  }

  // initialization of cluster window
  _layerWindowCreator = new SCellLayerWindowCreator( "VLR Properties" );

  // registered datatypes for qt
  qRegisterMetaType<cellHistory>( "cellHistory" );
  qRegisterMetaType<std::size_t>( "std::size_t" );

  // signals are not connected at the beginning
  _signalsConnected = false;
}

// ---------------------------------------------------------------------

SVLRBrowser::~SVLRBrowser()
{
  // delete window creator instance
  if( _layerWindowCreator )
    delete( _layerWindowCreator );

  SKernel::get()->removeNodeInfo( _vlrInfo.at(0) );
  SKernel::get()->removeNodeInfo( _lineageInfo );
  SKernel::get()->removeNodeInfo( _divisionInfo );
  SKernel::get()->removeNodeInfo( _averagedDivisionInfo );
  SKernel::get()->removeNodeInfo( _correlationAnalysisInfo );
  SKernel::get()->removeNodeInfo( _colorbarInfo );
}

// ---------------------------------------------------------------------

void SVLRBrowser::initProfile()
{
  _numData = 6;
  _profile->addEntry<std::size_t>( _numData, "Number of lineage data" );

  _lineageDataId.resize( _numData );
  for( std::size_t d = 0; d < _numData; d++ )
  {
    std::string dataStr = boost::lexical_cast<std::string>( d );
    _profile->addDatasetSelector( _lineageDataId[d],
                                  "Lineage Data " + dataStr, _lineageFilter );
  }

  _sDataId = 0;
  _profile->addEntry<unsigned int>( _sDataId, "Data ID used for cell visualization" );

  _loadLayerInformation = true;
  _profile->addCheckBox( _loadLayerInformation, "Load layer info from data set" );

  _loadDivisionTypeInformation = true;
  _profile->addCheckBox( _loadDivisionTypeInformation, "Load division type info from data set" );

  _sliderType = 0;
  _sliderTypes.clear();
  _sliderTypes.push_back( "Time Step" );
  _sliderTypes.push_back( "Cell Cycle" );
  _sliderTypes.push_back( "Registered Steps" );
  _profile->addComboBox( _sliderType, _sliderTypes, "Slider type" );

  _onlyMasterFiles = false;
  _profile->addCheckBox( _onlyMasterFiles, "Consider only master cell file" );

  _divisionLabels = false;
  //_profile->addCheckBox( _divisionLabels, "Show division labels" );

  _profile->addTab( "3D Cell/Surface Visualization" );

  _draw3DCells = true;
  _profile->addCheckBox( _draw3DCells, "Render 3D cells and/or layer surfaces" );

  _cellRadius = 6.;
  _profile->addEntry<double>( _cellRadius, "Cell radius" );

  // render type of 3D VLR
  // 0 -> Layering: spheres are colored based on layers plus delaunay triangulation
  // of each layer
  // 1 -> Cell Files: spheres are colored based on cell file assignment plus delaunay
  // triangulation of whole primordium
  _render3DType = 3;
  _render3DTypes.clear();
  _render3DTypes.push_back( "Layers" );
  _render3DTypes.push_back( "Cell Files" );
  _render3DTypes.push_back( "Lineage Trees" );
  _render3DTypes.push_back( "Pair Layers" );
  _render3DTypes.push_back( "Single Color" );
  _render3DTypes.push_back( "Life Duration" );
  _render3DTypes.push_back( "Quiescent Center" );
  _profile->addComboBox( _render3DType, _render3DTypes, "Sphere coloring type" );

  _divisionColorType = 0;
  _divisionColorTypes.clear();
  _divisionColorTypes.push_back( "Type" );
  _divisionColorTypes.push_back( "Layer" );
  _profile->addComboBox( _divisionColorType, _divisionColorTypes, "Division coloring type" );

  _renderArrows = false;
  _profile->addCheckBox( _renderArrows, "Render division arrows" );

  _renderTracks = false;
  _profile->addCheckBox( _renderTracks, "Render 3D lineage tracks (only cycle slider)" );

  _renderCompletePrimordium = true;
  _profile->addCheckBox( _renderCompletePrimordium, "Triangulation of complete primordium" );

  _shapeType = 1;
  _shapeTypes.clear();
  _shapeTypes.push_back( "Convex Hull" );
  _shapeTypes.push_back( "Alpha Shape" );
  _profile->addComboBox( _shapeType, _shapeTypes, "Shape type" );

  _renderCurvature = false;
  _profile->addCheckBox( _renderCurvature, "Render curvature information (only time slider)" );

  _curvatureType = 0;
  _curvatureTypes.clear();
  _curvatureTypes.push_back( "Mean" );
  _curvatureTypes.push_back( "Gaussian" );
  _profile->addComboBox( _curvatureType, _curvatureTypes, "Curvature type" );

  _wireframe = false;
  _profile->addCheckBox( _wireframe, "Wireframe mode of surface" );

  _profile->addTab( "Lineage Visualization" );

  _drawLineages = true;
  _profile->addCheckBox( _drawLineages, "Render 2D lineages" );

  _renderLineageLineType = 1;
  _renderLineageLineTypes.clear();
  _renderLineageLineTypes.push_back( "Single Color" );
  _renderLineageLineTypes.push_back( "Pair Layers" );
  _renderLineageLineTypes.push_back( "Life Duration" );
  _profile->addComboBox( _renderLineageLineType, _renderLineageLineTypes, "Line coloring type" );

  _renderLineageNodeType = 3;
  _renderLineageNodeTypes.clear();
  _renderLineageNodeTypes.push_back( "Single Color" );
  _renderLineageNodeTypes.push_back( "Pair Layers" );
  _renderLineageNodeTypes.push_back( "Life Duration" );
  _renderLineageNodeTypes.push_back( "Division Type" );
  _profile->addComboBox( _renderLineageNodeType, _renderLineageNodeTypes, "Node coloring type" );

  _divStop = 10;
  _profile->addEntry<std::size_t>( _divStop, "Render trees until division sequence" );

  _registerTrees = false;
  _profile->addCheckBox( _registerTrees, "Register trees based on the number of cells" );

  _renderTreeDivisionSequence = false;
  _profile->addCheckBox( _renderTreeDivisionSequence, "Render periclinal division sequences" );

  _sequenceAxisType = 0;
  _sequenceAxisTypes.clear();
  _sequenceAxisTypes.push_back( "Time steps" );
  _sequenceAxisTypes.push_back( "Cells" );
  _profile->addComboBox( _sequenceAxisType, _sequenceAxisTypes, "Axis type for Sequences" );

//  _lineageColorType = 0;
//  _lineageColorTypes.clear();
//  _lineageColorTypes.push_back( "Layer/DivisionType" );
//  _lineageColorTypes.push_back( "Area size" );
//  _profile->addComboBox( _lineageColorType, _lineageColorTypes, "Lineage color type" );

  _profile->addTab( "Division Analysis" );

  _renderCorrelationAnalysis = false;
  _profile->addCheckBox( _renderCorrelationAnalysis, "Render correlation analysis" );

  _renderDivisionAnalysis = false;
  _profile->addCheckBox( _renderDivisionAnalysis, "Render division analysis" );

  _gridResolution = 15;
  _profile->addEntry<std::size_t>( _gridResolution, "Resolution of grid" );

  _angleTypeDivAnalysis = 0;
  _angleTypesDivAnalysis.clear();
  _angleTypesDivAnalysis.push_back( "Centre of mass" );
  _angleTypesDivAnalysis.push_back( "Primordium centre" );
  _angleTypesDivAnalysis.push_back( "Surface normal" );
  _angleTypesDivAnalysis.push_back( "Master Height Axis" );
  _profile->addComboBox( _angleTypeDivAnalysis, _angleTypesDivAnalysis, "Angle types for division analysis" );

  _regSteps = 20;
  _profile->addEntry<std::size_t>( _regSteps, "Number of registerd steps" );


//  _profile->addTab( "Cell Layers" );

  // threshold of a cell for changing its layer or not
  // in other words: anticlinal/radial or periclinal division
  _changeLayerThreshold = 50.;
//  _profile->addEntry<double>( _changeLayerThreshold, "Angle threshold: Anticlinal/Radial or Periclinal" );

  // threshold for checking if the current division is an anticlinal
  // or radial division
  _radialAngleThreshold = 45.;
//  _profile->addEntry<double>( _radialAngleThreshold, "Angle threshold: Anticlinal or Radial" );

  _profile->addTab( "External Data Analysis" );

  _enableSlicing = false;
  _profile->addCheckBox( _enableSlicing, "Enable geometry slicing" );

  _projDir = 3;
  _projections.clear();
  _projections.push_back("x->");
  _projections.push_back("y->");
  _projections.push_back("z->");
  _projections.push_back("all->");
  _profile->addComboBox( _projDir, _projections, "Direction of raw projection as well as slicing" );

  _MIPfilename = "";
  _profile->addDirectory( _MIPfilename, "File path to raw MIP images" );

  boost::filesystem::path pV( "src/algos/algosCellSimilarityMeasure/voronoi" );
  boost::filesystem::path fpV = boost::filesystem::absolute( pV );
  //_voronoiDirectory = std::string( fpV.c_str() );
  _voronoiDirectory = "";
  _profile->addDirectory( _voronoiDirectory, "Voronoi directory" );

  _loadModel = false;
  _profile->addCheckBox( _loadModel, "Load Model Data" );

  _renderCellWalls = false;
  _profile->addCheckBox( _renderCellWalls, "Render Cell Walls (only Model Data)" );

  _cellWallDirectory = "";
  _profile->addFilename( _cellWallDirectory, "Cell Wall File" );

  _modelDataDirectory = "";
  _profile->addFilename( _modelDataDirectory, "Sample Model Data File" );
}

// ---------------------------------------------------------------------

bool SVLRBrowser::canUndo()
{
  return( true );
}

// ---------------------------------------------------------------------

void SVLRBrowser::undo()
{
  SKernel::get()->removeNodeInfo( _vlrInfo.at(0) );
  SKernel::get()->removeNodeInfo( _lineageInfo );
  SKernel::get()->removeNodeInfo( _divisionInfo );
  SKernel::get()->removeNodeInfo( _averagedDivisionInfo );
  SKernel::get()->removeNodeInfo( _divisionAnalysisInfo );
  SKernel::get()->removeNodeInfo( _correlationAnalysisInfo );
  SKernel::get()->removeNodeInfo( _infoBar );
}

// ---------------------------------------------------------------------

std::string SVLRBrowser::getMenuEntry() const
{
  return( std::string( "Algorithms/Data: Biology/Virtual Lateral Root/"
                       + getName() ) );
}

// ---------------------------------------------------------------------

std::string SVLRBrowser::algoName()
{
  return( "VLR Browser" );
}

// ---------------------------------------------------------------------

std::string SVLRBrowser::getName() const
{
  return( algoName() );
}

// ---------------------------------------------------------------------

void SVLRBrowser::execute( void* )
{
  try
  {
    this->disconnectSignals();
    this->print();

    // if model data is loaded then read the layer information from the file
    // without determing any layers or division type information
    if( _loadModel )
    {
      _loadLayerInformation = true;
      _loadDivisionTypeInformation = false;
      _renderCompletePrimordium = false;
      _renderCurvature = false;
      //_cellRadius = 0.05; // model based on bezier patches
      //_cellRadius = 0.17; // model based on single bezier surface
      _cellRadius = 6.; // model based on real data
      //_cellWallWidth = 0.01; // model based on bezier patches
      //_cellWallWidth = 0.02; // model based on single bezier surface
      _cellWallWidth = 1.; // model based on real data
      //_render3DType = 0;
    }

    // pre check if selected data id is smaller than max number of data sets
    if( _sDataId >= _numData )
    {
      std::cout << "Selected data id " << _sDataId
                << " have to be in range [0," << _numData-1
                << "]." << std::endl;
      return;
    }

    /* load lineage data set */

    _lineages.resize( _numData );
    _cellFiles.resize( _numData );
    _divisionScheme.resize( _numData );
    _layerValues.resize( _numData );
    _cellsPerTimeStep.resize( _numData );
    _numTotalTimesteps.resize( _numData );
    _numCellCycles.resize( _numData );
    _centresOfMass.resize( _numData );
    _divisionsPerCellCycle.resize( _numData );
    _rotMatrices.resize( _numData );
    _rotInverseMatrices.resize( _numData );

    _tmin = 1;
    _tmax = 0;
    _ccmax = 0;

    // loop over all data sets
    for( std::size_t d = 0; d < _numData; d++ )
    {
      /* loading lineage data */

      _lineages.at(d) = boost::dynamic_pointer_cast<SAdvancedLineageTreeCollection>(
            SDataManager::get()->getByIndex( _lineageDataId[d] ) );

      if( !_lineages.at(d) )
      {
        std::cout << "Loading data " << d << " failed! Aborting." << std::endl;
        return;
      }

      // determine maximum number of time steps for each dataset
      std::size_t timesteps = 1;
      std::size_t cellCycle = 1;
      for( SAdvancedLineageTreeCollection::const_iterator l = _lineages.at(d)->begin();
           l != _lineages.at(d)->end(); ++l )
      {
        for( auto tree = l->second->begin(); tree != l->second->end(); ++tree )
        {
          if( timesteps <= tree->timeStep )
            timesteps = tree->timeStep;

          if( cellCycle <= tree->getCellCycleId() )
            cellCycle = tree->getCellCycleId();
        }
      }
      _numTotalTimesteps.at(d) = timesteps;
      _numCellCycles.at(d) = cellCycle;

      // store the data set matrices
      SArray rotation = _lineages.at(d)->getRotation();
      osg::Matrix rotMat;
      rotMat.postMult( osg::Matrix::rotate( rotation[2]/180.*M_PI, osg::Vec3(0,0,1) ) );
      rotMat.postMult( osg::Matrix::rotate( rotation[1]/180.*M_PI, osg::Vec3(0,1,0) ) );
      rotMat.postMult( osg::Matrix::rotate( rotation[0]/180.*M_PI, osg::Vec3(1,0,0) ) );
      _rotMatrices.at(d) = rotMat;
      _rotInverseMatrices.at(d) = rotMat.inverse( rotMat );


      // determine max time step range for slider
      if( _lineages.at(d)->getMinimumTimestep() < _tmin )
        _tmin = _lineages.at(d)->getMinimumTimestep();

      if( _numTotalTimesteps.at(d) > (std::size_t)_tmax )
        _tmax = _numTotalTimesteps.at(d);

      if( _numCellCycles.at(d) > (std::size_t)_ccmax )
        _ccmax = _numCellCycles.at(d);
    }

    std::cout << "minT: " << _tmin << " maxT: " << _tmax << " maxCC: " << _ccmax << std::endl;

    /* adjust for tranformations on the dataset */


    if( _lineages.at(_sDataId)->applyRotation() )
    {
      // NOTE that either the dimension OR the center translation should be applied
      // else undefined positions are generated
      SArray center = _lineages.at(_sDataId)->getCenter();
      // The dimension is required to get a correct matching between VLR and the MIP of the raw data
      // NOTE that the dimension translation is only required for the MIP matching and nothing else
      SArray dimension = _lineages.at(_sDataId)->getDimension();
      // Note that we do not divide the z coordinate by 2 since we resampled the data set in this dimension
      // the z coordinate was also flipped, thus no minus sign
      //osg::Matrix rawCenterTranslation =
      //   osg::Matrix::translate( -dimension[0]/2., -dimension[1]/2., dimension[2] );
      // The center translation is applied in order to generate an almost good
      // matching between different data sets such that we get a registration of them
//      osg::Matrix rawCenterTranslation =
//          osg::Matrix::translate( -center[0], -center[1], -center[2] );
//      _vlrRoot = new osg::MatrixTransform( rawCenterTranslation );
      _vlrRoot = new osg::MatrixTransform( osg::Matrix::identity() );

      // then apply the data rotations
      SArray rotation = _lineages.at(_sDataId)->getRotation();
      _vlrRoot->postMult( osg::Matrix::rotate( rotation[2]/180.*M_PI, osg::Vec3(0,0,1) ) );
      _vlrRoot->postMult( osg::Matrix::rotate( rotation[1]/180.*M_PI, osg::Vec3(0,1,0) ) );
      _vlrRoot->postMult( osg::Matrix::rotate( rotation[0]/180.*M_PI, osg::Vec3(1,0,0) ) );
    }
    else
      _vlrRoot = new osg::MatrixTransform( osg::Matrix::identity() );

    _vlrRoot->setCullingActive( false );
    _vlrRoot->setName( "VLR" );

    /* rendering of center */


    // draw center line
    SArray c = _lineages.at(_sDataId)->getCenter();
    osg::Vec3 center( c[0], c[1], c[2] );
    // Note that the center is also transformed according to the
    // rotation matrix stored in the data set, thus the center information
    // is set in the original data set
    center = center * _vlrRoot->getMatrix();
    //SSimilarityMeasureGraphics::render3DCross( center, _vlrRoot );

    // setting of side view normal for division type check between anticlinal
    // or radial
    osg::Vec3 sideViewNormal = osg::Vec3( 0., -1., 0. );
    sideViewNormal = sideViewNormal * _vlrRoot->getInverseMatrix();
    sideViewNormal.normalize();

    _maxLayers = 0;
    // loop over all data sets
    for( std::size_t d = 0; d < _numData; d++ )
    {
      /* determine data type for storing cells per time step */

      _cellsPerTimeStep.at(d) = boost::shared_ptr<cellTimeVector>( new cellTimeVector() );
      _cellsPerTimeStep.at(d)->resize( _numTotalTimesteps.at(d) );
      _centresOfMass.at(d).resize( _numTotalTimesteps.at(d) );
      _divisionsPerCellCycle.at(d).resize( _numCellCycles.at(d) );
      for( auto l = _lineages.at(d)->begin(); l != _lineages.at(d)->end(); ++l )
      {
        SLineageTree *tree = l->second;

//        if( _onlyMasterFiles )
//        {
//          if( _cellFiles.at(d).at( tree->treeId ) != 0 )
//            continue;
//        }

        for( auto nodeIt = tree->begin(); nodeIt != tree->end(); ++nodeIt )
        {
          _cellsPerTimeStep.at(d)->at(nodeIt->timeStep-1).insert( *nodeIt );
          osg::Vec3 pos( nodeIt->getX(), nodeIt->getY(), nodeIt->getZ() );
          //pos = pos * _vlrRoot->getMatrix();
          _centresOfMass.at(d).at(nodeIt->timeStep-1) += pos;

          if( nodeIt->children.size() == 2 )
          {
            // store the division positions for each cell cycle
            std::size_t cellCycle = nodeIt->getCellCycleId();
            _divisionsPerCellCycle.at(d).at(cellCycle).insert( pos );
          }
        }
      }

      for( std::size_t c=0; c < _numTotalTimesteps.at(d); c++ )
      {
        if( _cellsPerTimeStep.at(d)->at(c).size() > 0 )
          _centresOfMass.at(d).at(c) /= _cellsPerTimeStep.at(d)->at(c).size();
      }

      /* load cell file information */

      if( _lineages.at(d)->trackGroupInfoIncluded() )
        _cellFiles.at(d) = _lineages.at(d)->getCellFiles();
      // if there is no cell file information then set the cell
      // coloing based on the cell layers
      else
      {
        std::cerr << "No cell file information is available!" << std::endl;
        _render3DType = 0;
      }

      /* load layer information */

      if( _lineages.at(d)->layerInfoIsIncluded() )
      {
        if( _loadLayerInformation )
        {
          _layerValues.at(d) = _lineages.at(d)->getCellLayers();

          // determine maximum layer value
          for( auto iterT = _layerValues.at(d).begin(); iterT != _layerValues.at(d).end();
               ++iterT )
          {
            for( auto iterC = iterT->begin(); iterC != iterT->end();
                 ++iterC )
            {
              int layer = iterC->second;
              if( layer > _maxLayers )
                _maxLayers = layer;
            }
          }
        }

        std::cout << "Max Layers: " << _maxLayers << std::endl;
      }
      else
        std::cerr << "No layer information is available!" << std::endl;

      /* load division type information */

      if( _lineages.at(d)->divisionInfoIsIncluded() )
      {
        if( _loadDivisionTypeInformation )
          _divisionScheme.at(d) = _lineages.at(d)->getDivisionType();
      }
      else
        std::cerr << "No division information is available!" << std::endl;
    }


    /* initialization of layers */


    // initialize layer colors and cell files colors
    this->setColorMaps();

    if( !_loadLayerInformation )
    {
      // initialize cell layers
      _cellLayers = boost::shared_ptr<SCellLayers>( new SCellLayers( _lineages.at(_sDataId),
                                                                     _changeLayerThreshold,
                                                                     _vlrRoot->getInverseMatrix(),
                                                                     _layerColors,
                                                                     5.0,
                                                                     0.3,
                                                                     _wireframe,
                                                                     _shapeType ) );
    }


    /* initialization of slider */

    {
      // initialize callback for time slider with a waiting thread
      boost::thread t( boost::bind( &SVLRBrowser::bindSlider, this, _sliderLabels.at( _sliderType ) ) );
      t.join();
    }

    /* layer generation w/o loaded layer information */

    // TODO: has to be refactored
    if( _draw3DCells && false )
    {
      // determine layers based on properties set above
      // generate layer information
      if( !_lineages.at(_sDataId)->layerInfoIsIncluded() || !_loadLayerInformation )
      {
        _cellLayers->determineLayers( _vlrRoot, _sliderCallback, _loadLayerInformation, !_renderCurvature );

        // TODO: the exporter should write the layer sequence into the file and
        // not the numbers of my algorithm
        // export layer information into file but only if it is not already included
        if( !_lineages.at(_sDataId)->layerInfoIsIncluded() )
          SLayerInformationIO::writeLayerInformation( _lineages.at(_sDataId)->getFullPath(),
                                                      _cellLayers->getCellLayersPerTimeStep() );
      }
      // else load from file
      else
      {
        if( _lineages.at(_sDataId)->interpolatedData() )
        {
          SLayerInformationIO::readLayerInformation( _lineages.at(_sDataId)->getFullPath(),
                                                     _cellLayers, _lineages.at(_sDataId)->layerColumnIndex() );
        }
        else
        {
          SLayerInformationIO::readLayerInformationPerTimeStep( _lineages.at(_sDataId)->getFullPath(),
                                                                _cellLayers, _lineages.at(_sDataId)->layerColumnIndex() );
        }

        if( !_loadModel )
          _cellLayers->determineLayers( _vlrRoot, _sliderCallback, _loadLayerInformation, !_renderCurvature );
      }

      _cellLayers->printTimeComplexities();
      _cellLayers->transformLayerData();
    }

    // TODO: refactoring
    // process division information: if it does not exist then the layer information is
    // used to generate it
//    for( std::size_t d = 0; d < _numData; d++ )
//    {
//      if( _lineages.at(_sDataId)->layerInfoIsIncluded() )
//        _layerValues.at(d) = _lineages.at(d)->getCellLayers();
//      else
//        std::cerr << "No layer information is available!" << std::endl;

//      if( !_lineages.at(d)->divisionInfoIsIncluded() || !_loadDivisionTypeInformation )
//      {
//        _divisionScheme.at(d).clear();

//        if( !_loadModel )
//          SLayerInformationIO::generateDivisionScheme( _lineages.at(d),
//                                                       _cellLayers,
//                                                       _divisionScheme.at(d),
//                                                       sideViewNormal,
//                                                       _radialAngleThreshold,
//                                                       _vlrRoot->getMatrix() );
//        else
//          SLayerInformationIO::generatePartialDivisionScheme( _lineages.at(d),
//                                                              cellLayers,
//                                                              _divisionScheme.at(d) );

//        if( !_lineages.at(d)->divisionInfoIsIncluded() && !_loadModel )
//        {
//          SLayerInformationIO::writeDivisionScheme( _lineages.at(d)->getFullPath(),
//                                                    _cellLayers->getCellsPerTimeStep(),
//                                                    _divisionScheme.at(d) );
//        }
//      }
//      else
//      {
//        SLayerInformationIO::readDivisionScheme( _lineages.at(d)->getFullPath(),
//                                                 _numTotalTimesteps.at(d),
//                                                 _divisionScheme.at(d) );
//      }
//    }


    // initialize geometry node maps
    _geomNodeMaps = boost::shared_ptr<SGeometryNodeMaps>( new SGeometryNodeMaps() );

    // wait for the rendering to complete everything since we
    // access the media slider which is also done below
    {
      boost::thread t1(  boost::bind( &SVLRBrowser::initializeSlider, this ) );
      t1.join();
    }

    if( _renderDivisionAnalysis )
    {
      boost::thread t( boost::bind( &SVLRBrowser::performDivisionAnalysis, this ) );
      t.join();

      if( _draw3DCells )
      {
        SSimilarityMeasureGraphics::renderCylindricalDivisionAnalysis(
              _lineages.at(_sDataId), _vlrRoot, _sliderCallback,
              _numTotalTimesteps.at( _sDataId ),
              _numCellCycles.at( _sDataId ), _onlyMasterFiles,
              _centresOfMass.at(_sDataId), _sliderType );
      }
    }

    if( _sliderType == 0 && _draw3DCells )
    {
      // wait for the rendering to complete everything since we
      // access the media slider which is also done below
      boost::thread t(  boost::bind( &SVLRBrowser::renderSingleTimeSteps, this ) );
      t.join();
    }
    else if( _sliderType == 1 && _draw3DCells )
    {
      // wait for the rendering to complete everything since we
      // access the media slider which is also done below
      boost::thread t(  boost::bind( &SVLRBrowser::generateDivisionArrows, this ) );
      t.join();
    }

    /* output detailed division information */

    for( std::size_t d = 0; d < _numData; d++ )
      SVLRDataAnalysisIO::printDivisionInformation( _lineages.at(d),
                                                    _numTotalTimesteps.at( d ),
                                                    _centresOfMass.at( d ) );

    // initialize tree visualizer
    SCellTreeVisualizer cTreeVis( _lineages, _numTotalTimesteps,
                                  _divisionScheme, _layerValues,
                                  _cellFiles, _rotMatrices,
                                  _centresOfMass, _layerColors,
                                  _onlyMasterFiles, _renderLineageLineType,
                                  _renderLineageNodeType );

    /* draw 2D lineages with color-coded layer information */

    if( _drawLineages )
    {
      // wait for the lineage rendering
      boost::thread t(  boost::bind( &SCellTreeVisualizer::generateLineageTrees, cTreeVis,
                                     _registerTrees, _divStop,
                                     _lineageLayerGroup,
                                     _lineageInfo, _lineageWindowId ) );
      t.join();
    }


    /* correlation analysis of different cell parameters */

    if( _renderCorrelationAnalysis )
    {
      boost::thread t(  boost::bind( &SVLRBrowser::performCorrelationAnalysis, this ) );
      t.join();
    }

    /* render MIPs and activate volume slicing for validation of cell layers */


    if( _draw3DCells )
    {
      _validation = boost::shared_ptr<SCellLayerValidation>(
            new SCellLayerValidation( _MIPfilename,
                                      _projDir,
                                      _enableSlicing,
                                      _lineages.at(_sDataId),
                                      _sliderCallback,
                                      _vlrRoot ) );
    }

    // generat layer information window in a separate waiting thread
    boost::thread t2( boost::bind( &SVLRBrowser::createLayerWindow, this ) );
    t2.join();


    /* render alpha shape or convex hull of the whole primordium */


    // Note that by default, the shapes are rendered and no output is generated
    if( _renderCompletePrimordium && _draw3DCells )
    {
      // if time slider is chosen then render a surface for each time step
      if( _sliderType == 0 )
      {
        SSimilarityMeasureGraphics::computeAndRenderSurfacesPerTimeStep( true, false, _shapeType,
                                                                     _sliderCallback,
                                                                     _cellsPerTimeStep.at(_sDataId),
                                                                     _vlrRoot, _wireframe, 15. );
      }
      // else render the surface based on the average time step of divisions according
      // to the cell cycles
      else
      {
        std::vector<int> ccVec =
            SSimilarityMeasureUtil::determineAveragedCellCycleTimeSteps( _lineages.at(_sDataId),
                                                                         _numCellCycles.at(_sDataId) );
        SSimilarityMeasureGraphics::computeAndRenderSurfacesPerTimeStep( true, false, _shapeType,
                                                                     _sliderCallback,
                                                                     _cellsPerTimeStep.at(_sDataId),
                                                                     ccVec, _vlrRoot, _wireframe, 15. );
      }
    }

    /* generate voronoi information */


    if( _voronoiDirectory != "" && _draw3DCells )
    {
      // export information for generating voronoi diagram externally
//      std::string path = _voronoiDirectory + "/" + _lineages.at(_sDataId)->getName();
//      SLayerInformationIO::writeCellPositions( path, _cellLayers );

      boost::thread t( boost::bind( &SVLRBrowser::generateVoronoiCells, this ) );
      t.join();
    }


    /* perform cell neighbor length analysis */


    // generate voronoi cells
    if( _voronoiDirectory != "" && _draw3DCells )
    {
      boost::thread t( boost::bind( &SVLRBrowser::lengthAnalysis, this ) );
      t.join();
    }


    /* perform curvature analysis */


    if( _draw3DCells && _renderCurvature )
    {
      boost::thread t( boost::bind( &SVLRBrowser::curvatureAnalysis, this ) );
      t.join();
    }

    /* generate tree division sequence */

    if( _renderTreeDivisionSequence )
    {
      boost::thread t( boost::bind( &SCellTreeVisualizer::generateTreeDivisionSequence, cTreeVis,
                                    _sequenceAxisType, _loadModel,
                                    _divisionInfo, _divisionWindowId ) );
      t.join();
    }

    /* generate tree division sequence */

    if( _loadModel && _modelDataDirectory != "" )
    {
      boost::thread t( boost::bind( &SCellTreeVisualizer::generateAveragedTreeDivisionSequence, cTreeVis,
                                    _modelDataDirectory, _averagedDivisionInfo,
                                    _averagedDivisionWindowId ) );
      t.join();
    }

    this->loadAndGenerateCellShapes();

    /* complete rendering */

    if( _draw3DCells )
    {
      SKernel::get()->render( _vlrRoot, _vlrInfo.at(0), _vlrWindowId.at(0) );

      osg::ref_ptr<osg::Node> colorbar = _lifeDurationColorMap.colorbar( true );
      colorbar->setName( "Color legend" );
      colorbar->getOrCreateStateSet()->setMode( GL_LIGHTING, osg::StateAttribute::OFF);
      SKernel::get()->render( colorbar, _infoBar, _vlrWindowId.at(0) );
      //_vlrRoot->addChild( colorbar );
    }

    SActions::signal_node_t signal = SActions::get( *_vlrRoot ).getNodeSignal( "Create screenshots" );
    _connections.push_back( signal->connect( boost::bind( &SVLRBrowser::generateFrames, this, _vlrRoot, false ) ) );

    //SVLRDataAnalysisIO::appendDivisionProperties( _lineages.at(_sDataId), _rotMatrices.at(_sDataId) );

    //this->performDivisionSequenceAnalysis();
    //this->performDivisionAngleAnalysis();

    //this->exportDivisionOrientations();

    SKernel::get()->setProgress( 1., "Done" );
  }
  CATCH_N_RETHROW( VException );
}

// ---------------------------------------------------------------------

void SVLRBrowser::setColorMaps()
{
  // initialize layer colors
  _layerColors = new osg::Vec4Array;

  int minCellFile = 10;
  int maxCellFile = -10;

  // get max number of cell files for the last time step
  for( std::map<int,int>::const_iterator mapIter =
       _cellFiles.at(_sDataId).begin(); mapIter !=
       _cellFiles.at(_sDataId).end(); mapIter++ )
  {
    if( mapIter->second < minCellFile )
      minCellFile = mapIter->second;

    if( mapIter->second > maxCellFile )
      maxCellFile = mapIter->second;
  }

  std::cout << "Cell files in [" << minCellFile << ", " << maxCellFile << "]." << std::endl;

  if( minCellFile < -4 || maxCellFile > 4 )
    std::cout << "Colors are only supported for cell files in [-4, 4]." << std::endl;

  std::size_t treeCounter = 1;
  for( SAdvancedLineageTreeCollection::const_iterator l = _lineages.at(_sDataId)->begin();
       l != _lineages.at(_sDataId)->end(); ++l )
  {
    SLineageTree *tree = l->second;

    int cellFile = _cellFiles.at(_sDataId).at( tree->treeId );
    bool master = cellFile == 0;

    if( master )
    {
      _treeIdMap[tree->treeId] = treeCounter;
      treeCounter++;
    }
  }

  // at the moment we only allow cell files between -3 and 3
  _cellFileColors[-4] = osg::Vec4( 153./255., 153./255., 153./255., 255./255. );
  _cellFileColors[-3] = osg::Vec4( 228./255.,  26./255.,  28./255., 255./255. );
  _cellFileColors[-2] = osg::Vec4(  55./255., 126./255., 184./255., 255./255. );
  _cellFileColors[-1] = osg::Vec4(  77./255., 175./255.,  74./255., 255./255. );
  _cellFileColors[0] =  osg::Vec4( 152./255.,  78./255., 163./255., 255./255. );
  _cellFileColors[1] =  osg::Vec4( 255./255., 127./255.,   0./255., 255./255. );
  _cellFileColors[2] =  osg::Vec4( 255./255., 255./255.,  51./255., 255./255. );
  _cellFileColors[3] =  osg::Vec4( 166./255.,  86./255.,  40./255., 255./255. );
  _cellFileColors[4] =  osg::Vec4( 247./255., 129./255., 191./255., 255./255. );

  // set colors for lineage trees
  _cellLineageColors[2] = osg::Vec4( 0./255.,   0./255., 255./255., 255./255. );
  _cellLineageColors[1] = osg::Vec4( 0./255.,   255./255., 0./255., 255./255. );
  _cellLineageColors[3] = osg::Vec4( 255./255., 0./255., 0./255., 255./255. );
  _cellLineageColors[4] = osg::Vec4( 255./255., 255./255., 0./255., 255./255. );
  _cellLineageColors[5] = osg::Vec4( 0./255., 255./255., 255./255., 255./255. );
  _cellLineageColors[6] = osg::Vec4( 255./255., 0./255., 255./255., 255./255. );
  _cellLineageColors[7] = osg::Vec4( 255./255., 128./255., 0./255., 255./255. );
  _cellLineageColors[8] = osg::Vec4( 255./255., 0./255., 128./255., 255./255. );
  _cellLineageColors[9] = osg::Vec4( 0./255., 255./255., 0./255., 255./255. );
  _cellLineageColors[10] = osg::Vec4( 0./255.,   0./255., 255./255., 255./255. );
  _cellLineageColors[11] = osg::Vec4( 255./255., 0./255., 0./255., 255./255. );
  _cellLineageColors[12] = osg::Vec4( 255./255., 255./255., 0./255., 255./255. );
  _cellLineageColors[13] = osg::Vec4( 0./255., 255./255., 255./255., 255./255. );
  _cellLineageColors[14] = osg::Vec4( 255./255., 0./255., 255./255., 255./255. );
  _cellLineageColors[15] = osg::Vec4( 255./255., 128./255., 0./255., 255./255. );
  _cellLineageColors[16] = osg::Vec4( 255./255., 0./255., 128./255., 255./255. );

  // colors for layers
  std::size_t numColors = 48;
  double cValues[144] = { 237., 255., 197., 253., 3., 4., 125., 230., 3.,
                          255., 178., 230., 90., 249., 255., 255., 250., 104.,
                          188., 188., 188., 255., 155., 85., 254., 68., 69.,
                          79., 133., 255., 71., 198., 129., 255., 226., 12.,
                          193., 67., 255., 180., 124., 221., 255., 112., 194.,
                          166., 206., 227., 31., 120., 180., 178., 223., 138.,
                          51., 160., 44., 251., 154., 153.,
                          227., 26., 28., 253., 191., 111., 255., 127., 0.,
                          202., 178., 214., 106., 61., 154., 255., 255., 153.,
                          177., 89., 40., 141., 211., 199., 255., 255., 179.,
                          190., 186., 218., 251., 128., 114., 128., 177., 211.,
                          253., 180., 98., 179., 222., 105., 252., 205., 229.,
                          217., 217., 217., 188., 128., 189., 204., 235., 197.,
                          255., 237., 111., 228., 26., 28., 55., 126., 184.,
                          77., 175., 74., 152., 78., 163., 255., 127., 0.,
                          255., 255., 51., 166., 86., 40., 247., 129., 191.,
                          153., 153., 153. };

  std::size_t i = 0;
  for( std::size_t c=0; c < numColors; c++, i+=3 )
    _layerColors->push_back( osg::Vec4( cValues[i]/255., cValues[i+1]/255., cValues[i+2]/255., 255./255. ) );

  // initialize area color map
  _areaColorMap.setBlueWhiteRed();
  _areaColorMap.scaleColors( _lineages.at(_sDataId)->getMinMaxArea().first,
                             _lineages.at(_sDataId)->getMinMaxArea().second );
  _areaColorMap.setNumberColorbarTicks( 8 );

  std::size_t min, max;
  std::map<std::size_t, std::size_t> lifeDuration;
  SSimilarityMeasureUtil::determineCellLifeDuration( _lineages.at(_sDataId),
                                                     _centresOfMass.at(_sDataId),
                                                     _cellFiles.at(_sDataId),
                                                     _rotMatrices.at(_sDataId),
                                                     lifeDuration,
                                                     min, max );
  _lifeDurationColorMap.setCoolToWarm();
  _lifeDurationColorMap.scaleColors( min, max );
  _lifeDurationColorMap.setNumberColorbarTicks( 2 );
}

// ---------------------------------------------------------------------

void SVLRBrowser::loadAndGenerateCellShapes()
{
  osg::ref_ptr<osg::Group> shapes = new osg::Group;
  std::vector< std::vector<osg::Vec3> > vertices;
  std::string fileName = "/home/necrolyte/Dropbox/ShapeVerticesLOD_20160427_T010.txt";
  SVLRDataAnalysisIO::readCellShapeVertices( fileName, vertices );

  SInteractiveDivisionArrows *is = new SInteractiveDivisionArrows();
  osg::ref_ptr<osg::Group> volumeGroup = new osg::Group;
  is->setName( "Alpha Cell Shapes" );
  is->addUpdateCallback( _sliderCallback );
  is->addChild( volumeGroup );
  unsigned int numCells = vertices.size();
  std::cout << "numCells: " << numCells << std::endl;
  unsigned int cellCounter = 0;
  for( auto iter1 = vertices.begin(); iter1<vertices.end(); iter1++ )
  {
    osg::Vec4 color = (*_layerColors)[cellCounter%_layerColors->size()];

    osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array;
    for( auto iter2 = iter1->begin(); iter2<iter1->end(); iter2++ )
      points->push_back( osg::Vec3( iter2->x(), iter2->y(), iter2->z() ) );

    boost::shared_ptr<SConvexHull> tri = boost::shared_ptr<SConvexHull>(
          new SConvexHull( points, 1, color, true, false ) );

    shapes->addChild( tri->getTriangles() );
    cellCounter++;
  }

  is->addDivisionArrows( 1, shapes );

  _vlrRoot->addChild( is );
}

// ---------------------------------------------------------------------

void SVLRBrowser::performCorrelationAnalysis()
{
  osg::ref_ptr<osg::Group> corGroup = new osg::Group;
  corGroup->setName( "ScatterplotAnalysis" );
  // correlation analysis of division data for each division type
  for( std::size_t divType = 0; divType < 3; divType++ )
  {
    std::size_t dim = 12;
    SCorrelationAnalysis corAnalysis( _lineages, dim );
    corAnalysis.generateDivisionData( divType, _tmax );
    osg::ref_ptr<osg::MatrixTransform> mat = new osg::MatrixTransform;
    mat->setMatrix( osg::Matrix::translate( osg::Vec3( (double)divType * (dim+8), 0., 0. ) ) );
    corAnalysis.renderScatterplotMatrix( mat );
    corGroup->addChild( mat );
  }

  // general correlation analysis for each nuclei
//  std::size_t dim = 8;
//  SCorrelationAnalysis corAnalysis( _lineages, dim );
//  corAnalysis.generateGeneralData( _tmax );
//  osg::ref_ptr<osg::MatrixTransform> mat = new osg::MatrixTransform;
//  mat->setMatrix( osg::Matrix::identity() );
//  corAnalysis.renderScatterplotMatrix( mat );
//  corGroup->addChild( mat );

//  std::string name = "/tmp/" + corGroup->getName() + ".svg";
//  SSVGRenderer svg( true );
//  svg.startImage( name );
//  svg.renderOSGNode( corGroup, 1., 1., 2. );

  SKernel::get()->render( corGroup, _correlationAnalysisInfo, _correlationAnalysisWindowId );
}

// ---------------------------------------------------------------------

void SVLRBrowser::performDivisionAnalysis()
{
  _infoBar = SNodeInfo( "Color legend", this, -1 );
  _divGroup = new osg::Group;
  _divGroup->getOrCreateStateSet()->setMode( GL_LIGHTING, osg::StateAttribute::OFF);

  double xOff[2];
  xOff[0] = 770.;
  xOff[1] = 530.;

  double yOffset = 0.;
  double tileSize;

  std::vector<SDivisionAnalysis> divAnalysisVec;
  std::pair<double, double> globalMinMax;
  globalMinMax.first = 1000.;
  globalMinMax.second = -1000.;
  std::vector< std::vector<osg::ref_ptr<osg::MatrixTransform> > > contourMat;
  for( std::size_t d = 0; d < _numData; d++ )
  {
    std::vector<osg::ref_ptr<osg::MatrixTransform> > contour;
    for( std::size_t v=0; v < 3; v++ )
    {
      contour.push_back( new osg::MatrixTransform );
      contour.back()->setMatrix( _rotMatrices.at(d) );
    }

    SDivisionAnalysis divAnalysis( _lineages.at(d),
                                   _numTotalTimesteps.at(d),
                                   _numCellCycles.at(d),
                                   _regSteps,
                                   _onlyMasterFiles,
                                   _sliderCallback,
                                   contour, _sliderType,
                                   _angleTypeDivAnalysis );

    contourMat.push_back( contour );

    std::pair<double, double> minMax = divAnalysis.getMinMaxValues();

    if( minMax.first < globalMinMax.first )
      globalMinMax.first = minMax.first;

    if( minMax.second >= globalMinMax.second )
      globalMinMax.second = minMax.second;

    divAnalysisVec.push_back( divAnalysis );
  }

  std::cout << "min: " << globalMinMax.first << " max: " << globalMinMax.second << std::endl;

  SColorMapTF cmtf( SColorMapTF::LINEAR );
  cmtf.setNumberOfTicks( 2 );
  std::map<double,SColor> &colorList = cmtf.getColorList();
  colorList.clear();
  colorList[0./7.] = osg::Vec4( 255./255., 255./255., 204./255., 1. );
  colorList[1./7.] = osg::Vec4( 255./255., 237./255., 160./255., 1. );
  colorList[2./7.] = osg::Vec4( 254./255., 217./255., 118./255., 1. );
  colorList[3./7.] = osg::Vec4( 254./255., 178./255., 76./255., 1. );
  colorList[4./7.] = osg::Vec4( 253./255., 141./255., 60./255., 1. );
  colorList[5./7.] = osg::Vec4( 252./255., 78./255., 42./255., 1. );
  colorList[6./7.] = osg::Vec4( 227./255., 26./255., 28./255., 1. );
  colorList[7./7.] = osg::Vec4( 177./255., 0./255., 38./255., 1. );
  //cm.setQualitative9();
  cmtf.setRange( globalMinMax.first, globalMinMax.second );

  osg::Node *colorbar = cmtf.colorbar( true );
  colorbar->setName( "Color legend" );
  colorbar->getOrCreateStateSet()->setMode( GL_LIGHTING, osg::StateAttribute::OFF);

  // scale colors of colormap
  SColorMap colorMap;
  colorMap.getColorList() = colorList;
  colorMap.scaleColors( globalMinMax.first, globalMinMax.second );

  for( std::size_t d = 0; d < _numData; d++ )
  {
    double xOffset = 0.;

    osg::ref_ptr<osg::MatrixTransform> dataGroup = new osg::MatrixTransform;
    dataGroup->postMult( osg::Matrix::translate( osg::Vec3( 0., -yOffset, 0. ) ) );

    // print data name
    osg::ref_ptr<osg::Geode> labels = new osg::Geode;
    osg::ref_ptr<osgText::Text> label = SLabelGenerator::getLabel(
          _lineages.at(d)->getName(), osg::Vec3( -380., 100., 0. ),
          osg::Vec4(0.2,0.2,0.2,1.), 20., osgText::Text::LEFT_TOP);

    label->setCharacterSizeMode( osgText::Text::OBJECT_COORDS );
    label->setAutoRotateToScreen( false );

    labels->addDrawable( label );
    dataGroup->addChild( labels );

    std::vector< std::set<Angles, setComp> > vec = divAnalysisVec.at(d).getDivisionProperties();
    std::vector< std::set<Angles, setComp> > regVec = divAnalysisVec.at(d).getRegisteredDivisionProperties();
    std::vector< std::vector< std::pair<double,double> > > bound = divAnalysisVec.at(d).getMinMaxBoundary();
    std::vector< osg::ref_ptr<osg::MatrixTransform> > gridGroup;
    gridGroup.resize(3);
    std::vector< osg::ref_ptr<osg::MatrixTransform> > gridAndContourGroup;
    gridAndContourGroup.resize(3);
    for( std::size_t v=0; v < 3; v++ )
    {
      gridGroup.at(v) = new osg::MatrixTransform;
      gridAndContourGroup.at(v) = new osg::MatrixTransform;
      SScalarGrid sg( _gridResolution, _numTotalTimesteps.at(d),
                      _numCellCycles.at(d), _regSteps, _sliderCallback,
                      bound.at(v), colorMap, v, _angleTypeDivAnalysis, _sliderType );

      // apply the same tile size to each viewing type
      if( v == 0 && d == 0 )
      {
        sg.createGrid();
        tileSize = sg.getTileSize();
      }
      else
      {
        sg.setTileSize( tileSize );
        sg.createGrid();
      }

      if( _sliderType == 2 )
        sg.setTileValues( regVec, gridGroup.at(v) );
      else
        sg.setTileValues( vec, gridGroup.at(v) );

      osg::Vec3 s( bound.at(v).at(0).first, bound.at(v).at(1).first, bound.at(v).at(2).first );
      osg::Vec3 e( bound.at(v).at(0).second, bound.at(v).at(1).second, bound.at(v).at(2).second );
      osg::Vec3 c = s + (e-s)/2.;

//      std::cout << "c: " << c << std::endl;

      gridAndContourGroup.at(v)->addChild( gridGroup.at(v) );
      gridAndContourGroup.at(v)->addChild( contourMat.at(d).at(v) );

      if( v == 0 )
      {
        // translate to origin
        gridAndContourGroup.at(v)->postMult( osg::Matrix::translate( osg::Vec3( -c[0], -c[1], 0. ) ) );
        gridAndContourGroup.at(v)->setName( "Front" );
        xOffset += xOff[0];//(e[0]-s[0]) + 25.;
      }
      else if( v == 1 )
      {
        // translate to origin
        contourMat.at(d).at(v)->postMult( osg::Matrix::rotate( -90./180.*M_PI, osg::Vec3(1,0,0) ) );
        gridAndContourGroup.at(v)->postMult( osg::Matrix::translate( osg::Vec3( -c[0], -c[2], 0. ) ) );
        gridAndContourGroup.at(v)->setName( "Side" );
        // translate by xOffset
        gridAndContourGroup.at(v)->postMult( osg::Matrix::translate( osg::Vec3( xOffset, 0., 0. ) ) );
        xOffset += xOff[1];//3.*(e[0]-s[0])/4. + 25.;
      }
      else
      {
        // translate to origin
        contourMat.at(d).at(v)->postMult( osg::Matrix::rotate( -90./180.*M_PI, osg::Vec3(0,0,1) ) );
        contourMat.at(d).at(v)->postMult( osg::Matrix::rotate( -90./180.*M_PI, osg::Vec3(1,0,0) ) );
        gridAndContourGroup.at(v)->postMult( osg::Matrix::translate( osg::Vec3( -c[1], -c[2], 0. ) ) );
        gridAndContourGroup.at(v)->setName( "Radial" );
        // translate by xOffset
        gridAndContourGroup.at(v)->postMult( osg::Matrix::translate( osg::Vec3( xOffset, 0., 0. ) ) );
        yOffset += (e[2]-s[2]) + 100.;
      }

      dataGroup->addChild( gridAndContourGroup.at(v) );
    }

    _divGroup->addChild( dataGroup );
  }

  _divGroup->setName( "DivisionAnalysis" );

  SKernel::get()->render( colorbar, _infoBar, _divisionAnalysisWindowId );

  SActions::signal_node_t signal1 = SActions::get( *_divGroup ).getNodeSignal( "Create screenshots" );
  _connections.push_back( signal1->connect( boost::bind( &SVLRBrowser::generateFrames, this, _divGroup, true ) ) );

  SKernel::get()->render( _divGroup, _divisionAnalysisInfo, _divisionAnalysisWindowId );
}

// ---------------------------------------------------------------------

void SVLRBrowser::performDivisionAngleAnalysis()
{
  bool med = false;
  bool chooseTheta = false;

  std::string fileName = "/home/necrolyte/Uni/LateralRootGrowth/211215/model/";
  std::vector< std::map<std::size_t, FileAngles> > realAngles;
  for( std::size_t d = 0; d < _numData; d++ )
  {
    std::string dataName = fileName + "withOrigPositions/DivisionAnalysis_";
    dataName += _lineages.at(d)->getName();
    dataName += ".csv";
    std::map<std::size_t, FileAngles> dataAngles;
    SVLRDataAnalysisIO::readAngleData( dataName, dataAngles );
    realAngles.push_back( dataAngles );
  }

  // combine the angle results of the real data
  std::vector< std::map< int, std::vector<double> > > rAngles;
  rAngles.resize( 300 );
  std::vector< std::map<int, double> > rAveragedValues;
  rAveragedValues.resize( 300 );

  for( auto iter = realAngles.begin(); iter != realAngles.end(); ++iter )
  {
    for( auto divIter = iter->begin(); divIter != iter->end(); ++divIter )
    {
      std::size_t numCells = divIter->first;
      int cellFile = divIter->second.cellFile;
      std::size_t divType = divIter->second.divType;

      if( cellFile < -2 || cellFile > 2 )
        continue;

      // only consider pericinal and radial divisions
      if( divType == 0 )
        continue;

      double theta = divIter->second.theta;
      double phi = divIter->second.phi;

      auto mapIter = rAngles.at( numCells ).find( cellFile );

      if( mapIter != rAngles.at( numCells ).end() )
      {
        if( chooseTheta )
          mapIter->second.push_back( theta );
        else
          mapIter->second.push_back( phi );
      }
      else
      {
        std::vector<double> angValues;
        if( chooseTheta )
          angValues.push_back( theta );
        else
          angValues.push_back( phi );

        rAngles.at( numCells ).insert( std::make_pair( cellFile, angValues ) );
      }
    }
  }

  // afterwards get the median of angle values for each cell file
  // and draw point
  // for each number of cells
  std::size_t vecI = 0;
  for( auto iter = rAngles.begin(); iter != rAngles.end(); ++iter, vecI++ )
  {
    for( auto fileIter = iter->begin(); fileIter != iter->end(); ++fileIter )
    {
      int cellFile = fileIter->first;
      std::vector<double> values = fileIter->second;
      std::sort( values.begin(), values.end() );

      double res;

      // median
      if( med )
      {
        if( values.size()%2 == 0 )
        {
          std::size_t i1 = (double)(values.size())/2. - 1;
          std::size_t i2 = (double)(values.size())/2.;
          res = 0.5 * ( values.at(i1) + values.at(i2) );
        }
        else
        {
          std::size_t i = (double)(values.size()+1)/2. - 1;
          res = values.at(i);
        }
      }
      // mean
      else
      {
        double sum = 0.;
        for( std::size_t v=0; v < values.size(); v++ )
          sum += values.at(v);

        res = sum/(double)values.size();
      }

      rAveragedValues.at( vecI ).insert( std::make_pair( cellFile, res ) );
    }
  }

  std::string exportName = "/tmp/RealDataAngleAnalysis";
  if( med )
    exportName += "Median";
  else
    exportName += "Mean";

  if( chooseTheta )
    exportName += "Theta";
  else
    exportName += "Phi";

  exportName += ".csv";
  SVLRDataAnalysisIO::exportAngleAnalysis( exportName, rAveragedValues );

  std::string modelName = fileName + "divisionAnglesRadialModelBesson-Dumais.csv";
  // first read angle data of model runs
  std::vector< std::map<std::size_t, FileAngles> > modelAngles;
  SVLRDataAnalysisIO::readModelAngleData( modelName, modelAngles );

  std::vector< std::map<int, double> > averagedValues;
  averagedValues.resize( 80 );

  // average the results of the model runs
  std::vector< std::map< int, std::vector<double> > > angles;
  angles.resize( 80 );

  // for each run
  for( auto iter = modelAngles.begin(); iter != modelAngles.end(); ++iter )
  {
    // for each division
    for( auto divIter = iter->begin(); divIter != iter->end(); ++divIter )
    {
      std::size_t numCells = divIter->first;
      int cellFile = divIter->second.cellFile;
      double theta = divIter->second.theta;

      auto mapIter = angles.at( numCells ).find( cellFile );

      if( mapIter != angles.at( numCells ).end() )
        mapIter->second.push_back( theta );
      else
      {
        std::vector<double> angValues;
        angValues.push_back( theta );
        angles.at( numCells ).insert( std::make_pair( cellFile, angValues ) );
      }
    }
  }

  // afterwards get the median of angle values for each cell file
  // and draw point
  // for each number of cells
  vecI = 0;
  for( auto iter = angles.begin(); iter != angles.end(); ++iter, vecI++ )
  {
    for( auto fileIter = iter->begin(); fileIter != iter->end(); ++fileIter )
    {
      int cellFile = fileIter->first;
      std::vector<double> values = fileIter->second;
      std::sort( values.begin(), values.end() );

      double res;

      // median
      if( med )
      {
        if( values.size()%2 == 0 )
        {
          std::size_t i1 = (double)(values.size())/2. - 1;
          std::size_t i2 = (double)(values.size())/2.;
          res = 0.5 * ( values.at(i1) + values.at(i2) );
        }
        else
        {
          std::size_t i = (double)(values.size()+1)/2. - 1;
          res = values.at(i);
        }
      }
      // mean
      else
      {
        double sum = 0.;
        for( std::size_t v=0; v < values.size(); v++ )
          sum += values.at(v);

        res = sum/(double)values.size();
      }

      averagedValues.at( vecI ).insert( std::make_pair( cellFile, res ) );
    }
  }

  exportName = "/tmp/ModelAngleAnalysis";
  if( med )
    exportName += "Median";
  else
    exportName += "Mean";

  exportName += ".csv";
  SVLRDataAnalysisIO::exportAngleAnalysis( exportName, averagedValues );
}

// ---------------------------------------------------------------------

void SVLRBrowser::performDivisionSequenceAnalysis()
{
  bool divTypePCA = false;
  std::map<std::string, std::vector<std::size_t> > divisionSequences;
  std::size_t numRadialDiv = 0;

  for( std::size_t d = 0; d < _numData; d++ )
  {
    std::map< std::pair<std::size_t, std::size_t>, std::size_t > divTypes;

//    if( divTypePCA )
//      SVLRDataAnalysisIO::readAngleData( _lineages.at(d), divTypes );

    std::size_t firstTime = 350;
    std::size_t firstId = 1;
    std::size_t firstFile = 0;
    // determine id and time step of first radial division
    for( SAdvancedLineageTreeCollection::const_iterator l = _lineages.at(d)->begin();
         l != _lineages.at(d)->end(); ++l )
    {
      for( SLineageTree::const_iterator iter = l->second->begin();
           iter != l->second->end(); ++iter )
      {
        if( iter->children.size() == 2 )
        {
          std::size_t id = iter->cellId;
          std::size_t t = iter->timeStep;
          bool found;
          std::size_t divType;
          if( divTypePCA )
          {
            auto divIter = divTypes.find( std::make_pair( id, t ) );
            found = divIter != divTypes.end();
            if( found )
              divType = divIter->second;
          }
          else
          {
            auto divIter = _divisionScheme.at(d).at(t-1).find( id );
            found = divIter != _divisionScheme.at(d).at(t-1).end();
            if( found )
              divType = divIter->second;
          }

          if( found )
          {
            if( divType == 2 && t < firstTime )
            {
              firstTime = t;
              firstId = id;
              firstFile = _cellFiles.at(d).at( iter->treeId );
            }
          }
        }
      }
    }

    for( SAdvancedLineageTreeCollection::const_iterator l = _lineages.at(d)->begin();
         l != _lineages.at(d)->end(); ++l )
    {
      for( SLineageTree::const_iterator iter = l->second->begin();
           iter != l->second->end(); ++iter )
      {
        if( iter->children.size() == 2 )
        {
          std::string seq = "";
          std::size_t id = iter->cellId;
          std::size_t t = iter->timeStep;

          bool found;
          std::size_t divType;
          if( divTypePCA )
          {
            auto divIter = divTypes.find( std::make_pair( id, t ) );
            found = divIter != divTypes.end();
            if( found )
              divType = divIter->second;
          }
          else
          {
            auto divIter = _divisionScheme.at(d).at(t-1).find( id );
            found = divIter != _divisionScheme.at(d).at(t-1).end();
            if( found )
              divType = divIter->second;
          }

          if( found )
          {
            // analyze division sequence before each radial division
            if( divType == 2 )
            {
              numRadialDiv++;
              auto iter2 = iter->parent;
              while( iter2->parent )
              {
                iter2 = iter2->parent;
                if( iter2->children.size() == 2 )
                {
                  std::size_t idp = iter2->cellId;
                  std::size_t tp = iter2->timeStep;

                  bool foundI;
                  std::size_t divTypep;
                  if( divTypePCA )
                  {
                    auto divIter = divTypes.find( std::make_pair( idp, tp ) );
                    foundI = divIter != divTypes.end();
                    if( foundI )
                      divTypep = divIter->second;
                  }
                  else
                  {
                    auto divIter = _divisionScheme.at(d).at(tp-1).find( idp );
                    foundI = divIter != _divisionScheme.at(d).at(tp-1).end();
                    if( foundI )
                      divTypep = divIter->second;
                  }

                  if( foundI )
                  {
                    switch( divTypep )
                    {
                    case 0: seq.insert(0, "A"); break;
                    case 1: seq.insert(0, "P"); break;
                    case 2: seq.insert(0, "R"); break;
                    }
                  }
                }
              }
            }
          }

          if( seq != "" )
          {
            int cF = _cellFiles.at(d).at( iter->treeId );
            auto seqIter = divisionSequences.find( seq );
            if( seqIter != divisionSequences.end() )
            {
              seqIter->second.at( cF + 2 )++;
            }
            else
            {
              std::vector<std::size_t> vecSeq;
              vecSeq.resize( 5, 0 );
              vecSeq.at( cF + 2 )++;
              divisionSequences.insert( std::make_pair( seq, vecSeq ) );
            }

            // check if this radial division is the first one
            if( id == firstId && firstTime == t )
              std::cout << "first radial sequence: " << seq
                        << " in cell file " << firstFile
                        << " for dataset " << _lineages.at(d)->getName() << std::endl;
          }
        }
      }
    }
  }

  std::string fileName = "/tmp/DivisionSequence.csv";
  // open file at the end of the stream and allow writing
  std::ofstream out( fileName.c_str(), std::ofstream::out );

  if( !out.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to write into data " + fileName );
    return;
  }

  out << "CellFile -2 -1 0 1 2\n";

  for( auto iter = divisionSequences.begin(); iter != divisionSequences.end(); ++iter )
  {
    out << iter->first;
    for( std::size_t f = 0; f < 5; f++ )
      out << " " << iter->second.at(f);

    out << "\n";
  }

  std::cout << "Total radial divisions: " << numRadialDiv << std::endl;

  out.close();
}

// ---------------------------------------------------------------------

void SVLRBrowser::generateDivisionArrows()
{
  SSimilarityMeasureGraphics::renderDivisions( _lineages.at(_sDataId),
                                               _vlrRoot, _sliderCallback,
                                               _layerColors, _layerValues.at(_sDataId),
                                               _renderTracks, _onlyMasterFiles,
                                               _divisionColorType, _renderArrows );
}

// ---------------------------------------------------------------------

void SVLRBrowser::curvatureAnalysis()
{
  std::vector<std::pair<double,double> > minMaxCurvatures;
  minMaxCurvatures.resize( 2, std::make_pair( boost::numeric::bounds<double>::highest(),
                                              boost::numeric::bounds<double>::lowest() ) );

  curvatureType type;

  if( _renderCurvature && _curvatureType == 0 )
    type = MEANCURVATURE;
  else if( _renderCurvature && _curvatureType == 1 )
    type = GAUSSIANCURVATURE;
  else
    type = NOCURVATURE;

  boost::shared_ptr<SCurvatureAnalysis> curvatureAnalysis =
      boost::shared_ptr<SCurvatureAnalysis>( new SCurvatureAnalysis( _sliderCallback,
                                                                     _cellsPerTimeStep.at( _sDataId ),
                                                                     type, minMaxCurvatures,
                                                                     _shapeType,
                                                                     _vlrRoot ) );

  for( std::size_t c = 0; c < 2; c++ )
  {
    std::cout << "MinMax: " << minMaxCurvatures.at(c).first << " "
              << minMaxCurvatures.at(c).second << std::endl;
  }

  SColorMap cm;
  SColorMap::ColorList &cl = cm.getColorList();
  cl.clear();
  cl[ 0.0/4.0 ] = SColor(0,0,1);
  cl[ 1.0/4.0 ] = SColor(0,1,1);
  cl[ 2.0/4.0 ] = SColor(0,1,0);
  cl[ 3.0/4.0 ] = SColor(1,1,0);
  cl[ 4.0/4.0 ] = SColor(1,0,0);
  cm.scaleColors( minMaxCurvatures.at(_curvatureType).first, minMaxCurvatures.at(_curvatureType).second );
  cm.setNumberColorbarTicks( 2 );
  osg::Node *colorbar = cm.colorbar( false );

  colorbar->setName( "Color legend" );
  colorbar->getOrCreateStateSet()->setMode( GL_LIGHTING, osg::StateAttribute::OFF);
  _colorbarInfo = SNodeInfo( "Color legend", this, -1 );

  SKernel::get()->render( colorbar, _colorbarInfo, _vlrWindowId.at(0) );
}

// ---------------------------------------------------------------------

void SVLRBrowser::lengthAnalysis()
{
  boost::shared_ptr<SCellAnalysis> analysis =
      boost::shared_ptr<SCellAnalysis>( new SCellAnalysis( _voronoi, _lineages.at(_sDataId),
                                                           _sliderCallback, _vlrRoot, true,
                                                           _sliderType == 1 ) );

  bool exportAnalysis = false;

  if( exportAnalysis )
  {
    std::string path = "/home/necrolyte/Uni/LateralRootGrowth/Status_200214/GrowthAnalysisInterpolated"
        + _lineages.at(_sDataId)->getName() + ".csv";

    //analysis->writeAnalysisResults( path );
    analysis->writeSingleDisplacements( path );
  }

//  std::string path = "/home/necrolyte/scifer/src/algos/algosCellSimilarityMeasure/DescriptorInformation/LengthParameters/lengthParameters"
//      + _lineages.at(_sDataId)->getName();
//  path += ".txt";

//  analysis->writeLengthParameters( path );
}

// ---------------------------------------------------------------------

void SVLRBrowser::generateFrames( osg::ref_ptr<osg::Group> renderGroup,
                                  const bool exportAsSVG )
{
  // NOTE that contiuous rendering must be activated for the window in which
  // the osg node is assigned such that the callback of the screenshot method
  // is called

  SScreenshot::ScreenshotSettings sss;
  sss.doNotUseCurrentView = false;
  sss.resolution = 2000;
  sss.enableLighting = false;

  std::size_t maxTimeStep;
  if( _sliderType == 0 )
    maxTimeStep = _tmax;
  else if( _sliderType == 1 )
    maxTimeStep = _ccmax;
  else
    maxTimeStep = _regSteps;

  std::string digitsString;
  std::string dataName = _lineages.at(_sDataId)->getName();
  boost::filesystem::path dir( "/tmp/" + dataName );
  boost::filesystem::create_directory( dir );

  for( std::size_t t = 0; t < maxTimeStep; t++ )
  {
    SKernel::get()->setSliderPosition( _sliderLabels.at( _sliderType ), t+1 );
    SKernel::get()->setProgress( (double)t/(double)(maxTimeStep - 1),
                                 "Creating screenshots" );

    if( t < 10 )
      digitsString = "00";
    else if( t < 100 )
      digitsString = "0";
    else //if( t < 1000 )
      digitsString = "";

    std::string fileName = renderGroup->getName() + "_" + digitsString + boost::lexical_cast<std::string>(t);
    SScreenshot::takeScreenshot( renderGroup, "/tmp/" + dataName, fileName, "png", sss );

    if( exportAsSVG )
    {
      int i = SKernel::get()->getSliderPosition();
      std::string name = "/tmp/" + dataName + "/" + renderGroup->getName() + "_" + std::to_string( i ) + ".svg";
      SSVGRenderer svg( true );
      svg.startImage( name );
      svg.renderOSGNode( renderGroup, 1., 1., 2. );
    }
  }

  SKernel::get()->setProgress( 1, "Done" );
}

// ---------------------------------------------------------------------

void SVLRBrowser::generateVoronoiCells()
{
  std::string path = _voronoiDirectory + "/" + _lineages.at(_sDataId)->getName();

  _voronoi = boost::shared_ptr<SVoronoiDiagram>( new SVoronoiDiagram( path, _lineages.at(_sDataId),
                                                                      _sliderCallback,
                                                                      _cellsPerTimeStep.at(_sDataId) ) );

  _vlrRoot->addChild( _voronoi->getVoronoi() );
}

// ---------------------------------------------------------------------

void SVLRBrowser::initializeSlider()
{
  if( _sliderType == 0 )
  {
    SKernel::get()->resetSlider( _sliderLabels.at( _sliderType ), _tmin, _tmax, _tfoc );
    // set the initial position of the slider to 1
    SKernel::get()->setSliderPosition( _sliderLabels.at( _sliderType ), 1 );
  }
  else if( _sliderType == 1 )
  {
    SKernel::get()->resetSlider( _sliderLabels.at( _sliderType ), _tmin, _ccmax, _tfoc );
    // set the initial position of the slider to 1
    SKernel::get()->setSliderPosition( _sliderLabels.at( _sliderType ), 1 );
  }
  else if( _sliderType == 2 )
  {
    SKernel::get()->resetSlider( _sliderLabels.at( _sliderType ), _tmin, _regSteps, _tfoc );
    SKernel::get()->setSliderPosition( _sliderLabels.at( _sliderType ), 1 );
  }
}

// ---------------------------------------------------------------------

void SVLRBrowser::apply()
{
  if( !_enableSlicing )
    std::cout << "Slicing is not updated since it is not enabled in the dialog." << std::endl;

  if( _MIPfilename == "" )
    std::cout << "MIP is not updated since it is not enabled in the dialog." << std::endl;

  // if both are disabled then do nothing at all
  if( !_enableSlicing && _MIPfilename == "" )
    return;

  // get index of new projection
  int index;

  if( _layerWindow )
  {
    QMetaObject::invokeMethod( _layerWindow,
                               "getProjection",
                               Q_ARG( int&, index ) );
  }

  _validation->update( index );

  SKernel::get()->requestRedraw();
}

// ---------------------------------------------------------------------

void SVLRBrowser::createLayerWindow()
{
  // at last enable layer window
  if( !_signalsConnected )
  {
    _layerWindowCreator->init();
    while( !_layerWindowCreator->widget() );

    // get pointer to descriptor widget for accessing several methods
    _layerWindow = dynamic_cast<SCellLayerWindow*>( _layerWindowCreator->widget() );

    // handling of signals with descriptor widget
    _connections.push_back( _layerWindowCreator->onSaveLayersClicked(
                              boost::bind( &SVLRBrowser::saveLayers, this ) ) );

    _connections.push_back( _layerWindowCreator->onCreateScreenshotsClicked(
                              boost::bind( &SVLRBrowser::generateFrames, this, _vlrRoot, false ) ) );

    _connections.push_back( _layerWindowCreator->onApplyClicked(
                              boost::bind( &SVLRBrowser::apply, this ) ) );

    _signalsConnected = true;
  }

  if( _layerWindow )
  {
    QMetaObject::invokeMethod( _layerWindow,
                               "setProjections",
                               Q_ARG( const std::list<std::string>, _projections ),
                               Q_ARG( const int, _projDir ) );

    QMetaObject::invokeMethod( _layerWindow, "show" );
  }
}

// ---------------------------------------------------------------------

void SVLRBrowser::setCellFeatures( const int cellId,
                                   const int timeStep,
                                   const int treeId,
                                   const int layer,
                                   const cellHistory cellHist,
                                   const double xValue,
                                   const double yValue,
                                   const double zValue,
                                   const int cellFile,
                                   const std::size_t cellDivisionType )
{
  if( _layerWindow )
  {
    if( cellHist.size() != 0 )
    {
      QMetaObject::invokeMethod( _layerWindow,
                                 "setCellFeatures",
                                 Q_ARG( const int, cellId ),
                                 Q_ARG( const int, timeStep ),
                                 Q_ARG( const int, treeId ),
                                 Q_ARG( const int, layer ),
                                 Q_ARG( const cellHistory, cellHist ),
                                 Q_ARG( const double, xValue ),
                                 Q_ARG( const double, yValue ),
                                 Q_ARG( const double, zValue ),
                                 Q_ARG( const int, cellFile ),
                                 Q_ARG( const std::size_t, cellDivisionType ) );
    }
    else
    {
      QMetaObject::invokeMethod( _layerWindow,
                                 "setCellFeatures",
                                 Q_ARG( const int, cellId ),
                                 Q_ARG( const int, timeStep ),
                                 Q_ARG( const int, treeId ),
                                 Q_ARG( const int, layer ),
                                 Q_ARG( const double, xValue ),
                                 Q_ARG( const double, yValue ),
                                 Q_ARG( const double, zValue ),
                                 Q_ARG( const int, cellFile ),
                                 Q_ARG( const std::size_t, cellDivisionType ) );
    }
  }
}

// ---------------------------------------------------------------------

void SVLRBrowser::saveLayers()
{
  // TODO
  if( _draw3DCells && false )
    SLayerInformationIO::overwriteLayerInformation( _lineages.at(_sDataId)->getFullPath(),
                                                    _cellLayers->getCellLayersPerTimeStep() );
}

// ---------------------------------------------------------------------

void SVLRBrowser::bindSlider( const std::string &sliderName )
{
  _sliderCallback = new SSliderCallback();
  _connections.push_back( SKernel::get()->onSliderUpdate( sliderName,
                                                          boost::bind( &SSliderCallback::updateStep,
                                                                       _sliderCallback, _1 ) ) );
}

// ---------------------------------------------------------------------

void SVLRBrowser::onSliderUpdate( int val )
{
  SKernel::get()->requestRedraw();
}

// ---------------------------------------------------------------------

boost::signals2::connection SVLRBrowser::onHighlightChanged( const slotHighlightUpdate& slot )
{
  return( _signalHighlightUpdate.connect( slot ) );
}

// ---------------------------------------------------------------------

void SVLRBrowser::highlightCell( osg::Node* node, const NodeType::NodeType nT )
{
  if( !_firstClick )
    this->selectNode( node, false, nT );

  this->selectNode( node, true, nT );

  if( _firstClick )
    _firstClick = false;

  _signalHighlightCell( node );
}

// ---------------------------------------------------------------------

void SVLRBrowser::selectNode( osg::Node *node,
                              const bool select,
                              const NodeType::NodeType nT )
{
  osg::Vec4 color;

  // selection color
  if( select )
    color = _selectionColor;
  else
    color = _lastClickedColor;

  SHighlightCellVisitor hCV( color, nT );
  if( select )
    node->accept( hCV );
  else
  {
    osg::Geode *geode = new osg::Geode;
    if( nT == NodeType::NODE2D )
      geode->addDrawable( _lastClicked2DNode );
    else
      geode->addDrawable( _lastClicked3DNode );

    geode->accept( hCV );
  }

  // store last color
  if( select )
  {
    _lastClickedColor = hCV.getOriginColor();

    if( nT == NodeType::NODE2D )
      _lastClicked2DNode = dynamic_cast<SSelectableEllipse*>( node->asGeode()->getDrawable(0) );
    else
      _lastClicked3DNode = dynamic_cast<SSelectableSphere*>( node->asGeode()->getDrawable(0) );
  }

  SSelectableEllipse *ellipse = 0;
  SSelectableSphere *sphere = 0;

  if( nT == NodeType::NODE2D )
  {
    if( select )
      ellipse = dynamic_cast<SSelectableEllipse*>( node->asGeode()->getDrawable(0) );
    else
      ellipse = _lastClicked2DNode;
  }
  else
  {
    if( select )
      sphere = dynamic_cast<SSelectableSphere*>( node->asGeode()->getDrawable(0) );
    else
      sphere = _lastClicked3DNode;
  }

  // get the node type which should be updated
  // in the other window
  NodeType::NodeType type;

  if( nT == NodeType::NODE2D )
    type = NodeType::NODE3D;
  else
    type = NodeType::NODE2D;

  SHighlightCellVisitor hCV2( color, type );
  osg::Geode *geode = new osg::Geode;

  // this is a special case which can occur when the follwing cases are true:
  // (a) we deselect the previous node
  // (b) the current node type is a 2D node in the lineage
  // (c) if it is a 2D node and if the ellipse pointer is zero
  // (d) if we render no move nodes
  // if we render no move nodes then it could be that the previous 2D node
  // was not rendered and then the ellipse pointer is zero, however the
  // 3D node has to be deselected which is done with the call below
  bool no2DNode = !select && !ellipse && nT == NodeType::NODE2D;
  // deselect previous 3D node
  if( no2DNode )
  {
    cellInfo cI = std::make_pair( _clickedCellId, _clickedTimeStep );
    geode->addDrawable( _geomNodeMaps->get3DGeometry( cI ) );
  }

  if( ellipse || sphere )
  {
    std::string name;

    if( nT == NodeType::NODE2D )
      name = ellipse->getName();
    else
      name = sphere->getName();

    std::vector<std::string> idT;

    boost::split( idT, name, boost::is_any_of(" ") );
    int cellId = boost::lexical_cast<int>( idT[1] );
    int timeStep;
    int dataId = _sDataId;
    if( nT == NodeType::NODE2D )
      timeStep = boost::lexical_cast<int>( idT[2] );
    else
      timeStep = SKernel::get()->getSliderPosition();

    if( select )
    {
      _clickedTimeStep = timeStep;
      _clickedCellId = cellId;

      // get tree id and layer value
      int layer;

      // get node pointer of selected node
      SLineageTree *selectedNode = this->getTreeIdOfSelectedCell( dataId, layer );

      osg::Vec3 pos( selectedNode->getX(), selectedNode->getY(), selectedNode->getZ() );
      pos = pos * _vlrRoot->getMatrix();

      // find cell file info
      std::map<int,int>::const_iterator cFIter =
          _cellFiles.at(_sDataId).find(selectedNode->treeId);

      // find division type info
      std::map<std::size_t,int>::const_iterator cDTIter =
          _divisionScheme.at(_sDataId).at(timeStep-1).find(cellId);

      // if the division type does not exist then assign a
      // value of 3 to indicate its non-existance
      std::size_t divType;
      if( cDTIter != _divisionScheme.at(_sDataId).at(timeStep-1).end() )
        divType = cDTIter->second;
      else
        divType = 3;

      // only with the drawing of the cell visualization
      // we also generate the cell layers
      // TODO
      if( _draw3DCells && false )
      {
        // get cell layer history of selected cell
        cellHistoryMap::const_iterator lHIter =
            _cellLayers->getCellHistories()->at(timeStep-1).find( selectedNode );

        // if found then set everything
        if( lHIter != _cellLayers->getCellHistories()->at(timeStep-1).end() &&
            cFIter != _cellFiles.at(_sDataId).end() )
        {
          cellHistory cellHist = lHIter->second;
          this->setCellFeatures( cellId, timeStep, selectedNode->treeId, layer, cellHist,
                                 pos[0], pos[1], pos[2], cFIter->second, divType );
        }
      }
      // else we just set the time step, cell id, tree id and layer information
      else
      {
        if( cFIter != _cellFiles.at(_sDataId).end() )
        {
          cellHistory cellHist;
          this->setCellFeatures( cellId, timeStep, selectedNode->treeId, layer, cellHist,
                                 pos[0], pos[1], pos[2], cFIter->second, divType );
        }
      }

    }
    else
    {
      timeStep = _clickedTimeStep;
      cellId = _clickedCellId;
    }

    cellInfo cI = std::make_pair( cellId, timeStep );

    if( type == NodeType::NODE2D )
    {
      geode->addDrawable( _geomNodeMaps->get2DGeometry( cI ) );
      // store last 2D cell
      if( select )
        _lastClicked2DNode = _geomNodeMaps->get2DGeometry( cI );
    }
    else
    {
      geode->addDrawable( _geomNodeMaps->get3DGeometry( cI ) );
      // store last 3D cell
      if( select )
        _lastClicked3DNode = _geomNodeMaps->get3DGeometry( cI );
    }

    // emit signal to set new highlight position of arrows
    cellDataInfo cDI;
    cDI.push_back( cellId );
    cDI.push_back( timeStep );
    cDI.push_back( dataId );
    _signalHighlightUpdate( cDI );
  }

  geode->accept( hCV2 );

//  SSVGRenderer svg( true );
//  svg.startImage( "/tmp/lineageTrees.svg" );
//  svg.renderOSGNode( _lineageLayerGroup, 1., 1., 0.25 );
}

// ---------------------------------------------------------------------

void SVLRBrowser::updateLayerGroups( osg::Node* node,
                                     const int newLayer,
                                     const int oldLayer )
{
  osg::ref_ptr<osg::Group> layerGroup;

  // get layer group
  for( std::size_t c = 0;c<_vlrRoot->asMatrixTransform()->getNumChildren();c++ )
  {
    if( _vlrRoot->getChild(c)->getName() == "Cell Layers" )
      layerGroup = _vlrRoot->getChild(c)->asGroup();
  }

  if( layerGroup->getName() == "" )
    std::cout << "No cell layers were found!" << std::endl;

  SSelectableSphere *sphere = dynamic_cast<SSelectableSphere*>( node->asGeode()->getDrawable(0) );

  // for each layer group
  for( std::size_t c = 0;c<layerGroup->getChild(oldLayer)->asGroup()->getNumChildren();c++ )
  {
    SSelectableSphere *curSphere = dynamic_cast<SSelectableSphere*>( layerGroup->getChild(oldLayer)->asGroup()->
                                                                     getChild(c)->asTransform()->asMatrixTransform()->
                                                                     getChild(0)->asGeode()->getDrawable(0) );

    if( curSphere )
    {
      // update layer assignment of 3D nodes in groups
      if( curSphere->getName() == sphere->getName() )
      {
        layerGroup->getChild( newLayer )->asGroup()->addChild( layerGroup->getChild(oldLayer)->asGroup()->getChild(c) );
        layerGroup->getChild( oldLayer )->asGroup()->removeChild( c, 1 );
      }
    }
  }
}

// ---------------------------------------------------------------------

void SVLRBrowser::assignLayer( osg::Node* node, const int layer,
                               const NodeType::NodeType nT )
{
  std::cout << "assign layer " << layer << std::endl;

  std::string name = node->asGeode()->getDrawable(0)->getName();
  std::vector<std::string> idT;

  boost::split( idT, name, boost::is_any_of(" ") );
  int cellId = boost::lexical_cast<int>( idT[1] );
  int timeStep;
  if( nT == NodeType::NODE2D )
    timeStep = boost::lexical_cast<int>( idT[2] );
  else
    timeStep = SKernel::get()->getSliderPosition();

  std::cout << "assign id: " << cellId << " t: " << timeStep << std::endl;

  unsigned int oldLayerValue;

  // assign new layer value to selected node if necessary
  bool needsUpdate = _cellLayers->assignNewLayerValue( layer, cellId, timeStep, oldLayerValue );

  if( !needsUpdate )
    return;

  // re-assign 3D nodes to layers
  this->updateLayerGroups( node, layer, oldLayerValue );

  // clear all data after selected time step
  _cellLayers->clearLayerInformation( timeStep );
  // re-determine layers
  SInteractiveSurface *is = 0;
  SInteractiveNormals *in = 0;

  for( std::size_t c = 0;c<_vlrRoot->asMatrixTransform()->getNumChildren();c++ )
  {
    if( !in )
      in = dynamic_cast<SInteractiveNormals*>( _vlrRoot->getChild(c) );
    if( !is )
      is = dynamic_cast<SInteractiveSurface*>( _vlrRoot->getChild(c) );

    if( is && in )
      break;
  }

  if( is && in )
    _cellLayers->updateLayers( timeStep, is, in );
  else
    std::cout << "cast failed" << std::endl;

  // update 3D nodes
  this->updateLayerColorsForCells();
  // update 2D nodes
  if( _drawLineages )
  {
    SColorCellVisitor ccv2D( *(_cellLayers->getCellLayersPerTimeStep()), _layerColors );
    _lineageLayerGroup->accept( ccv2D );
  }

  SKernel::get()->setProgress( 1., "Done" );

  _signalLayerAssignment( node );
}

// ---------------------------------------------------------------------

void SVLRBrowser::setDivisionType( osg::Node* node, const int divisionType,
                                   const NodeType::NodeType nT )
{
  if( divisionType != 0 && divisionType != 2 )
    return;

  std::string name = node->asGeode()->getDrawable(0)->getName();
  std::vector<std::string> idT;

  boost::split( idT, name, boost::is_any_of(" ") );
  int cellId = boost::lexical_cast<int>( idT[1] );
  int timeStep;
  if( nT == NodeType::NODE2D )
    timeStep = boost::lexical_cast<int>( idT[2] );
  else
    timeStep = SKernel::get()->getSliderPosition();

  // get current division type of cell and check if update is required
  std::map<std::size_t,int>::iterator findIter = _divisionScheme.at(_sDataId).at(timeStep-1).find( cellId );

  if( findIter != _divisionScheme.at(_sDataId).at(timeStep-1).end() )
  {
    if( divisionType == findIter->second )
    {
      std::cout << "Selected division type is already set." << std::endl;
      return;
    }
    else
      findIter->second = divisionType;
  }
  // else the selected cell is not a dividing cell
  else
  {
    std::cout << "The selected cell is not a dividing cell." << std::endl;
    return;
  }

  this->highlightCell( node, nT );

  if( divisionType == 0 )
    std::cout << "Set anticlinal division type for id: " << cellId << " and t: " << timeStep << std::endl;
  else
    std::cout << "Set radial division type for id: " << cellId << " and t: " << timeStep << std::endl;
}

// ---------------------------------------------------------------------

SLineageTree* SVLRBrowser::getTreeIdOfSelectedCell( int &dataid, int &layer )
{
  // loop over all data sets
  for( std::size_t d = 0; d < _numData; d++ )
  {
    for( SAdvancedLineageTreeCollection::const_iterator l = _lineages.at(d)->begin();
         l != _lineages.at(d)->end(); ++l )
    {
      // get lineage id
      int id = l->first;

      SLineageTree *tree = (*_lineages.at(d))[id];

      for( SLineageTree::iterator iter = tree->begin(); iter != tree->end(); ++iter )
      {
        if( iter->cellId == _clickedCellId &&
            iter->timeStep == _clickedTimeStep )
        {
          NodeFeatureInfo cellLayers = _lineages.at(d)->getCellLayers();
          layer = cellLayers.at(iter->timeStep-1).at(iter->cellId);
          dataid = d;
          return *iter;
        }
      }
    }
  }

  layer = dataid = 0;
  return 0;
}

// ---------------------------------------------------------------------

void SVLRBrowser::updateLayerColorsForCells()
{
  for( std::vector< std::pair<SInteractiveCell*,const SLineageTree*> >::iterator
       iter = _interactiveCellVector.begin(); iter != _interactiveCellVector.end(); ++iter )
  {
    SSelectableSphere *sphere = dynamic_cast<SSelectableSphere*>( iter->first->asTransform()->asMatrixTransform()->
                                                                  getChild(0)->asGeode()->getDrawable(0));

    if( sphere )
    {
      //unsigned int layerValue = _cellLayers->determineLayerValue( iter->second );
      unsigned int layerValue = _layerValues.at( _sDataId ).at( iter->second->timeStep-1 ).at( iter->second->cellId );
      osg::Vec4 color = (*_layerColors)[layerValue%_layerColors->size()];

      SHighlightCellVisitor hCV( color, NodeType::NODE3D );
      osg::Geode *geode = new osg::Geode;
      geode->addDrawable( sphere );
      geode->accept( hCV );

      // if the selected node to be assigned to a new layer
      // was the laste node highlighted then we store its
      // new layer color as its previous one
      if( iter->second->cellId == _clickedCellId &&
          iter->second->timeStep == _clickedTimeStep )
        _lastClickedColor = color;
    }
  }
}

// ---------------------------------------------------------------------

boost::signals2::connection SVLRBrowser::onHighlightCell( const slotHighlight& slot )
{
  return( _signalHighlightCell.connect( slot ) );
}

// ---------------------------------------------------------------------

boost::signals2::connection SVLRBrowser::onLayerAssignment( const slotLayerAssignment& slot )
{
  return( _signalLayerAssignment.connect( slot ) );
}

// ---------------------------------------------------------------------

void SVLRBrowser::renderSingleTimeSteps()
{
  int minCellFile = 10;
  int maxCellFile = -10;

  CellWalls cellWalls;
  SInteractiveSurface *is = new SInteractiveSurface();
  std::vector <osg::ref_ptr<osg::Group> > wallGroup;

  if( _renderCellWalls && _loadModel )
  {
    osg::ref_ptr<osg::Group> cellWallGroup = new osg::Group;
    is->setName( "Cell Walls" );
    is->addUpdateCallback( _sliderCallback );
    is->addChild( cellWallGroup );

    // read cell wall information
    SLayerInformationIO::readCellWallInformation( _cellWallDirectory, cellWalls );
    for( std::size_t t=0; t<cellWalls.size(); t++ )
    {
      osg::ref_ptr<osg::Group> group = new osg::Group;
      wallGroup.push_back( group );
    }

    _vlrRoot->addChild( is );
  }

  // group for cell layers/files or lineage trees
  osg::ref_ptr<osg::Group> layerGroup = new osg::Group;
  layerGroup->setName( "Cells" );
  boost::shared_ptr<NodeFeatureInfo> cellLayers;

  // data type holding life duration information
  std::map<std::size_t, std::size_t> lifeDuration;

  // layer based coloring
  if( _render3DType == 0 )
  {
    cellLayers = _cellLayers->getCellsPerTimeStepInLayer();

    for( std::size_t i = 0; i < _cellLayers->getMaxNumLayers();i++ )
    {
      osg::ref_ptr<osg::Group> subLayerGroup = new osg::Group;
      subLayerGroup->setName( "Layer " + boost::lexical_cast<std::string>( i ) );
      subLayerGroup->setUpdateCallback( _sliderCallback );
      layerGroup->addChild( subLayerGroup );
    }
  }
  // cell files
  else if( _render3DType == 1 )
  {
    // get max number of cell files for the last time step
    for( std::map<int,int>::const_iterator mapIter =
         _cellFiles.at(_sDataId).begin(); mapIter !=
         _cellFiles.at(_sDataId).end(); mapIter++ )
    {
      if( mapIter->second < minCellFile )
        minCellFile = mapIter->second;

      if( mapIter->second > maxCellFile )
        maxCellFile = mapIter->second;
    }

    for( int i = minCellFile; i <= maxCellFile; i++ )
    {
      osg::ref_ptr<osg::Group> subLayerGroup = new osg::Group;
      subLayerGroup->setName( "Cell File " + boost::lexical_cast<std::string>( i ) );
      subLayerGroup->setUpdateCallback( _sliderCallback );
      layerGroup->addChild( subLayerGroup );
    }
  }
  // lineage tree based coloring
  else if( _render3DType == 2 )
  {
    for( std::size_t i = 0; i < _lineages.at(_sDataId)->size();i++ )
    {
      osg::ref_ptr<osg::Group> subLayerGroup = new osg::Group;
      subLayerGroup->setName( "Lineage " + boost::lexical_cast<std::string>( i ) );
      subLayerGroup->setUpdateCallback( _sliderCallback );
      layerGroup->addChild( subLayerGroup );
    }
  }
  // layering coloring type for which both children are updated
  else if( _render3DType == 3 )
  {
    std::cout << "maxLayers: " << _maxLayers << std::endl;

    for( std::size_t i = 0; i < _maxLayers + 1;i++ )
    {
      osg::ref_ptr<osg::Group> subLayerGroup = new osg::Group;
      subLayerGroup->setName( "Layer " + boost::lexical_cast<std::string>( i ) );
      subLayerGroup->setUpdateCallback( _sliderCallback );
      layerGroup->addChild( subLayerGroup );
    }
  }
  // single color
  else if( _render3DType == 4 )
  {
    osg::ref_ptr<osg::Group> subLayerGroup = new osg::Group;
    subLayerGroup->setUpdateCallback( _sliderCallback );
    layerGroup->addChild( subLayerGroup );
  }
  // life duration
  else if( _render3DType == 5 )
  {
    std::size_t min, max;
    SSimilarityMeasureUtil::determineCellLifeDuration( _lineages.at(_sDataId),
                                                       _centresOfMass.at(_sDataId),
                                                       _cellFiles.at(_sDataId),
                                                       _rotMatrices.at(_sDataId),
                                                       lifeDuration,
                                                       min, max );

    std::cout << "min: " << min << " max: " << max << std::endl;
    osg::ref_ptr<osg::Group> subLayerGroup = new osg::Group;
    subLayerGroup->setUpdateCallback( _sliderCallback );
    layerGroup->addChild( subLayerGroup );
  }
  // quiescent center
  else if( _render3DType == 6 )
  {
    osg::ref_ptr<osg::Group> subLayerGroup = new osg::Group;
    subLayerGroup->setUpdateCallback( _sliderCallback );
    layerGroup->addChild( subLayerGroup );
  }

  _vlrRoot->addChild( layerGroup );

  int numTrees = 0;

  for( auto l = _lineages.at(_sDataId)->begin(); l != _lineages.at(_sDataId)->end(); ++l )
  {
    SLineageTree *tree = l->second;

    if( _onlyMasterFiles )
    {
      if( _cellFiles.at(_sDataId).at( tree->treeId ) != 0 )
        continue;
    }

    SLineageTree *currCell;
    std::stack<SLineageTree*> stack;
    stack.push( tree->getRoot() );

    while( !stack.empty() )
    {
      currCell = stack.top();
      stack.pop();

      int layerOrCellFileValue;
      osg::Vec4 color;

      if( _render3DType == 0 )
      {
        if( !_loadModel )
        {
          if( _cellLayers->getMaxNumLayers() != 0 )
            layerOrCellFileValue = _cellLayers->determineLayerValue( currCell );
          else
            layerOrCellFileValue = cellLayers->at(currCell->timeStep-1).at(currCell->cellId);

          color = osg::Vec4( (*_layerColors)[layerOrCellFileValue%_layerColors->size()] );
        }
        else
        {
          NodeFeatureInfo lay = _lineages.at(_sDataId)->getCellLayers();
          layerOrCellFileValue = lay.at(currCell->timeStep-1).at(currCell->cellId);
          color = osg::Vec4( (*_layerColors)[(layerOrCellFileValue-1)%_layerColors->size()] );
        }
      }
      // layer coloring by default else the cells are colored based on their cell file info
      else if( _render3DType == 1 )
      {
        layerOrCellFileValue = _cellFiles.at(_sDataId).at( currCell->treeId );
        color = _cellFileColors[layerOrCellFileValue];
        // substract min value such the values always starts with zero in order
        // to use it as a child index for the layerGroup
        layerOrCellFileValue -= minCellFile;
      }
      // coloring based on lineage tree assignment
      else if( _render3DType == 2 )
      {
        layerOrCellFileValue = _treeIdMap[currCell->treeId];
        color = _cellLineageColors[layerOrCellFileValue];
        layerOrCellFileValue--;
      }
      else if( _render3DType == 3 )
      {
        layerOrCellFileValue = _layerValues.at(_sDataId).at(currCell->timeStep-1).at(currCell->cellId) - 1;
        color = osg::Vec4( (*_layerColors)[(layerOrCellFileValue)%_layerColors->size()] );
      }
      else if( _render3DType == 4 )
      {
        layerOrCellFileValue = 0;
        color = osg::Vec4( 0., 1., 0., 1. );
      }
      else if( _render3DType == 5 )
      {
        layerOrCellFileValue = 0;
        auto iter = lifeDuration.find( currCell->cellId );
        if( iter != lifeDuration.end() )
        {
          std::size_t life = iter->second;
          color = _lifeDurationColorMap.mapToColor( life );
        }
        else
          color = osg::Vec4( 1., 1., 1., 0.4 );
      }
      else if( _render3DType == 6 )
      {
        layerOrCellFileValue = 0;
        int cellFi = _cellFiles.at(_sDataId).at( currCell->treeId );

        if( cellFi == 1 && currCell->treeId == 7 &&
            currCell->timeStep-1 > 200 )
          color = osg::Vec4( 1., 0., 0., 1. );
        else
          color = osg::Vec4( 1., 1., 1., 0.2 );
      }

      SInteractiveCell* cell = new SInteractiveCell( _cellRadius,
                                                     currCell->cellId,
                                                     currCell->timeStep,
                                                     color );

      layerGroup->getChild( layerOrCellFileValue )->asGroup()->addChild( cell );

      // add this interactive cell to the vector for later color update
      _interactiveCellVector.push_back( std::make_pair( cell, currCell ) );

      SActions::signal_node_t signal2 = SActions::get( *cell ).getNodeSignal( "Highlight cell" );
      _connections.push_back( signal2->connect(
                                boost::bind( &SVLRBrowser::highlightCell, this, _1, NodeType::NODE3D ) ) );

      // TODO: set max layers and refactor assignLayer method
      if( _lineageColorType == 0 && false )
      {
        for( std::size_t i=0; i < _cellLayers->getMaxNumLayers();i++ )
        {
          SActions::signal_node_t signal1 = SActions::get( *cell ).getNodeSignal(
                "Assign layer " + boost::lexical_cast<std::string>(i) );
          _connections.push_back( signal1->connect(
                                    boost::bind( &SVLRBrowser::assignLayer, this, _1, i, NodeType::NODE3D ) ) );
        }
      }

      while( currCell != 0 )
      {
        SSelectableSphere *sphere = static_cast<SSelectableSphere*>( cell->asMatrixTransform()->getChild(0)->
                                                                     asGeode()->getDrawable( 0 ) );

        cellInfo cI = std::make_pair( currCell->cellId, currCell->timeStep );

        _geomNodeMaps->insertCellInfoGeometry( cI, sphere );

        cell->addPosition( currCell->timeStep,
                           osg::Vec3( currCell->getX(), currCell->getY(),
                                      currCell->getZ() ) );

        // render cell walls for each cell and time step
        if( _renderCellWalls && _loadModel )
        {
          std::map<std::size_t,std::vector<osg::Vec3> >::const_iterator iter =
              cellWalls.at(currCell->timeStep-1).find( currCell->cellId );

          if( iter != cellWalls.at(currCell->timeStep-1).end() )
          {
            osg::Vec3 pos1 = iter->second.front();
            for( std::size_t n=1; n<iter->second.size(); n++ )
            {
              osg::Vec3 pos2 = iter->second.at(n);

              SSimilarityMeasureGraphics::addCylinderBetweenPoints( pos1, pos2,
                                                                _cellWallWidth, color,
                                                                wallGroup.at(currCell->timeStep-1) );

              pos1 = pos2;
            }

            // render last cell wall joint
            SSimilarityMeasureGraphics::addCylinderBetweenPoints( iter->second.back(), iter->second.front(),
                                                              _cellWallWidth, color,
                                                              wallGroup.at(currCell->timeStep-1));
          }
        }

        if( currCell->children.size() == 2 )
        {
          // add entry for setting division type info only if this current cell is a dividing cell
          // only two allowed division types: anticlinal and radial
          // that do not change the layering
          SActions::signal_node_t signal3 = SActions::get( *cell ).getNodeSignal(
                "Set to anticlinal division" );
          _connections.push_back( signal3->connect(
                                    boost::bind( &SVLRBrowser::setDivisionType, this, _1, 0, NodeType::NODE3D ) ) );

          SActions::signal_node_t signal4 = SActions::get( *cell ).getNodeSignal(
                "Set to radial division" );
          _connections.push_back( signal4->connect(
                                    boost::bind( &SVLRBrowser::setDivisionType, this, _1, 2, NodeType::NODE3D ) ) );

          stack.push( currCell->children[0] );
          stack.push( currCell->children[1] );
          break;
        }
        else if( currCell->children.size() == 1 )
          currCell = currCell->children[0];
        else
          currCell = 0;
      }

      cell->update( 1 );
    }

    ++numTrees;
  }

  if( _renderCellWalls && _loadModel )
  {
    for( std::size_t t=0; t<wallGroup.size(); t++ )
      is->addTriangulation( t+1, 0, wallGroup.at(t) );
  }
}

// ---------------------------------------------------------------------
