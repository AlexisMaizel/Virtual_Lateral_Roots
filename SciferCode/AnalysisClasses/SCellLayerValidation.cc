#include "SCellLayerValidation.hh"

/**
  @file   SCellLayerValidation.cc
  @brief  Contains class for validating layer results
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include <kernel/SKernel.hh>

#include <boost/numeric/conversion/bounds.hpp>
#include <boost/limits.hpp>

// ---------------------------------------------------------------------

SCellLayerValidation::SCellLayerValidation( const std::string &MIPFilename,
                                            const int projDir,
                                            const bool useSlicing,
                                            boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                                            osg::ref_ptr<SSliderCallback> sliderCallback,
                                            osg::ref_ptr<osg::MatrixTransform> addToThisGroup )
  : _lineages( lineages ),
    _projDir( projDir ),
    _MIPFilename( MIPFilename ),
    _useSlicing( useSlicing ),
    _vlrGroup( addToThisGroup )
{
  // generate bounding box information of data in order
  // to set the clipping planes and MIPs correctly
  this->generateBBox();

  if( _MIPFilename != "" )
  {
    // initialize and set MIPs of the raw data set
    _MIPValidation = boost::shared_ptr<SMIPValidation>(
          new SMIPValidation( _MIPFilename, _bbox, sliderCallback, _lineages ) );
    {
      // wait for the MIP init
      boost::thread t(  boost::bind( &SCellLayerValidation::renderMIPs, this, _projDir ) );
      t.join();
    }
  }

  if( _useSlicing )
  {
    // initialize interactive clipping plane
    {
      // wait for the clipping init
      boost::thread t(  boost::bind( &SCellLayerValidation::renderClippingPlane, this, _projDir ) );
      t.join();
    }
  }
}

// ---------------------------------------------------------------------

SCellLayerValidation::~SCellLayerValidation()
{
  delete _bbox;
}

// ---------------------------------------------------------------------

void SCellLayerValidation::generateBBox()
{
  // check and set the minimum and maximum bbox values
  // TODO: at the moment this does only fit to the MIPs if no rotation is applied
  double min = boost::numeric::bounds<double>::lowest();
  double max = boost::numeric::bounds<double>::highest();
  _bbox = new osg::BoundingBox;
  _bbox->set( max, max, max, min, min, min );

  for( SAdvancedLineageTreeCollection::const_iterator l = _lineages->begin();
       l != _lineages->end(); ++l )
  {
    SLineageTree *tree = (*_lineages)[l->first];

    for( SLineageTree::iterator treeIter = tree->begin();
         treeIter != tree->end(); ++treeIter )
    {
      osg::Vec3 pos( treeIter->getX(), treeIter->getY(), treeIter->getZ() );
      // apply rotation to position
      pos = pos * _vlrGroup->getMatrix();

      if( pos[0] < _bbox->xMin() )
        _bbox->xMin() = pos[0];
      if( pos[0] > _bbox->xMax() )
        _bbox->xMax() = pos[0];

      if( pos[1] < _bbox->yMin() )
        _bbox->yMin() = pos[1];
      if( pos[1] > _bbox->yMax() )
        _bbox->yMax() = pos[1];

      if( pos[2] < _bbox->zMin() )
        _bbox->zMin() = pos[2];
      if( pos[2] > _bbox->zMax() )
        _bbox->zMax() = pos[2];
    }
  }
}

// ---------------------------------------------------------------------

void SCellLayerValidation::update( const int index )
{
  if( index == _projDir )
  {
    std::cout << "Projection already set." << std::endl;
    return;
  }
  else
  {
    std::string dir;
    switch( index )
    {
    case 0: dir = "X"; break;
    case 1: dir = "Y"; break;
    case 2: dir = "Z"; break;
    case 3: dir = "All"; break;
    }

    std::cout << "Projection direction set to " << dir << "." << std::endl;
    _projDir = index;
  }

  // remove all previous osg nodes related to MIP and slicing
  this->removeProjections();

  // redraw MIP from new projection direction in a waiting thread
  if( _MIPFilename != "" )
  {
    boost::thread t(  boost::bind( &SCellLayerValidation::renderMIPs, this, _projDir ) );
    t.join();
  }

  // redraw clipping plane within a waiting thread
  if( _useSlicing )
  {
    boost::thread t(  boost::bind( &SCellLayerValidation::renderClippingPlane, this, _projDir ) );
    t.join();
  }
}

// ---------------------------------------------------------------------

void SCellLayerValidation::removeProjections()
{
  // here we remove all stored osg nodes which are
  // related to the MIP or slicing if it is enabled
  // Note that we have to loop over them for three
  // times in order to access the correct child
  // in the root group

  // remove clipping plane
  if( _useSlicing )
  {
    for( std::size_t p = 0; p < 3; p++ )
    {
      std::string subName = "Slicing";
      for( std::size_t c=0;c< _vlrGroup->getNumChildren();c++ )
      {
        if( _vlrGroup->getChild(c)->getName().find( subName )
            != std::string::npos )
        {
          // before removing the node from the root node
          // we first have to deactivate the clipping for it
          // since this one will not be removed otherwise
          SClippingPlane *clipPlane =
              dynamic_cast<SClippingPlane*>( _vlrGroup->getChild(c) );

          if( clipPlane )
            clipPlane->disableClipping( _vlrGroup );

          _vlrGroup->removeChild( c, 1 );
          break;
        }
      }
    }
  }

  // remove MIP texture
  if( _MIPFilename != "" )
  {
    for( std::size_t p = 0; p < 3; p++ )
    {
      std::string subName = "MIPs";
      for( std::size_t c=0;c< _vlrGroup->getNumChildren();c++ )
      {
        if( _vlrGroup->getChild(c)->getName().find( subName )
            != std::string::npos )
        {
          _vlrGroup->removeChild( c, 1 );
          break;
        }
      }
    }
  }
}

// ---------------------------------------------------------------------

void SCellLayerValidation::initializeClippingPlane( const int projDir )
{
  osg::ref_ptr<SClippingPlane> plane = new SClippingPlane( projDir, _vlrGroup,
                                                           _bbox, _lineages );
  _vlrGroup->addChild( plane );
}

// ---------------------------------------------------------------------

void SCellLayerValidation::renderClippingPlane( const int projDir )
{
  // single projection
  if( projDir != 3 )
    this->initializeClippingPlane( projDir );
  // all projections
  else
  {
    // for all projections
    for( std::size_t p = 0; p < 3; p++ )
    {
      boost::thread t(  boost::bind( &SCellLayerValidation::initializeClippingPlane, this, p ) );
      t.join();
    }

    // avoid the clipping planes to be clipped by each other
    std::string subName = "Slicing";
    std::vector<SClippingPlane*> planes;
    for( std::size_t c = 0; c < _vlrGroup->getNumChildren(); c++ )
    {
      if( _vlrGroup->getChild(c)->getName().find( subName ) != std::string::npos )
      {
        SClippingPlane *clipPlane =
            dynamic_cast<SClippingPlane*>( _vlrGroup->getChild(c) );

        if( clipPlane )
          planes.push_back( clipPlane );
      }
    }

    for( std::size_t p = 0; p < planes.size(); p++ )
    {
      if( p == 0 )
      {
        planes.at(p)->disableRectClipping( planes.at(1)->getRenderedPlane() );
        planes.at(p)->disableRectClipping( planes.at(2)->getRenderedPlane() );
      }
      else if( p == 1 )
      {
        planes.at(p)->disableRectClipping( planes.at(0)->getRenderedPlane() );
        planes.at(p)->disableRectClipping( planes.at(2)->getRenderedPlane() );
      }
      else if( p == 2 )
      {
        planes.at(p)->disableRectClipping( planes.at(0)->getRenderedPlane() );
        planes.at(p)->disableRectClipping( planes.at(1)->getRenderedPlane() );
      }
    }
  }
}

// ---------------------------------------------------------------------

void SCellLayerValidation::renderMIPs( const int projDir )
{
  // render single MIPs for the chosen projection direction
  if( projDir != 3 )
    _MIPValidation->mapRawData( projDir, _vlrGroup );
  // else render all projections
  else
  {
    for( std::size_t i = 0; i < 3; i++ )
      _MIPValidation->mapRawData( i, _vlrGroup );
  }
}

// ---------------------------------------------------------------------
