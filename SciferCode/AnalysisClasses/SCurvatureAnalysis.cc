#include "SCurvatureAnalysis.hh"

/**
  @file   SCurvatureAnalysis.cc
  @brief  Contains class for analyzing the curvature of VLR data
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SInteractiveSurface.hh"

#include <boost/foreach.hpp>

// ---------------------------------------------------------------------

SCurvatureAnalysis::SCurvatureAnalysis( osg::ref_ptr<SSliderCallback> sliderCallback,
                                        boost::shared_ptr<cellTimeVector> cTS,
                                        const curvatureType cType,
                                        std::vector<std::pair<double,double> > &minMaxCurvatures,
                                        const int shapeType,
                                        osg::ref_ptr<osg::MatrixTransform> addToThisGroup )
{
  SInteractiveSurface *is = new SInteractiveSurface();
  osg::ref_ptr<osg::Group> volumeGroup = new osg::Group;
  is->setName( "Curvature" );
  is->addUpdateCallback( sliderCallback );
  is->addChild( volumeGroup );

  for( std::size_t t = 0; t < cTS->size(); t++ )
  {
    osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array;
    std::map<osg::Vec3,std::pair<std::size_t,std::size_t> > cells;

    BOOST_FOREACH( const SLineageTree *tree, cTS->at(t) )
    {
      points->push_back( osg::Vec3( tree->getX(), tree->getY(), tree->getZ() ) );
      cells[ osg::Vec3( tree->getX(), tree->getY(), tree->getZ() ) ] =
          std::make_pair( tree->timeStep, tree->cellId );
    }

    // convex hull
    if( shapeType == 0 )
    {
      boost::shared_ptr<SConvexHull> tri = boost::shared_ptr<SConvexHull>(
            new SConvexHull( points, t, cType, minMaxCurvatures, cells ) );

      is->addTriangulation( t+1, 0, tri->getTriangles() );

      //tri->exportCurvatureValuesOfCells();
    }
    // alpha shape
    else if( shapeType == 1 )
    {
      boost::shared_ptr<SAlphaShape> tri = boost::shared_ptr<SAlphaShape>(
            new SAlphaShape( points, t, cType, minMaxCurvatures, cells ) );

      is->addTriangulation( t+1, 0, tri->getTriangles() );

      //tri->exportCurvatureValuesOfCells();
    }
  }

  addToThisGroup->addChild( is );
}

// ---------------------------------------------------------------------
