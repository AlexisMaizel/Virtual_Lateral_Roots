#include "SCorrelationAnalysis.hh"

#include "SVLRDataAnalysisIO.hh"

// ---------------------------------------------------------------------

SCorrelationAnalysis::SCorrelationAnalysis( std::vector< boost::shared_ptr<SAdvancedLineageTreeCollection> > lineages,
                                            const std::size_t dim )
  : _lineages( lineages ),
    _dim( dim ),
    _numPoints( 0 ),
    _dataType( 0 ),
    _loadData( false )
{
}

// ---------------------------------------------------------------------

void SCorrelationAnalysis::readDivisionData( const std::string &fileName,
                                             const std::size_t dataId )
{
  std::ifstream in( fileName.c_str(), std::ofstream::in );

  if( !in.is_open() )
  {
    THROW_EXCEPTION( VInvalidFilenameException, "Unable to open " + fileName );
    return;
  }

  std::string line;
  getline( in, line );

  while( in.good() )
  {
    getline( in, line );
    if( line == "" )
      break;

    divData dd;

    std::string dummy;
    std::size_t id, t, divType;
    double VAngle, distVA, primAngle, distPrim, surAngle, distSur;
    std::stringstream lineStream( line );
    lineStream >> id >> t;

    for( std::size_t e = 0; e < 4; e++ )
      lineStream >> dummy;

    lineStream >> VAngle >> distVA >> primAngle >> distPrim >> surAngle >> distSur;

    dd.VLRAngle = VAngle;
    dd.distDivToCentre = distVA;
    dd.primAngle = primAngle;
    dd.distDivToMainRoot = distPrim;
    dd.surfaceAngle = surAngle;
    dd.distDivToSurface = distSur;

    for( std::size_t e = 0; e < 10; e++ )
      lineStream >> dummy;

    if( dummy.at(0) == '\"' )
    {
      while( *dummy.rbegin() != '\"' )
        lineStream >> dummy;
    }

    lineStream >> divType;
    dd.divType = divType;
    _readDivisionData.at(dataId).insert( std::make_pair( std::make_pair( id, t-1 ), dd ) );
  }

  in.close();
}

// ---------------------------------------------------------------------

void SCorrelationAnalysis::generateDivisionData( const std::size_t divType,
                                                 const std::size_t maxTimeSteps )
{
  _dataType = 0;
//  if( _loadData )
//  {
//    _readDivisionData.resize( _lineages.size() );
//    _dim += 6;

//    for( std::size_t i=0; i < _lineages.size(); i++ )
//    {
//      std::string name = "/home/necrolyte/Uni/LateralRootGrowth/091115/DivisionOrientations_Late_";
//      name += _lineages.at(i)->getName() + ".csv";
//      this->readDivisionData( name, i );
//    }
//  }

  _data.resize( _dim );
  for( std::size_t d = 0; d < _lineages.size(); d++ )
  {
    NodeFeatureInfo divScheme = _lineages.at(d)->getDivisionType();
    std::map<int,int> cellFiles = _lineages.at(d)->getCellFiles();
    std::map< std::pair<std::size_t, std::size_t>, std::size_t > divTypes;

//    if( _loadData )
//      SVLRDataAnalysisIO::readAngleData( _lineages.at(d), divTypes );

    SArray rotation = _lineages.at(d)->getRotation();
    osg::Matrix rotMat;
    rotMat.postMult( osg::Matrix::rotate( rotation[2]/180.*M_PI, osg::Vec3(0,0,1) ) );
    rotMat.postMult( osg::Matrix::rotate( rotation[1]/180.*M_PI, osg::Vec3(0,1,0) ) );
    rotMat.postMult( osg::Matrix::rotate( rotation[0]/180.*M_PI, osg::Vec3(1,0,0) ) );

    std::vector<std::size_t> numCells;
    std::vector<double> height_LR;
    numCells.resize( maxTimeSteps, 0 );
    height_LR.resize( maxTimeSteps, -5000. );

    SArray ce = _lineages.at(d)->getCenter();
    osg::Vec3 primCen( ce[0], ce[1], ce[2] );
    primCen = primCen * rotMat;

    for( SAdvancedLineageTreeCollection::const_iterator l = _lineages.at(d)->begin();
         l != _lineages.at(d)->end(); ++l )
    {
      for( SLineageTree::const_iterator iter = l->second->begin();
           iter != l->second->end(); ++iter )
      {
        numCells.at( iter->timeStep-1 )++;

        osg::Vec3 pos( iter->getX(), iter->getY(), iter->getZ() );
        pos = pos * rotMat;
        // substract center of root
        pos = pos - primCen;
        if( pos.z() > height_LR.at( iter->timeStep-1 ) )
          height_LR.at( iter->timeStep-1 ) = pos.z();
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
          std::size_t t = iter->timeStep;
          std::size_t id = iter->cellId;
          bool found;
          bool chosenDivType;
          if( _loadData )
          {
            auto divIter = divTypes.find( std::make_pair( id, t ) );
            found = divIter != divTypes.end();
            chosenDivType = divIter->second == divType;

            //std::cout << "id: " << id << " t: " << t << "dT: " << iter2->second << std::endl;
          }
          else
          {
            auto iter2 = divScheme.at(t-1).find( id );
            found = iter2 != divScheme.at(t-1).end();
            chosenDivType = iter2->second == divType;
          }

          // only process specific division type
          if( found && chosenDivType )
          {
            osg::Vec3 divPos( iter->getX(), iter->getY(), iter->getZ() );
            osg::Vec3 prevDivPos = this->getPreviousDivPos( *iter );

            std::size_t cc = iter->getCellCycleId() + 1;

            osg::Vec3 c1( iter->children[0]->getX(),
                iter->children[0]->getY(),
                iter->children[0]->getZ() );

            osg::Vec3 c2( iter->children[1]->getX(),
                iter->children[1]->getY(),
                iter->children[1]->getZ() );

            divPos = divPos * rotMat;
            prevDivPos = prevDivPos * rotMat;
            c1 = c1 * rotMat;
            c2 = c2 * rotMat;

            osg::Vec3 dir1 = c2 - c1;
            osg::Vec3 dir2 = c1 - c2;
            osg::Vec3 growthDir = divPos - prevDivPos;
            dir1.normalize();
            dir2.normalize();
            growthDir.normalize();

            // compute the two angles between previous growth direction
            // and division orientation
            double angle1 = acos( dir1 * growthDir )*180./M_PI;
            //double angle2 = acos( dir2 * growthDir )*180./M_PI;

            double chosenAngle = angle1;
            osg::Vec3 chosenDir = dir1;
            if( angle1 > 90. )
            {
              chosenAngle = 180. - angle1;
              chosenDir = dir2;
            }

            _data.at(0).push_back( chosenDir[0] );
            _data.at(1).push_back( growthDir[0] );
            _data.at(2).push_back( chosenDir[1] );
            _data.at(3).push_back( growthDir[1] );
            _data.at(4).push_back( chosenDir[2] );
            _data.at(5).push_back( growthDir[2] );
            _data.at(6).push_back( chosenAngle );
            _data.at(7).push_back( numCells.at( t-1 ) );
            _data.at(8).push_back( t );
            _data.at(9).push_back( cc );
            _data.at(10).push_back( height_LR.at( t-1 ) );
            _data.at(11).push_back( cellFiles.at( iter->treeId ) );

            if( _loadData )
            {
              auto iter = _readDivisionData.at(d).find( std::make_pair( id, t-1 ) );
              if( iter != _readDivisionData.at(d).end() )
              {
                divData dd = iter->second;
                _data.at(12).push_back( dd.VLRAngle );
                _data.at(13).push_back( dd.distDivToCentre );
                _data.at(14).push_back( dd.primAngle );
                _data.at(15).push_back( dd.distDivToMainRoot );
                _data.at(16).push_back( dd.surfaceAngle );
                _data.at(17).push_back( dd.distDivToSurface );
              }
            }

            _numPoints++;
          }
        }
      }
    }
  }
}

// ---------------------------------------------------------------------

void SCorrelationAnalysis::generateGeneralData( const std::size_t maxTimeSteps )
{
  _dataType = 1;

  _data.resize( _dim );
  for( std::size_t d = 0; d < _lineages.size(); d++ )
  {
    std::map<int,int> cellFiles = _lineages.at(d)->getCellFiles();

    SArray rotation = _lineages.at(d)->getRotation();
    osg::Matrix rotMat;
    rotMat.postMult( osg::Matrix::rotate( rotation[2]/180.*M_PI, osg::Vec3(0,0,1) ) );
    rotMat.postMult( osg::Matrix::rotate( rotation[1]/180.*M_PI, osg::Vec3(0,1,0) ) );
    rotMat.postMult( osg::Matrix::rotate( rotation[0]/180.*M_PI, osg::Vec3(1,0,0) ) );

    std::vector<std::size_t> numCells;
    std::vector<double> height_LR;
    numCells.resize( maxTimeSteps, 0 );
    height_LR.resize( maxTimeSteps, -5000. );

    SArray ce = _lineages.at(d)->getCenter();
    osg::Vec3 primCen( ce[0], ce[1], ce[2] );
    primCen = primCen * rotMat;

    for( SAdvancedLineageTreeCollection::const_iterator l = _lineages.at(d)->begin();
         l != _lineages.at(d)->end(); ++l )
    {
      for( SLineageTree::const_iterator iter = l->second->begin();
           iter != l->second->end(); ++iter )
      {
        numCells.at( iter->timeStep-1 )++;

        osg::Vec3 pos( iter->getX(), iter->getY(), iter->getZ() );
        pos = pos * rotMat;
        // substract center of root
        pos = pos - primCen;
        if( pos.z() > height_LR.at( iter->timeStep-1 ) )
          height_LR.at( iter->timeStep-1 ) = pos.z();
      }
    }

    for( SAdvancedLineageTreeCollection::const_iterator l = _lineages.at(d)->begin();
         l != _lineages.at(d)->end(); ++l )
    {
      for( SLineageTree::const_iterator iter = l->second->begin();
           iter != l->second->end(); ++iter )
      {
        std::size_t t = iter->timeStep;
        std::size_t cc = iter->getCellCycleId() + 1;
        osg::Vec3 pos( iter->getX(), iter->getY(), iter->getZ() );
        pos = pos * rotMat;
        pos = pos - primCen;

        _data.at(0).push_back( pos[0] );
        _data.at(1).push_back( pos[1] );
        _data.at(2).push_back( pos[2] );
        _data.at(3).push_back( numCells.at( t-1 ) );
        _data.at(4).push_back( t );
        _data.at(5).push_back( cc );
        _data.at(6).push_back( height_LR.at( t-1 ) );
        _data.at(7).push_back( cellFiles.at( iter->treeId ) );

        _numPoints++;
      }
    }
  }
}

// ---------------------------------------------------------------------

void SCorrelationAnalysis::renderScatterplotMatrix( osg::ref_ptr<osg::Group> addToThisGroup )
{
  std::vector<std::string> featureNames;
  featureNames.resize( _dim );

  if( _dataType == 0 )
  {
    featureNames.at(0) = "X";
    featureNames.at(1) = "GX";
    featureNames.at(2) = "Y";
    featureNames.at(3) = "GY";
    featureNames.at(4) = "Z";
    featureNames.at(5) = "GZ";
    featureNames.at(6) = "divAngle";
    featureNames.at(7) = "NumCells";
    featureNames.at(8) = "T";
    featureNames.at(9) = "CC";
    featureNames.at(10) = "Height";
    featureNames.at(11) = "CellFile";
    if( _loadData )
    {
      featureNames.at(12) = "VLRAngle";
      featureNames.at(13) = "|Div-Center|";
      featureNames.at(14) = "primAngle";
      featureNames.at(15) = "|Div-Root|";
      featureNames.at(16) = "surfaceAngle";
      featureNames.at(17) = "|Div-Surface|";
    }
  }
  else
  {
    featureNames.at(0) = "X";
    featureNames.at(1) = "Y";
    featureNames.at(2) = "Z";
    featureNames.at(3) = "NumCells";
    featureNames.at(4) = "T";
    featureNames.at(5) = "CC";
    featureNames.at(6) = "Height";
    featureNames.at(7) = "CellFile";
  }

  osg::ref_ptr<SScatterplotMatrix> sm = new SScatterplotMatrix( _data, _numPoints,
                                                                false, featureNames,
                                                                0.15, true );
  sm->setName( "ScatterplotMatrix" );
  sm->setTransformationMatrix( addToThisGroup->asTransform()->asMatrixTransform()->getMatrix() );
  addToThisGroup->addChild( sm );
}

// ---------------------------------------------------------------------

osg::Vec3 SCorrelationAnalysis::getPreviousDivPos( const SLineageTree *node )
{
  if( node->parent )
  {
    const SLineageTree *prev = node->parent;

    while( prev && prev->children.size() != 2 )
    {
      if( prev->parent )
        prev = prev->parent;
      // else return the position of the root node
      else
        return osg::Vec3( prev->getX(), prev->getY(), prev->getZ() );
    }

    // return the position of the previous division node
    return osg::Vec3( prev->getX(), prev->getY(), prev->getZ() );
  }
  // else return the position of the current node
  else
    return osg::Vec3( node->getX(), node->getY(), node->getZ() );
}

// ---------------------------------------------------------------------
