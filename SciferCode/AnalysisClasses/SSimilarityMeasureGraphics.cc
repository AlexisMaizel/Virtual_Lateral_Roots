#include "SSimilarityMeasureGraphics.hh"

#include "graphics/SSelectableSphere.hh"
#include "graphics/SLabelGenerator.hh"
#include "graphics/SRectangle.hh"

#include "SArrow3DAngled.hh"
#include "SArrow3DCylindrical.hh"

#include "SLayerInformationIO.hh"

#include "SInteractiveDivisionArrows.hh"
#include "SInteractiveSurface.hh"
#include "SInteractiveNormals.hh"

#include "kernel/SKernel.hh"

#include <osg/LineWidth>
#include <osg/Material>

#include <osgUtil/SmoothingVisitor>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

namespace SSimilarityMeasureGraphics
{
//----------------------------------------------------------------------------

void drawTimeLineLabels( osg::ref_ptr<osg::Geode> labelGeode,
                         const double stepSize,
                         const double maxDepth,
                         const double coordStart,
                         const double coordMax,
                         double maxXCoord,
                         const bool vertical,
                         const double labelSize,
                         const bool drawTimeLines )
{
  std::vector<std::string> labels;

  double coord = 0.;
  double startLabel = 0.;

  if( maxXCoord == 0. )
    maxXCoord = maxDepth;

  double stepCoord = stepSize * (maxXCoord/maxDepth);

  for( std::size_t i = 0; i <= (std::size_t)(maxDepth/stepSize) + 1; i++ )
  {
    labels.push_back( boost::lexical_cast<std::string>( startLabel ) );

    osg::Vec3 pos;

    if( vertical )
      pos = osg::Vec3( coordStart, -coord, 1. );
    else
      pos = osg::Vec3( coord, coordStart, 1. );

    osg::ref_ptr<osgText::Text> label = SLabelGenerator::getLabel(
          labels[i], pos, osg::Vec4(0.3,0.3,0.3,1.), labelSize );
    if( vertical )
      label->setAlignment( osgText::Text::RIGHT_CENTER );
    else
      label->setAlignment( osgText::Text::CENTER_CENTER );

    label->setCharacterSizeMode( osgText::Text::OBJECT_COORDS );
    label->setAutoRotateToScreen( false );

    labelGeode->addDrawable( label );

    // draw line for time step
    osg::Vec3 rectCoords[2];
    if( vertical )
    {
      rectCoords[0] = osg::Vec3( pos[0], pos[1], -3. );
      rectCoords[1] = osg::Vec3( pos[0] + coordMax, pos[1], -3. );
    }
    else
    {
      rectCoords[0] = osg::Vec3( pos[0], pos[1], -3. );
      rectCoords[1] = osg::Vec3( pos[0], pos[1] + coordMax, -3. );
    }

    if( drawTimeLines )
      SSimilarityMeasureGraphics::drawLine( labelGeode, rectCoords[0],
                                        rectCoords[1], 2.,
                                        osg::Vec4( 0.7, 0.7, 0.7, 1. ) );

    coord += (double)stepCoord;
    startLabel += (double)stepSize;
  }
}

//----------------------------------------------------------------------------

void drawLine( osg::ref_ptr<osg::Geode> lineGeode,
               const osg::Vec3 &pos1, const osg::Vec3 &pos2,
               const double lineWidth, const osg::Vec4 &color )
{
  osg::Vec3 rectCoords[2];
  rectCoords[0] = pos1;
  rectCoords[1] = pos2;
  osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array( 2, rectCoords );
  osg::ref_ptr<osg::Geometry> line = new osg::Geometry;
  line->setVertexArray( vertices );

  // color
  osg::ref_ptr<osg::Vec4Array> colorRect = new osg::Vec4Array;
  colorRect->push_back( color );
  line->setColorArray( colorRect.get() );
  line->setColorBinding( osg::Geometry::BIND_OVERALL );

  // line properties
  line->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::LINES, 0, 2) );
  osg::ref_ptr<osg::LineWidth> linewidth = new osg::LineWidth();
  linewidth->setWidth( lineWidth );
  line->getOrCreateStateSet()->setAttributeAndModes( linewidth, osg::StateAttribute::ON );
  line->getOrCreateStateSet()->setMode( GL_LIGHTING, osg::StateAttribute::ON  );
  line->getOrCreateStateSet()->setMode( GL_LINE_SMOOTH, osg::StateAttribute::ON  );
  line->getOrCreateStateSet()->setMode( GL_BLEND, osg::StateAttribute::ON  );

  lineGeode->addDrawable( line );
}

//----------------------------------------------------------------------------

void drawSequenceRectangle( osg::ref_ptr<osg::Geode> geode,
                            const int divType,
                            const double size,
                            const osg::Vec2 &pos )
{
  osg::Vec4 color;

  DivisionType::DivisionType divT = (DivisionType::DivisionType)divType;

  switch( divT )
  {
  case DivisionType::ANTICLINAL:
    color = osg::Vec4( 1., 0., 0., 1. );
    break;
  case DivisionType::PERICLINAL:
    color = osg::Vec4( 0., 1., 0., 1. );
    break;
  case DivisionType::RADIAL:
    color = osg::Vec4( 0., 0., 1., 1. );
    break;
  }

  osg::ref_ptr<SRectangle> rect = new SRectangle(
        pos[0], pos[1], 1., size, size, color, color, 0. );

  // vertices
  osg::Vec3 rectCoords[5];
  rectCoords[0] = osg::Vec3( pos[0]-size/2., pos[1]-size/2., 1.1 );
  rectCoords[1] = osg::Vec3( pos[0]+size/2., pos[1]-size/2., 1.1 );
  rectCoords[2] = osg::Vec3( pos[0]+size/2., pos[1]+size/2., 1.1 );
  rectCoords[3] = osg::Vec3( pos[0]-size/2., pos[1]+size/2., 1.1 );
  rectCoords[4] = osg::Vec3( pos[0]-size/2., pos[1]-size/2., 1.1 );
  osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array( 5, rectCoords );
  osg::ref_ptr<osg::Geometry> lineRect = new osg::Geometry;
  lineRect->setVertexArray( vertices );

  // color
  osg::ref_ptr<osg::Vec4Array> colorRect = new osg::Vec4Array;
  colorRect->push_back(osg::Vec4( 0., 0., 0., 1. ) );
  lineRect->setColorArray( colorRect.get() );
  lineRect->setColorBinding( osg::Geometry::BIND_OVERALL );

  // line properties
  lineRect->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::LINE_STRIP, 0, 5) );
  osg::ref_ptr<osg::LineWidth> linewidth = new osg::LineWidth();
  linewidth->setWidth( 4. );
  lineRect->getOrCreateStateSet()->setAttributeAndModes( linewidth, osg::StateAttribute::ON );

  geode->addDrawable( rect );
  geode->addDrawable( lineRect );
}


//----------------------------------------------------------------------------

void renderDivisionSequence( const SLineageTree *node,
                               const NodeFeatureInfo &cellDivisionScheme,
                               const std::size_t startTime,
                               osg::ref_ptr<osg::Group> group )
{
  SLineageTree *newNode = new SLineageTree( *node );

  std::vector<int> divisionSequence;

  osg::ref_ptr<osg::Geode> geode = new osg::Geode;

  double size = 4.;

  double initXPos = newNode->getX();
  double initYPos = 3.;//-newNode->getY() - 5.;

  while( newNode )
  {
    if( newNode->children.size() == 2 )
    {
      // divScheme is either:
      // 0 -> periclinal
      // 1 -> anticlinal
      // 2 -> radial
      int divType = cellDivisionScheme.at(newNode->timeStep-1).at( newNode->cellId );
      divisionSequence.push_back( divType );
    }

    newNode = newNode->parent;

    if( newNode )
    {
      if( newNode->timeStep == startTime - 1 )
        break;
    }
  }

  int j = 1;
  for( int i=divisionSequence.size()-1; i >= 0; i--, j++ )
  {
    SSimilarityMeasureGraphics::drawSequenceRectangle( geode, divisionSequence[i], size,
                                                   osg::Vec2( initXPos, initYPos - j*size ) );
  }

  group->addChild( geode );
}

// ---------------------------------------------------------------------

void drawComparisonRectangle( osg::ref_ptr<osg::Geode> geode,
                              const double sizeX,
                              const double sizeY,
                              const osg::Vec3 &pos,
                              double ratio,
                              const std::size_t appearances )
{
  if( ratio > 1. )
    ratio = 1.;

  if( ratio < 0. )
    ratio = 0.;

  double rectXSize1 = sizeX*ratio;
  double rectXSize2 = sizeX*(1.-ratio);
  double newRectPos1 = pos[0]-sizeX/2. + rectXSize1/2.;
  double newRectPos2 = newRectPos1 + rectXSize1/2. + rectXSize2/2.;

  osg::Vec4 firstColor1( 255./255., 255./255., 178./255., 1. );
  osg::Vec4 firstColor2( 189./255., 0./255., 38./255., 1. );
  osg::Vec4 color1 = firstColor1 * ratio + firstColor2 * (1.-ratio);
  osg::Vec4 color2 = firstColor1 * (1.-ratio) + firstColor2 * ratio;
  osg::ref_ptr<SRectangle> rect1 = new SRectangle(
        newRectPos1, pos[1], 1., rectXSize1, sizeY, color1, color1, 0. );

  osg::ref_ptr<SRectangle> rect2 = new SRectangle(
        newRectPos2, pos[1], 1., rectXSize2, sizeY, color2, color2, 0. );

  osg::Vec3 textPos;
  std::string label;
  if( ratio > 0.5 )
  {
    textPos = osg::Vec3( newRectPos1, pos[1], 2. );
    label = std::to_string( ratio*100. );
  }
  else
  {
    textPos = osg::Vec3( newRectPos2, pos[1], 2. );
    label = std::to_string( (1.-ratio)*100. );
  }

  label = label.substr( 0, 4 );
  osg::ref_ptr<osgText::Text> labelText = SLabelGenerator::getLabel(
        label, textPos, osg::Vec4( 0., 0., 0., 1.), 4.,
        osgText::Text::CENTER_CENTER );
  labelText->setCharacterSizeMode( osgText::Text::OBJECT_COORDS );
  labelText->setAutoRotateToScreen( false );

  std::string appLabel = std::to_string( appearances );
  osg::Vec3 appTextPos( pos[0] + sizeX/2. + sizeX/3. + 2.5, pos[1] - sizeY/2., 2. );
  osg::ref_ptr<osgText::Text> appText = SLabelGenerator::getLabel(
        appLabel, appTextPos, osg::Vec4( 0., 0., 0., 1.), 3.,
        osgText::Text::LEFT_BOTTOM );
  appText->setCharacterSizeMode( osgText::Text::OBJECT_COORDS );
  appText->setAutoRotateToScreen( false );

  geode->addDrawable( labelText );
  geode->addDrawable( appText );
  geode->addDrawable( rect1 );
  geode->addDrawable( rect2 );

  osg::Vec3 rectCoords[2];
  rectCoords[0] = osg::Vec3( pos[0]-sizeX/2. + rectXSize1, pos[1]-sizeY/2.-0.2, 2. );
  rectCoords[1] = osg::Vec3( pos[0]-sizeX/2. + rectXSize1, pos[1]+sizeY/2.+0.2, 2. );

  SSimilarityMeasureGraphics::drawLine( geode, rectCoords[0],
                                    rectCoords[1], 3.,
                                    osg::Vec4( 0., 0., 0., 1. ) );
}

// ---------------------------------------------------------------------

void drawSequenceRectangle( osg::ref_ptr<osg::Geode> geode,
                            const std::string &label,
                            const double labelSize,
                            const double sizeX,
                            const double sizeY,
                            const osg::Vec3 &pos,
                            const osg::Vec4 &color )
{
  osg::ref_ptr<SRectangle> rect = new SRectangle(
        pos[0], pos[1], 1., sizeX, sizeY, color, color, 0. );

  osg::ref_ptr<osgText::Text> labelText = SLabelGenerator::getLabel(
        label, osg::Vec3( pos[0], pos[1], 5. ), osg::Vec4(0.,0.,0.,1.), labelSize,
        osgText::Text::CENTER_CENTER );
  labelText->setCharacterSizeMode( osgText::Text::OBJECT_COORDS );
  labelText->setAutoRotateToScreen( false );

  geode->addDrawable( labelText );
  geode->addDrawable( rect );
}

// ---------------------------------------------------------------------

void drawDoubleSequenceRectangle( osg::ref_ptr<osg::Geode> geode,
                                  const std::string &label1,
                                  const std::string &label2,
                                  const double labelSize,
                                  const double sizeX,
                                  const double sizeY,
                                  const osg::Vec3 &pos,
                                  const osg::Vec4 &color1,
                                  const osg::Vec4 &color2 )
{
  osg::ref_ptr<SRectangle> rect1 = new SRectangle(
        pos[0], pos[1]+sizeY/2., 10.5, sizeX, sizeY, color1, color1, 0. );

  osg::ref_ptr<SRectangle> rect2 = new SRectangle(
        pos[0], pos[1]-sizeY/2., 10.5, sizeX, sizeY, color2, color2, 0. );

  osg::ref_ptr<osgText::Text> labelText1 = SLabelGenerator::getLabel(
        label1, osg::Vec3( pos[0], pos[1]+sizeY/2., 11. ), osg::Vec4(0.,0.,0.,1.), labelSize,
        osgText::Text::CENTER_CENTER );

  osg::ref_ptr<osgText::Text> labelText2 = SLabelGenerator::getLabel(
        label2, osg::Vec3( pos[0], pos[1]-sizeY/2., 11. ), osg::Vec4(0.,0.,0.,1.), labelSize,
        osgText::Text::CENTER_CENTER );
  labelText1->setCharacterSizeMode( osgText::Text::OBJECT_COORDS );
  labelText1->setAutoRotateToScreen( false );
  labelText2->setCharacterSizeMode( osgText::Text::OBJECT_COORDS );
  labelText2->setAutoRotateToScreen( false );

  geode->addDrawable( labelText1 );
  geode->addDrawable( labelText2 );
  geode->addDrawable( rect1 );
  geode->addDrawable( rect2 );
}

// ---------------------------------------------------------------------

void drawComparisonBox( osg::ref_ptr<osg::Geode> geode,
                        const osg::Vec3 &leftTop,
                        const osg::Vec3 &rightBottom,
                        const std::string &compLabel,
                        const std::string &appLabel,
                        const double labelSize )
{
  double perc = std::atof( compLabel.c_str() );

  osg::Vec4 color;
  if( perc >= 50. )
    color = osg::Vec4( 0., 1., 0., 1. );
  else
    color = osg::Vec4( 1., 0., 0., 1. );

  osg::ref_ptr<SRectangle> rect = new SRectangle( leftTop, rightBottom,
        osg::Vec4( 1., 1., 1., 1. ), 0.5, color );

  osg::Vec3 midPos( (rightBottom[0] + leftTop[0])/2., (leftTop[1] + rightBottom[1])/2., 12. );

  osg::ref_ptr<osgText::Text> labelText = SLabelGenerator::getLabel(
        compLabel.substr( 0, 4 ) + "% / " + appLabel,
        midPos, osg::Vec4(0.,0.,0.,1.), labelSize, osgText::Text::CENTER_CENTER );

//  osg::ref_ptr<osgText::Text> labelText1 = SLabelGenerator::getLabel(
//        compLabel.substr( 0, 4 )+"%", osg::Vec3( leftTop[0] + 0.5, rightBottom[1] + 1., 12. ),
//        osg::Vec4(0.,0.,0.,1.), labelSize, osgText::Text::LEFT_CENTER );

//  osg::ref_ptr<osgText::Text> labelText2 = SLabelGenerator::getLabel(
//        appLabel, osg::Vec3( rightBottom[0] - 0.5, leftTop[1] - 1., 12. ),
//        osg::Vec4(0.,0.,0.,1.), labelSize, osgText::Text::RIGHT_CENTER );

  labelText->setCharacterSizeMode( osgText::Text::OBJECT_COORDS );
  labelText->setAutoRotateToScreen( false );

  geode->addDrawable( labelText );
  geode->addDrawable( rect );
}

// ---------------------------------------------------------------------

void addCylinderBetweenPoints( const osg::Vec3& p1, const osg::Vec3& p2,
                               double rad, const osg::Vec4& color,
                               osg::ref_ptr<osg::Group> addToThisGroup,
                               const std::string &drawableName )
{
  double h = (p1-p2).length();
  osg::Vec3 c( (p1.x()+p2.x()) / 2, (p1.y()+p2.y()) / 2, (p1.z()+p2.z()) / 2 );

  // This is the default direction for the cylinders to face in OpenGL
  osg::Vec3 z = osg::Vec3( 0, 0, 1 );

  // The distance (as vector) between the two points
  osg::Vec3 p = p1 - p2;

  // The cross product (the axis of rotation)
  osg::Vec3 t = z ^ p;

  double angle = acos( (z * p) / p.length() );

  osg::ref_ptr<osg::Cylinder> cylinder = new osg::Cylinder( c, rad, h );
  cylinder->setRotation( osg::Quat( angle, osg::Vec3(t.x(), t.y(), t.z()) ) );
  osg::ref_ptr<osg::ShapeDrawable> cylinderDrawable =
      new osg::ShapeDrawable( cylinder );
  cylinderDrawable->setColor( color );
  osg::TessellationHints *hints = new osg::TessellationHints;
  hints->setTessellationMode( osg::TessellationHints::USE_TARGET_NUM_FACES );
  hints->setTargetNumFaces( 15 );
  cylinderDrawable->setTessellationHints( hints );
  cylinderDrawable->setName( drawableName );

  osg::ref_ptr<osg::Geode> geode = new osg::Geode;
  geode->addDrawable( cylinderDrawable );
  addToThisGroup->addChild( geode );
}

//----------------------------------------------------------------------------//

void addArrowBetweenPoints( const osg::Vec3& p1, const osg::Vec3& p2,
                            const osg::Vec3& normal, const osg::Vec4& color,
                            const osg::Vec4& tipColor,
                            osg::ref_ptr<osg::Group> addToThisGroup )
{
  double arrowParam1 = 2.5;
  double arrowParam2 = 0.3;

  osg::Geode *arrow = new SArrow3DCylindrical( arrowParam1, arrowParam2, color, tipColor );

  osg::MatrixTransform *mt = new osg::MatrixTransform(
        SArrow3DAngled::getTransformationMatrix( p1, p2, normal ) );
  mt->addChild( arrow );
  addToThisGroup->addChild( mt );
}

//----------------------------------------------------------------------------//

void addSimpleArrowBetweenPoints( const osg::Vec3& p1, const osg::Vec3& p2,
                                  const osg::Vec4& color, const osg::Vec4& tipColor,
                                  osg::Group *addToThisGroup, const double radius )
{
  double arrowParam2 = 0.15;

  osg::Vec3 dir = p2 - p1;
  dir.normalize();
  osg::Vec3 normal = -dir ^ dir;
  normal.normalize();

  osg::Geode *arrow = new SArrow3DCylindrical( radius, arrowParam2, color, tipColor );

  osg::MatrixTransform *mt = new osg::MatrixTransform(
        SArrow3DAngled::getTransformationMatrix( p1, p2, normal ) );
  mt->addChild( arrow );
  addToThisGroup->addChild( mt );
}

// ---------------------------------------------------------------------

void renderCylindricalDivisionAnalysis( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                                        osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                                        osg::ref_ptr<SSliderCallback> sliderCallback,
                                        const std::size_t numTimeSteps,
                                        const std::size_t numCellCycles,
                                        const bool onlyMasterCellFile,
                                        const std::vector<osg::Vec3> &centresOfMass,
                                        const std::size_t sliderType )
{
  std::size_t totalSteps;
  osg::ref_ptr<SInteractiveNormals> in = new SInteractiveNormals();
  osg::ref_ptr<SInteractiveDivisionArrows> ida = new SInteractiveDivisionArrows();
  if( sliderType == 0 )
  {
    in->addUpdateCallback( sliderCallback );
    in->setName( "Divisions" );
    addToThisGroup->addChild( in );
    totalSteps = numTimeSteps;
  }
  else
  {
    ida->addUpdateCallback( sliderCallback );
    ida->setName( "Divisions" );
    addToThisGroup->addChild( ida );
    totalSteps = numCellCycles;
  }

  osg::Matrix rotMat = addToThisGroup->getMatrix();

  // VLR center
  SArray ce = lineages->getCenter();
  osg::Vec3 VLRCenter( ce[0], ce[1], ce[2] );
  VLRCenter = VLRCenter * rotMat;

  osg::Vec3 domeTip;
  // center is the dome tip of the last time step
  if( lineages->getName() == "120830" )
    domeTip = osg::Vec3( 277.2, 215.1, 72.3 );
  else if( lineages->getName() == "121204" )
    domeTip = osg::Vec3( 347.0, 300.1, 41.3 );
  else if( lineages->getName() == "121211" )
    domeTip = osg::Vec3( 278.2, 251.8, 145.6 );
  else if( lineages->getName() == "130508" )
    domeTip = osg::Vec3( 248.6, 144.5, 337.9 );
  else if( lineages->getName() == "130607" )
    domeTip = osg::Vec3( 160.1, 388.0, 298.0 );

  const double rhoEps = 6.;
  const double thetaEps = 5.;

  std::map<int,int> cellFiles = lineages->getCellFiles();
  std::vector<osg::ref_ptr<osg::Group> > groupVec;
  for( std::size_t c = 0; c < totalSteps; c++ )
    groupVec.push_back( new osg::Group );

  std::map< std::pair<std::size_t, std::size_t>, std::string> info;

  for( SAdvancedLineageTreeCollection::const_iterator l = lineages->begin();
       l != lineages->end(); ++l )
  {
    SLineageTree *tree = (*lineages)[l->first];

    if( onlyMasterCellFile && cellFiles.find( tree->cellId )->second != 0 )
      continue;

    for( SLineageTree::const_iterator nodeIt = tree->begin();
         nodeIt != tree->end(); ++nodeIt )
    {
      if( nodeIt->children.size() != 2 )
        continue;

      double m[2],rho[2],theta[2];
      osg::Vec3 divPos( nodeIt->getX(), nodeIt->getY(), nodeIt->getZ() );
      std::size_t step;
      if( sliderType == 0 )
        step = nodeIt->timeStep-1;
      else
        step = nodeIt->getCellCycleId();

      // consider centre of mass for the daughter time step
      osg::Vec3 cen( centresOfMass.at( nodeIt->timeStep ) );
      cen = cen * rotMat;

      std::vector<osg::Vec3> daughters;

      // first child
      SLineageTree::const_iterator iter1 = nodeIt->children[0];
      while( iter1->children.size() != 2 )
      {
        if( iter1->children.size() == 0 )
          break;
        else
          ++iter1;
      }

      osg::Vec3 c1 = osg::Vec3( iter1->getX(), iter1->getY(), iter1->getZ() );
      daughters.push_back( c1 );

      // second child
      SLineageTree::const_iterator iter2 = nodeIt->children[1];
      while( iter2->children.size() != 2 )
      {
        if( iter2->children.size() == 0 )
          break;
        else
          ++iter2;
      }

      osg::Vec3 c2 = osg::Vec3( iter2->getX(), iter2->getY(), iter2->getZ() );
      daughters.push_back( c2 );

      for( std::size_t c = 0; c < 2; c++ )
      {
        //        osg::Vec3 p( nodeIt->children[c]->getX(),
        //                     nodeIt->children[c]->getY(),
        //                     nodeIt->children[c]->getZ() );
        osg::Vec3 p = daughters.at(c);
        p = p * rotMat;
        //p = p - cen;
        //p = p - domeTip;
        p = p - VLRCenter;
        m[c] = p.y();
        rho[c] = std::sqrt( p.x()*p.x() + p.z()*p.z() );
        theta[c] = std::atan2( p.z(), p.x() ) * 180./M_PI;
      }

      double rhoDiff = std::fabs( rho[0] - rho[1] );
      double thetaDiff = std::fabs( theta[0] - theta[1] );

      if( thetaDiff > 180. )
        thetaDiff = 360. - thetaDiff;

      osg::Vec4 color;
      if( rhoDiff < rhoEps && thetaDiff < thetaEps )
        color = osg::Vec4( 1., 0., 0., 1. );
      else
        color = osg::Vec4( 0.8, 0.8, 0.8, 1. );

      std::string inf = "t: " + std::to_string(nodeIt->timeStep) +
          " rho: " + std::to_string(rhoDiff) + " theta: " + std::to_string(thetaDiff);
      info.insert( std::make_pair( std::make_pair( nodeIt->timeStep, nodeIt->cellId ), inf ) );

      //      osg::Vec3 c1( nodeIt->children[0]->getX(),
      //          nodeIt->children[0]->getY(),
      //          nodeIt->children[0]->getZ() );

      //      osg::Vec3 c2( nodeIt->children[1]->getX(),
      //          nodeIt->children[1]->getY(),
      //          nodeIt->children[1]->getZ() );

      SSimilarityMeasureGraphics::addCylinderBetweenPoints( daughters.at(0), daughters.at(1), 2., color,
                                                            groupVec.at(step), inf );

//      SSimilarityMeasureGraphics::addCylinderBetweenPoints( c1, c2, 2., color,
//                                                        groupVec.at(step) );
    }
  }

  for( std::size_t c = 0; c < totalSteps; c++ )
  {
    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    if( sliderType == 0 )
    {
      osg::Vec3 pos = centresOfMass.at(c);
      osg::ref_ptr<SSelectableSphere> sphere = new SSelectableSphere( pos, osg::Vec4( 0., 0., 1., 1. ), 5. );
      geode->addDrawable( sphere );
    }

    osg::ref_ptr<SSelectableSphere> domeSphere = new SSelectableSphere( domeTip * addToThisGroup->getInverseMatrix(), osg::Vec4( 0., 1., 0., 1. ), 5. );
    geode->addDrawable( domeSphere );
    osg::ref_ptr<SSelectableSphere> VLRSphere = new SSelectableSphere( VLRCenter * addToThisGroup->getInverseMatrix(), osg::Vec4( 0., 1., 1., 1. ), 5. );
    geode->addDrawable( VLRSphere );
    groupVec.at(c)->addChild( geode );
  }

  if( sliderType == 0 )
  {
    for( std::size_t c = 0; c < groupVec.size(); c++ )
      in->addNormals( c+1, groupVec.at(c) );
  }
  else
  {
    for( std::size_t c = 0; c < groupVec.size(); c++ )
      ida->addDivisionArrows( c+1, groupVec.at(c) );
  }
}

// ---------------------------------------------------------------------

void renderDivisions( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                      osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                      osg::ref_ptr<SSliderCallback> sliderCallback,
                      osg::ref_ptr<osg::Vec4Array> layerColors,
                      const NodeFeatureInfo &layerValues,
                      const bool renderTracks,
                      const bool onlyMasterCellFile,
                      const std::size_t divisionColorType,
                      const bool renderArrows )
{
  int numTrees = 0;

  // initialize callback for division arrows over several cell cycles
  osg::ref_ptr<SInteractiveDivisionArrows> ida = new SInteractiveDivisionArrows();
  ida->addUpdateCallback( sliderCallback );
  if( renderArrows )
    ida->setName( "Division Arrows" );
  else
    ida->setName( "Division Cylinders" );
  addToThisGroup->addChild( ida );

  std::map<int,int> cellFiles = lineages->getCellFiles();
  NodeFeatureInfo divisionTypes = lineages->getDivisionType();

  // for each cell cycle store all division arrows of
  // all lineage trees
  std::vector<osg::ref_ptr<osg::Group> > cycleGroupVector;

  // compute the timestep of the first division over all lineages
  std::vector<unsigned int> firstDivTS;

  for( SAdvancedLineageTreeCollection::const_iterator l = lineages->begin();
       l != lineages->end(); ++l )
  {
    SLineageTree *tree = (*lineages)[l->first];

    // only process these trees which belong to the master cell file
    if( onlyMasterCellFile )
    {
      if( cellFiles.find( tree->cellId )->second != 0 )
        continue;
    }

    for( SLineageTree::const_iterator nodeIt = tree->begin();
         nodeIt != tree->end(); ++nodeIt )
    {
      if( nodeIt->children.size() != 2 )
        continue;

      unsigned int nd = nodeIt->getCellCycleId();

      if( firstDivTS.size() <= nd )
        firstDivTS.resize( nd+1, boost::numeric::bounds<unsigned int>::highest() );

      if( firstDivTS[nd] > nodeIt->timeStep )
        firstDivTS[nd] = nodeIt->timeStep;

      nd++;
    }

    ++numTrees;
  }

  std::vector<osg::Vec4> divisionColors;
  // anticlinal
  divisionColors.push_back( osg::Vec4( 1., 0., 0., 1. ) );
  // periclinal
  divisionColors.push_back( osg::Vec4( 0., 1., 0., 1. ) );
  // radial
  divisionColors.push_back( osg::Vec4( 0., 0., 1., 1. ) );

  // init the division glyphs
  numTrees = 0;
  for( SAdvancedLineageTreeCollection::const_iterator l = lineages->begin();
       l != lineages->end(); ++l )
  {
    SLineageTree *tree = (*lineages)[l->first];

    // only process these trees which belong to the master cell file
    if( onlyMasterCellFile )
    {
      if( cellFiles.find( tree->cellId )->second != 0 )
        continue;
    }

    osg::ref_ptr<osg::Group> linGroup = new osg::Group;
    linGroup->setName( "Lineage " + boost::lexical_cast<std::string>( l->first ) );
    osg::ref_ptr<osg::Group> tracks = new osg::Group;
    tracks->setName( "Cell tracks" );
    linGroup->addChild( tracks );

    if( renderTracks )
      addToThisGroup->addChild( linGroup );

    std::vector< osg::ref_ptr<osg::Group> > arrowGroup;

    for( SLineageTree::const_iterator nodeIt = tree->begin();
         nodeIt != tree->end(); ++nodeIt )
    {
      // handle divisions
      if( nodeIt->children.size() == 2 )
      {
        // get positions of the children
        osg::Vec3 n;
        osg::Vec3 p0(nodeIt->children[0]->getX(),
            nodeIt->children[0]->getY(),
            nodeIt->children[0]->getZ());
        osg::Vec3 p1(nodeIt->getX(),
                     nodeIt->getY(),
                     nodeIt->getZ());
        osg::Vec3 p2(nodeIt->children[1]->getX(),
            nodeIt->children[1]->getY(),
            nodeIt->children[1]->getZ());
        osg::Vec3 d1 = p0-p1;
        d1.normalize();
        osg::Vec3 d2 = p2-p1;
        d2.normalize();
        n = d1 ^ d2;
        n.normalize();

        // compute the number of divisions before this one and the elapsed
        // time since the last division
        std::size_t nDivs = 0, dt = 0;
        SLineageTree const* nIt = *nodeIt;
        while( nIt != 0 )
        {
          if( nIt->children.size() > 1 )
          {
            nDivs++;
            if( dt == 0 )
              dt = nodeIt->timeStep - nIt->timeStep;
          }

          nIt = nIt->parent;
        }

        // render cylinder connecting consecutive divisions
        nIt = *nodeIt;
        while( nIt->parent != 0 )
        {
          nIt = nIt->parent;

          if( (nIt->parent && nIt->parent->children.size() > 1) || nIt->parent == 0 )
          {
            osg::Vec3 p(nIt->getX(), nIt->getY(), nIt->getZ());
            if( renderTracks )
            {
              SSimilarityMeasureGraphics::addCylinderBetweenPoints( p, p1, 0.2,
                                                                osg::Vec4( 220/255., 220/255., 220/255., 1. ),
                                                                tracks );
            }

            break;
          }
        }

        while( nDivs > arrowGroup.size() )
        {
          osg::ref_ptr<osg::Group> group = new osg::Group;
          group->setName( "Cell division " + boost::lexical_cast<std::string>(arrowGroup.size()+1) );
          arrowGroup.push_back( group );
        }

        int layerValue = layerValues.at(nodeIt->timeStep-1).at(nodeIt->cellId) - 1;
        osg::Vec4 color, colorFirstChild, colorSecondChild;
        if( divisionColorType == 0 )
          color = divisionColors.at(divisionTypes.at(nodeIt->timeStep-1).find(nodeIt->cellId)->second);
        else
          color = (*layerColors)[layerValue%layerColors->size()];

        layerValue = layerValues.at(nodeIt->children[0]->timeStep-1).at(nodeIt->children[0]->cellId) - 1;
        if( divisionColorType == 0 )
          colorFirstChild = color;
        else
          colorFirstChild = (*layerColors)[layerValue%layerColors->size()];

        layerValue = layerValues.at(nodeIt->children[1]->timeStep-1).at(nodeIt->children[1]->cellId) - 1;
        if( divisionColorType == 0 )
          colorSecondChild = color;
        else
          colorSecondChild = (*layerColors)[layerValue%layerColors->size()];

        if( renderArrows )
        {
          SSimilarityMeasureGraphics::addArrowBetweenPoints( p0, p1, n, color, colorFirstChild, arrowGroup[nDivs-1] );

          if( arrowGroup[nDivs-1]->getNumChildren() > 0 )
            arrowGroup[nDivs-1]->getChild( arrowGroup[nDivs-1]->getNumChildren()-1)
                ->setName( "Cell " + boost::lexical_cast<std::string>(nodeIt->cellId) );

          SSimilarityMeasureGraphics::addArrowBetweenPoints( p2, p1, n, color, colorSecondChild, arrowGroup[nDivs-1] );
        }
        else
          SSimilarityMeasureGraphics::addCylinderBetweenPoints( p0, p2, 2., color, arrowGroup[nDivs-1] );

        if( renderTracks )
        {
          SSimilarityMeasureGraphics::addCylinderBetweenPoints( p0, p1, 0.2,
                                                            osg::Vec4( 220/255., 220/255., 220/255., 1.),
                                                            tracks );
          SSimilarityMeasureGraphics::addCylinderBetweenPoints( p2, p1, 0.2,
                                                            osg::Vec4( 220/255., 220/255., 220/255., 1.),
                                                            tracks );
        }
      }
    }

    ++numTrees;

    // after the whole tree is traversed, add the rendered division
    // arrows to the corresponding cycle groups

    // first resize the cycle vector if required
    while( cycleGroupVector.size() < arrowGroup.size() )
    {
      osg::ref_ptr<osg::Group> newGroup = new osg::Group;
      cycleGroupVector.push_back( newGroup );
    }

    // add the stored arrows
    for( std::size_t c = 0; c < arrowGroup.size(); c++ )
      cycleGroupVector.at(c)->addChild( arrowGroup.at(c) );
  }

  // after traversing all trees in the data set
  // add the division arrows to the cycle callback
  for( std::size_t c = 0; c < cycleGroupVector.size(); c++ )
    ida->addDivisionArrows( c+1, cycleGroupVector.at(c) );
}


// ---------------------------------------------------------------------

void render3DCross( const osg::Vec3 &center,
                     osg::ref_ptr<osg::MatrixTransform> addToThisGroup )
{
  osg::ref_ptr<osg::MatrixTransform> lineGroup = new osg::MatrixTransform;
  lineGroup->setName( "Root center" );

  double lineLength = 150.;
  double rad = 1.0;

  // line in x direction
  osg::ref_ptr<osg::Group> groupLine1 = new osg::Group;
  groupLine1->setName( "Line in x direction" );
  SSimilarityMeasureGraphics::addCylinderBetweenPoints( osg::Vec3( center[0] - lineLength/2., center[1], center[2] ),
                                            osg::Vec3( center[0] + lineLength/2., center[1], center[2] ),
                                            rad, osg::Vec4( 0., 0., 0., 1.),
                                            groupLine1 );

  lineGroup->addChild( groupLine1 );

  // line in y direction
  osg::ref_ptr<osg::Group> groupLine2 = new osg::Group;
  groupLine2->setName( "Line in y direction" );
  SSimilarityMeasureGraphics::addCylinderBetweenPoints( osg::Vec3( center[0], center[1] - lineLength/2., center[2] ),
                                            osg::Vec3( center[0], center[1] + lineLength/2., center[2] ),
                                            rad, osg::Vec4( 0., 0., 0., 1.),
                                            groupLine2 );

  lineGroup->addChild( groupLine2 );

  // line in z direction
  osg::ref_ptr<osg::Group> groupLine3 = new osg::Group;
  groupLine3->setName( "Line in z direction" );
  SSimilarityMeasureGraphics::addCylinderBetweenPoints( osg::Vec3( center[0], center[1], center[2] - lineLength/2. ),
                                            osg::Vec3( center[0], center[1], center[2] + lineLength/2. ),
                                            rad, osg::Vec4( 0., 0., 0., 1.),
                                            groupLine3 );

  lineGroup->addChild( groupLine3 );

  // avoid the rotation of the center lines
  // by setting the matrix to the inverse one
  lineGroup->setMatrix( addToThisGroup->getInverseMatrix() );

  addToThisGroup->addChild( lineGroup );
}


// ---------------------------------------------------------------------

void drawTriangle( const osg::ref_ptr<osg::Vec3Array> points,
                   const osg::Vec4 &color,
                   osg::ref_ptr<osg::Group> addToThisGroup )
{
  osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
  geom->setVertexArray( points );
  geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::TRIANGLES, 0, points->size() ));

  // assign color via material
  osg::ref_ptr<osg::StateSet> state = geom->getOrCreateStateSet();
  osg::ref_ptr<osg::Material> material = new osg::Material;
  material->setDiffuse( osg::Material::FRONT_AND_BACK, color );
  state->setAttribute( material.get() );

  osg::ref_ptr<osg::Geode> geode = new osg::Geode;
  geode->addDrawable( geom );

  // generate normal of triangle
  osgUtil::SmoothingVisitor sv;
  geode->accept( sv );

  addToThisGroup->addChild( geode );
}

// ---------------------------------------------------------------------

void drawPlane( const osg::ref_ptr<osg::Vec3Array> points,
                const osg::Vec4 &color,
                osg::ref_ptr<osg::Group> addToThisGroup )
{
  osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
  geom->setVertexArray( points );
  geom->addPrimitiveSet( new osg::DrawArrays( osg::PrimitiveSet::POLYGON, 0, points->size() ));
  osg::ref_ptr<osg::Vec4Array> colorA = new osg::Vec4Array;
  colorA->push_back( color );
  geom->setColorArray( colorA );

  // assign color via material
  osg::ref_ptr<osg::StateSet> state = geom->getOrCreateStateSet();
  osg::ref_ptr<osg::Material> material = new osg::Material;
  material->setDiffuse( osg::Material::FRONT_AND_BACK, color );
  state->setAttribute( material.get() );

  osg::ref_ptr<osg::Geode> geode = new osg::Geode;
  geode->addDrawable( geom );

  // generate normal of triangle
  osgUtil::SmoothingVisitor sv;
  geode->accept( sv );

  addToThisGroup->addChild( geode );
}


// ---------------------------------------------------------------------

void computeAndRenderSurfacesPerTimeStep( const bool render,
                                          const bool output,
                                          const int type,
                                          osg::ref_ptr<SSliderCallback> sliderCallback,
                                          boost::shared_ptr<cellTimeVector> cTS,
                                          osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                                          const bool wireframe,
                                          const double surfaceOffset )
{
  SInteractiveSurface *is = new SInteractiveSurface();
  osg::ref_ptr<osg::Group> volumeGroup = new osg::Group;
  if( type == 1 )
    is->setName( "Alpha Shape" );
  else
    is->setName( "Convex Hull" );
  is->addUpdateCallback( sliderCallback );
  is->addChild( volumeGroup );

  std::vector< std::pair<double,std::size_t> > volumes;

  for( std::size_t t = 0; t < cTS->size(); t++ )
  {
    SKernel::get()->setProgress( (double)t/(double)(cTS->size() - 1),
                                 "Surface generation" );

    osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array;

    BOOST_FOREACH( const SLineageTree *tree, cTS->at(t) )
    {
      points->push_back( osg::Vec3( tree->getX(), tree->getY(), tree->getZ() ) );
    }

    osg::Vec4 color( 0., 1., 0., 1. );

    if( wireframe )
      color = osg::Vec4( 0., 0., 0., 1. );

    if( type == 1 )
    {
      boost::shared_ptr<SAlphaShape> tri = boost::shared_ptr<SAlphaShape>(
            new SAlphaShape( points, t, color, render, wireframe, surfaceOffset ) );

      if( output )
        volumes.push_back( std::make_pair( tri->computeVolume(), points->size() ) );

      if( render )
        is->addTriangulation( t+1, 0, tri->getTriangles() );
    }
    else
    {
      boost::shared_ptr<SConvexHull> tri = boost::shared_ptr<SConvexHull>(
            new SConvexHull( points, t, color, render, wireframe ) );

      if( output )
        volumes.push_back( std::make_pair( tri->computeVolume(), points->size() ) );

      if( render )
        is->addTriangulation( t+1, 0, tri->getTriangles() );
    }
  }

  if( output )
    SLayerInformationIO::writeVolume( "/tmp/volume.csv", volumes );

  if( render )
    addToThisGroup->addChild( is );
}

// ---------------------------------------------------------------------

void computeAndRenderSurfacesPerTimeStep( const bool render,
                                          const bool output,
                                          const int type,
                                          osg::ref_ptr<SSliderCallback> sliderCallback,
                                          boost::shared_ptr<cellTimeVector> cTS,
                                          const std::vector<int> &averagedCellCycleTimeSteps,
                                          osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                                          const bool wireframe,
                                          const double surfaceOffset )
{
  SInteractiveDivisionArrows *is = new SInteractiveDivisionArrows();
  osg::ref_ptr<osg::Group> volumeGroup = new osg::Group;
  if( type == 1 )
    is->setName( "Alpha Shape" );
  else
    is->setName( "Convex Hull" );
  is->addUpdateCallback( sliderCallback );
  is->addChild( volumeGroup );

  std::vector< std::pair<double,std::size_t> > volumes;

  for( std::size_t i = 0; i < averagedCellCycleTimeSteps.size(); i++ )
  {
    std::size_t t = averagedCellCycleTimeSteps.at(i);

    osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array;

    BOOST_FOREACH( const SLineageTree *tree, cTS->at(t) )
    {
      points->push_back( osg::Vec3( tree->getX(), tree->getY(), tree->getZ() ) );
    }

    osg::Vec4 color( 0., 1., 0., 1. );

    if( wireframe )
      color = osg::Vec4( 0., 0., 0., 1. );

    if( type == 1 )
    {
      boost::shared_ptr<SAlphaShape> tri = boost::shared_ptr<SAlphaShape>(
            new SAlphaShape( points, t, color, render, wireframe, surfaceOffset ) );

      if( output )
        volumes.push_back( std::make_pair( tri->computeVolume(), points->size() ) );

      if( render )
        is->addDivisionArrows( i+1, tri->getTriangles() );
    }
    else
    {
      boost::shared_ptr<SConvexHull> tri = boost::shared_ptr<SConvexHull>(
            new SConvexHull( points, t, color, render, wireframe ) );

      if( output )
        volumes.push_back( std::make_pair( tri->computeVolume(), points->size() ) );

      if( render )
        is->addDivisionArrows( i+1, tri->getTriangles() );
    }
  }

  if( output )
    SLayerInformationIO::writeVolume( "/tmp/volume.csv", volumes );

  if( render )
    addToThisGroup->addChild( is );
}

// ---------------------------------------------------------------------

void computeAndRenderAlphaShapesPerTimeStep( const bool render,
                                             cellSet cellsAtLastTimestep,
                                             osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                                             const bool wireframe )
{
  osg::ref_ptr<osg::Group> volumeGroup = new osg::Group;
  volumeGroup->setName( "Alpha Shape" );

  osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array;

  BOOST_FOREACH( const SLineageTree *tree, cellsAtLastTimestep )
      points->push_back( osg::Vec3( tree->getX(), tree->getY(), tree->getZ() ) );

  int t = (*cellsAtLastTimestep.begin())->timeStep;

  osg::Vec4 color( 0., 1., 0., 1. );

  if( wireframe )
    color = osg::Vec4( 0., 0., 0., 1. );

  boost::shared_ptr<SAlphaShape> tri = boost::shared_ptr<SAlphaShape>(
        new SAlphaShape( points, t, color, render, wireframe ) );

  volumeGroup->addChild( tri->getTriangles() );

  if( render )
    addToThisGroup->addChild( volumeGroup );
}

// ---------------------------------------------------------------------

void computeAndRenderConvexHullsPerTimeStep( const bool render,
                                             cellSet cellsAtLastTimestep,
                                             osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                                             const bool wireframe )
{
  osg::ref_ptr<osg::Group> volumeGroup = new osg::Group;
  volumeGroup->setName( "Convex Hull" );

  osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array;

  BOOST_FOREACH( const SLineageTree *tree, cellsAtLastTimestep )
      points->push_back( osg::Vec3( tree->getX(), tree->getY(), tree->getZ() ) );

  int t = (*cellsAtLastTimestep.begin())->timeStep;

  osg::Vec4 color( 0., 1., 0., 1. );

  if( wireframe )
    color = osg::Vec4( 0., 0., 0., 1. );

  boost::shared_ptr<SConvexHull> tri = boost::shared_ptr<SConvexHull>(
        new SConvexHull( points, t, color, render, wireframe ) );

  volumeGroup->addChild( tri->getTriangles() );

  if( render )
    addToThisGroup->addChild( volumeGroup );
}

// ---------------------------------------------------------------------

void renderArrow( const osg::Vec3 &p1,
                  const osg::Vec3 &normal,
                  osg::ref_ptr<osg::Group> timeStepGroup,
                  const osg::Vec4 &color,
                  const double arrowParam1,
                  const double arrowParam2,
                  const double arrowLength )
{
  osg::Vec3 dir = normal;
  osg::ref_ptr<osg::Geode> arrow = new SArrow3DCylindrical( arrowParam1, arrowParam2, color );

  osg::Vec3 p2 = p1 + dir * arrowLength;

  osg::Vec3 n = dir ^ -dir;

  osg::MatrixTransform *mt = new osg::MatrixTransform(
        SArrow3DAngled::getTransformationMatrix( p2, p1, n ) );
  mt->addChild( arrow );
  timeStepGroup->addChild( mt );
}

}
