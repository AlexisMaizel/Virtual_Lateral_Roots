#ifndef SSIMILARITYMEASUREGRAPHICS_HH
#define SSIMILARITYMEASUREGRAPHICS_HH

#include "SSimilarityMeasureHeader.hh"

#include "graphics/SColorMap.hh"

#include "SSliderCallback.hh"

#include "SAlphaShape.hh"
#include "SConvexHull.hh"

#include <osg/MatrixTransform>

namespace SSimilarityMeasureGraphics
{
void drawTimeLineLabels( osg::ref_ptr<osg::Geode> labelGeode,
                         const double stepSize,
                         const double maxDepth,
                         const double coordStart,
                         const double coordMax,
                         double maxXCoord = 0.,
                         const bool vertical = true,
                         const double labelSize = 10.,
                         const bool drawTimeLines = true );

void drawLine( osg::ref_ptr<osg::Geode> lineGeode,
               const osg::Vec3 &pos1, const osg::Vec3 &pos2,
               const double lineWidth, const osg::Vec4 &color );

void drawSequenceRectangle( osg::ref_ptr<osg::Geode> geode,
                            const int divType,
                            const double size,
                            const osg::Vec2 &pos );

void drawSequenceRectangle( osg::ref_ptr<osg::Geode> geode,
                            const std::string &label,
                            const double labelSize,
                            const double sizeX,
                            const double sizeY,
                            const osg::Vec3 &pos,
                            const osg::Vec4 &color );

void drawDoubleSequenceRectangle( osg::ref_ptr<osg::Geode> geode,
                                  const std::string &label1,
                                  const std::string &label2,
                                  const double labelSize,
                                  const double sizeX,
                                  const double sizeY,
                                  const osg::Vec3 &pos,
                                  const osg::Vec4 &color1,
                                  const osg::Vec4 &color2 );

void drawComparisonBox( osg::ref_ptr<osg::Geode> geode,
                        const osg::Vec3 &leftTop,
                        const osg::Vec3 &rightBottom,
                        const std::string &compLabel,
                        const std::string &appLabel,
                        const double labelSize );

void renderDivisionSequence( const SLineageTree *node,
                               const NodeFeatureInfo &cellDivisionScheme,
                               const std::size_t startTime,
                               osg::ref_ptr<osg::Group> group );

void drawComparisonRectangle( osg::ref_ptr<osg::Geode> geode,
                              const double sizeX,
                              const double sizeY,
                              const osg::Vec3 &pos,
                              double ratio,
                              const std::size_t appearances );

void addCylinderBetweenPoints( const osg::Vec3& p1, const osg::Vec3& p2,
                               double rad, const osg::Vec4& color,
                               osg::ref_ptr<osg::Group> addToThisGroup,
                               const std::string &drawableName = "Drawable" );

void addArrowBetweenPoints( const osg::Vec3& p1, const osg::Vec3& p2,
                            const osg::Vec3& normal, const osg::Vec4& color,
                            const osg::Vec4& tipColor,
                            osg::ref_ptr<osg::Group> addToThisGroup );

void addSimpleArrowBetweenPoints( const osg::Vec3& p1, const osg::Vec3& p2,
                                  const osg::Vec4& color, const osg::Vec4& tipColor,
                                  osg::Group *addToThisGroup, const double radius = 2.5 );

void renderCylindricalDivisionAnalysis( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                                        osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                                        osg::ref_ptr<SSliderCallback> sliderCallback,
                                        const std::size_t numTimeSteps,
                                        const std::size_t numCellCycles,
                                        const bool onlyMasterCellFile,
                                        const std::vector<osg::Vec3> &centresOfMass,
                                        const std::size_t sliderType );

void renderDivisions( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                      osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                      osg::ref_ptr<SSliderCallback> sliderCallback,
                      osg::ref_ptr<osg::Vec4Array> layerColors,
                      const NodeFeatureInfo &layerValues,
                      const bool renderTracks,
                      const bool onlyMasterCellFile,
                      const std::size_t divisionColorType,
                      const bool renderArrows );

void render3DCross( const osg::Vec3 &center,
                     osg::ref_ptr<osg::MatrixTransform> addToThisGroup );

void drawTriangle( const osg::ref_ptr<osg::Vec3Array> points,
                   const osg::Vec4 &color,
                   osg::ref_ptr<osg::Group> addToThisGroup );

void drawPlane( const osg::ref_ptr<osg::Vec3Array> points,
                const osg::Vec4 &color,
                osg::ref_ptr<osg::Group> addToThisGroup );

void computeAndRenderSurfacesPerTimeStep( const bool render,
                                          const bool output,
                                          const int type,
                                          osg::ref_ptr<SSliderCallback> sliderCallback,
                                          boost::shared_ptr<cellTimeVector> cTS,
                                          osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                                          const bool wireframe,
                                          const double surfaceOffset );

void computeAndRenderSurfacesPerTimeStep( const bool render,
                                          const bool output,
                                          const int type,
                                          osg::ref_ptr<SSliderCallback> sliderCallback,
                                          boost::shared_ptr<cellTimeVector> cTS,
                                          const std::vector<int> &averagedCellCycleTimeSteps,
                                          osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                                          const bool wireframe,
                                          const double surfaceOffset );

void computeAndRenderAlphaShapesPerTimeStep( const bool render,
                                             cellSet cellsAtLastTimestep,
                                             osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                                             const bool wireframe );

void computeAndRenderConvexHullsPerTimeStep( const bool render,
                                             cellSet cellsAtLastTimestep,
                                             osg::ref_ptr<osg::MatrixTransform> addToThisGroup,
                                             const bool wireframe );

void renderArrow( const osg::Vec3 &p1,
                  const osg::Vec3 &normal,
                  osg::ref_ptr<osg::Group> timeStepGroup,
                  const osg::Vec4 &color = osg::Vec4( 1., 0., 0., 1. ),
                  const double arrowParam1 = 1.8,
                  const double arrowParam2 = 0.4,
                  const double arrowLength = 28. );

}

#endif // SSIMILARITYMEASUREGRAPHICS_HH
