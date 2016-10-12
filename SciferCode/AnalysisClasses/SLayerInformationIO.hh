#ifndef SLAYERINFORMATIONIO_HH
#define SLAYERINFORMATIONIO_HH

/**
  @file   SLayerInformationIO.hh
  @brief  Contains class for I/O operations of layer information
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include "SCellLayers.hh"

#include <string>

typedef std::vector< std::map<std::size_t, std::vector<osg::Vec3> > > CellWalls;

namespace SLayerInformationIO
{
  void writeVolumes( const std::string &filename,
                     const boost::shared_ptr<cellLayerVector> &cellLayersPerTimeStep,
                     const std::vector< std::vector<double> > &layerVolumes,
                     const std::size_t numLayers );

  void writeVolume( const std::string &filename,
                    const std::vector< std::pair<double,std::size_t> > &volumes );

  void writeCellPositions( const std::string &filename,
                           const boost::shared_ptr<SCellLayers> cellLayers );

  // WARNING: This method should only be called if no cell layer
  // information was already appended in the data as a column else it will
  // erase the division type column or anything else at the back
  void writeLayerInformation( const std::string &filename,
                              const boost::shared_ptr<cellLayerVector> cellLayers );

  // WARNING: This method should only be called if some other column
  // information is appended after the layer information (which is
  // for example the division type column)
  void overwriteLayerInformation( const std::string &filename,
                                  const boost::shared_ptr<cellLayerVector> cellLayers );

  // This method is only valid to be called if the raw data includes
  // the information of each individual time step. If only start and
  // end positions are given then the method 'readLayerInformation'
  // have to be called
  void readLayerInformationPerTimeStep( const std::string &filename,
                                        boost::shared_ptr<SCellLayers> cellLayers,
                                        const int layerColumnIndex );

  void readLayerInformation( const std::string &filename,
                             boost::shared_ptr<SCellLayers> cellLayers,
                             const int layerColumnIndex );

  void generateDivisionScheme( const boost::shared_ptr<SAdvancedLineageTreeCollection> lineageTrees,
                               const NodeFeatureInfo &cellLayers,
                               NodeFeatureInfo &divisionScheme,
                               const osg::Vec3 &sideViewNormal,
                               const double radialAngleThreshold,
                               const osg::Matrix &rotMat );

  void generateDivisionScheme( const boost::shared_ptr<SAdvancedLineageTreeCollection> lineageTrees,
                               const boost::shared_ptr<SCellLayers> &cellLayers,
                               NodeFeatureInfo &divisionScheme,
                               const osg::Vec3 &sideViewNormal,
                               const double radialAngleThreshold,
                               const osg::Matrix &rotMat );

  // if this method is chosen than the division scheme is determined based only on the layer information
  // which means that only anticlinal and periclinal divisions are distinguished; this method is usefull for
  // 2D data like it is generated in the 2D model case
  void generatePartialDivisionScheme( const boost::shared_ptr<SAdvancedLineageTreeCollection> lineageTrees,
                                      const NodeFeatureInfo &cellLayers,
                                      NodeFeatureInfo &divisionScheme );

  void writeDivisionScheme( const std::string &filename,
                            const boost::shared_ptr<cellTimeVector> cTV,
                            const NodeFeatureInfo &divisionScheme );

  void overwriteDivisionScheme( const std::string &filename,
                                const NodeFeatureInfo &divisionScheme );

  void readDivisionScheme( const std::string &filename,
                           const std::size_t numTimeSteps,
                           NodeFeatureInfo &divisionScheme );

  std::string readDivisionType( const std::string &line );

  // assignment of layer value for the case of having only start and end
  // positions of the cells
  int assignLayerValue( const std::size_t tStart, const std::size_t tEnd,
                        const int id, const std::string &line,
                        const boost::shared_ptr<cellTimeVector> ctv,
                        boost::shared_ptr<cellLayerVector> clv,
                        const int layerColumnIndex );

  // assignment of layer value for the case of having all time steps of the cells
  int assignLayerValue( const int id,
                        const std::size_t timeStep,
                        const std::string &line,
                        const boost::shared_ptr<cellTimeVector> ctv,
                        boost::shared_ptr<cellLayerVector> clv,
                        const int layerColumnIndex );

  // read cell wall information if available
  void readCellWallInformation( const std::string &filename,
                                CellWalls &cellWalls );
}

#endif // SLAYERINFORMATIONIO_HH
