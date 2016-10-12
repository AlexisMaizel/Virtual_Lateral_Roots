#ifndef SVLRDATAANALYSISIO_HH
#define SVLRDATAANALYSISIO_HH

#include "SSimilarityMeasureHeader.hh"

struct Divs
{
  osg::Vec3 divPosOrig;
  osg::Vec3 divPos;
  osg::Vec3 dauPos1Orig;
  osg::Vec3 dauPos2Orig;
  osg::Vec3 dauPos1;
  osg::Vec3 dauPos2;
};

struct FileAngles
{
  int cellFile;
  double theta;
  double phi;
  std::size_t divType;
};

namespace SVLRDataAnalysisIO
{
void printDivisionInformation( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                               const std::size_t numTotalTimesteps,
                               const std::vector<osg::Vec3> &centresOfMass );

void exportCompressedDataset( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                              const NodeFeatureInfo &layerValues );

void appendDivisionProperties( boost::shared_ptr<SAdvancedLineageTreeCollection> lineages,
                               const osg::Matrix &rotMat );

void readAngleData( const std::string &fileName,
                    std::map<std::size_t, FileAngles> &angles );

void readModelAngleData( const std::string &fileName,
                         std::vector< std::map<std::size_t, FileAngles> > &angles );

void exportAngleAnalysis( const std::string &fileName,
                          const std::vector< std::map<int, double> > &averagedValues );

void readCellShapeVertices( const std::string &fileName,
                            std::vector< std::vector<osg::Vec3> > &vertices );

}

#endif // SVLRDATAANALYSISIO_HH
