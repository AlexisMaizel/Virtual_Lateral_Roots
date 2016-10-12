#ifndef SCORRELATIONANALYSIS_HH
#define SCORRELATIONANALYSIS_HH

#include "dataSet/SAdvancedLineageTreeCollection.hh"

#include "graphics/SScatterplotMatrix.hh"

struct divData
{
  std::size_t divType;
  double VLRAngle;
  double distDivToCentre;
  double primAngle;
  double distDivToMainRoot;
  double surfaceAngle;
  double distDivToSurface;
};

class SCorrelationAnalysis
{
public:
  SCorrelationAnalysis( std::vector< boost::shared_ptr<SAdvancedLineageTreeCollection> > lineages,
                        const std::size_t dim );

  void readDivisionData( const std::string &fileName,
                         const std::size_t dataId );

  void generateDivisionData( const std::size_t divType,
                             const std::size_t maxTimeSteps );

  void generateGeneralData( const std::size_t maxTimeSteps );

  osg::Vec3 getPreviousDivPos( const SLineageTree *node );

  void renderScatterplotMatrix( osg::ref_ptr<osg::Group> addToThisGroup );

private:

  std::vector< boost::shared_ptr<SAdvancedLineageTreeCollection> > _lineages;

  std::vector< std::vector<double> > _data;

  std::size_t _dim;

  std::size_t _numPoints;

  std::size_t _divType;

  bool _loadData;

  std::size_t _dataType;

  std::vector< std::map< std::pair<std::size_t, std::size_t>, divData > > _readDivisionData;
};

#endif // SCORRELATIONANALYSIS_HH
