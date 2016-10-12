#ifndef __SLoadArabidopsisLineages_hh
#define __SLoadArabidopsisLineages_hh

#include "kernel/SAlgorithm.hh"
#include "dataSet/SAdvancedLineageTreeCollection.hh"

/**
  * @class Algorithm for loading 3D cell lineages from a text file
  * The data is given in a table in "Alexis file format". Radius, precursors and
  * color are ignored for now.
  */
class SLoadArabidopsisLineages : public SAlgorithm
{
public:
  SLoadArabidopsisLineages();
  ~SLoadArabidopsisLineages();

  bool canUndo();
  void undo();

  std::string getMenuEntry() const;
  std::string getName() const;
  static  std::string algoName();
    
  void execute( void* );

protected:
  void initProfile();

  /// Routine loading the data from the file into an SLineageTreeCollection
  void loadLineages();
  /// \brief Sort children according to x-coordinate.
  void sortChildren( SLineageTree *root );

  /// get cell file from current line
  int getCellFile( const std::string &line );

  int getSequenceLayerValue( const std::string &line );

  /// get layer value from current line
  int getLayerValue( const std::string &line );

  /// get division type from current line
  int getDivisionType( const std::string &line );

  /// check data content which is at the moment only the layer
  /// and division type column since these can be overwritten
  /// in an algorithm
  void checkDataContent( const std::string &line );

  /// interpolation of cell movements between cell divisions
  /// or root and a cell division
  void interpolateCellMovements( const SArray &p1, const SArray &p2,
                                 const int prTime, const int time,
                                 SLineageTree *prNode, SLineageTree *node,
                                 const int id, const int treeId );

  std::string _filename;
  boost::shared_ptr<SAdvancedLineageTreeCollection> _data;
  bool _want2Dcopy;
  bool _applyRotations;

  /// true if cell movements between cell divisions
  /// should be interpolated
  bool _interpolateMovements;

  /// store the index of the TrackGroup column when traversing a single line from backwards
  unsigned int _indexTrackGroup;

  /// store the index of the Layer column when traversing a single line from backwards
  unsigned int _indexLayer;

  /// store the index of the DivisionType column when traversing a single line from backwards
  unsigned int _indexDivisionType;

private:
  /**
   * Check (and correct) the time step order in the given cell object
   * @param bottom Leaf node within a (linear) subtree of nodes with the same
   *        cell id
   * @param lnum Number of the line in which bottom is defined. For user info.
   */
  void correctObjectHistory( SLineageTree* bottom, unsigned int lnum );
};

#endif
