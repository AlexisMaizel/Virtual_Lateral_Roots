#ifndef ModelHeader_HH
#define ModelHeader_HH

/**
  @file   ModelHeader.h
  @brief  Contains all relevant headers and typedefs for the model
  @author Jens Fangerau <jens.fangerau@iwr.uni-heidelberg.de>
*/

#include <vve.h>

#include <util/parms.h>
#include <util/palette.h>
#include <util/assert.h>

#include <geometry/geometry.h>
#include <geometry/area.h>
#include <geometry/intersection.h>

#include <math.h>
#include <fstream>
#include <limits>

#include <QTextStream>
#include <stdio.h>

#include "bezier.h"
#include "surface.h"
#include "RealSurface.h"

static QTextStream out(stdout);

using util::norm;

namespace SurfaceType
{
  enum type { BEZIER, REALDATA };
};

namespace DivisionType
{
  // the numbers are chosen depending on the index of colors in the color palette
  enum type { ANTICLINAL=3, PERICLINAL=1, RADIAL=2, NONE=5 };
};

struct JunctionContent
{
  SurfacePoint sp;
  std::size_t id;
  TrianglePoint tp;
};

struct CellContent
{
  SurfacePoint sp;
  TrianglePoint tp;
  std::size_t id;
  std::size_t parentId;
  std::size_t treeId;
  std::size_t timeStep;
  Point3d center;
  // initial area of cell right after division
  double initialArea;
  // changing area of cell
  double area;
  // set of precursor cell ids
  std::set<std::size_t> precursors;
  // previous division angle according to x-axis
  double previousAngle;
  // division angle according to x-axis
  double angle;
  // longest wall length right after division
  double initialLongestWallLength;
  // longest wall length for division based on longest wall
  double longestWallLength;
  // origin division type of cell
  // the numbers are chosen based on the color map in scifer
  // 3 -> anticlinal -> red
  // 1 -> periclinal -> green
  // 2 -> radial (which is not used at the moment since this makes no sense in 2D) -> blue
  // 5 -> none which is only valid for the initial cells at the beginning -> cyan
  DivisionType::type divType;
  // layer value for current cell
  // the numbers are chosen based on the color map in scifer
  // layer 1 -> yellow (9)
  // layer 2 -> blue (10)
  // layer 3 -> magenta (11)
  // layer 4 -> red (12) etc.
  std::size_t layerValue;
  // corresponding cell cycle
  std::size_t cellCycle;
  // previous division direction of cell
  Point3d previousDivDir;
  // current division direction of cell
  Point3d divDir;
  // convex hull of cell
  std::vector<Point2d> convexHull;
};

struct WallContent
{
};

class MyModel;

typedef std::set< std::pair<std::size_t,std::size_t> > idPairSet;

typedef tissue::Tissue<MyModel, CellContent, JunctionContent, WallContent, graph::_EmptyEdgeContent, graph::_EmptyEdgeContent, graph::_EmptyEdgeContent, false> MyTissue;

typedef MyTissue::junction_cell_edge junction_cell_edge;
typedef MyTissue::const_junction_cell_edge const_junction_cell_edge;
typedef MyTissue::junction junction;
typedef MyTissue::wall wall;
typedef MyTissue::const_wall const_wall;
typedef MyTissue::wall_arc wall_arc;
typedef MyTissue::cell cell;
typedef MyTissue::cell_arc cell_arc;
typedef MyTissue::cell_edge cell_edge;
typedef MyTissue::const_cell_edge const_cell_edge;
typedef MyTissue::cell_junction_edge cell_junction_edge;
typedef MyTissue::const_cell_junction_edge const_cell_junction_edge;
typedef MyTissue::wall_graph wall_graph;
typedef MyTissue::complex_graph complex_graph;
typedef MyTissue::cell_graph cell_graph;

typedef util::Colorf Colorf;

#endif // ModelHeader_HH