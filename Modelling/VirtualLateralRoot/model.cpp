#include <vve.h>

#include <util/parms.h>
#include <util/palette.h>
#include <util/assert.h>

#include <geometry/geometry.h>
#include <geometry/area.h>
#include <geometry/intersection.h>

#include <math.h>
#include <fstream>

#include <QTextStream>
#include <stdio.h>

#include "bezier.h"
#include "surface.h"

static QTextStream out(stdout);

using util::norm;

struct JunctionContent
{
  SurfacePoint sp;
};

struct CellContent
{
  SurfacePoint sp;
	std::size_t id;
  std::size_t treeId;
  std::size_t timeStep;
  double initialArea;
  double area;
  std::set<std::size_t> precursors;
  Point3d center;
};

struct WallContent
{
};

class MyModel;

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

class MyModel : public Model 
{
  Q_OBJECT
public: 
  util::Parms parms;
  util::Palette palette;
	Surface lateralRoot;
  MyTissue T;

  double dt;
  double divisionArea;
  bool useRatio;
  double divisionAreaRatio;
  int bgColor;
  int cellInitWalls;
  int stepPerView;
  int initialConstellation;
  
  std::size_t _idCounter;
  std::size_t _time;
  std::string _fileName;

  //----------------------------------------------------------------
  
  void readParms()
  {
    // read the parameters here
    parms("Main", "Dt", dt);
    parms("Main", "DivisionArea", divisionArea);
    parms("Main", "UseRatio", useRatio);
    parms("Main", "DivisionAreaRatio", divisionAreaRatio);
    parms("Main", "CellInitWalls", cellInitWalls);
    parms("Main", "InitialConstellation", initialConstellation);

    parms("View", "StepPerView", stepPerView);
    parms("View", "BackgroundColor", bgColor);


    T.readParms(parms, "Tissue");
    T.readViewParms(parms, "TissueView");
  }

  //----------------------------------------------------------------
  
  // Here, reread the files when they are modified
  void modifiedFiles( const std::set<std::string>& filenames )
  {
    forall(const std::string& fn, filenames)
    {
      if(fn == "pal.map")
        palette.reread();
      else if(fn == "view.v")
        readParms();
    }
  }

  //----------------------------------------------------------------
  
  MyModel(QObject *parent) : Model(parent), parms("view.v"), palette("pal.map"), 
    lateralRoot(parms, "Surface"), T(palette, this),
    _idCounter(1), _time(1), _fileName( "/tmp/model.csv" )
  {
    readParms();
    // Registering the configuration files
    registerFile("pal.map");
    registerFile("view.v");

    this->initExportFile( _fileName );
    
    lateralRoot.GrowStep(0);
    if( initialConstellation == 0 )
      this->initOneCell();
    else if( initialConstellation == 1 )
      this->initLateralRoot();
    
    setStatus();
  }

  //----------------------------------------------------------------
  
  void initOneCell()
  {
    std::vector<junction> vs;

    for(int i = 0; i < cellInitWalls; i++) {
      junction v; 
      T.W.insert(v);
      vs.push_back(v);

      double initu, initv, s = (double)i/cellInitWalls;
      if(s < .25) {
        initu = 0;
        initv = s * 4;
      } else if(s < .5) {
        initu = (s - .25) * 4;
        initv = 1.0;
      } else if(s < .75) {
        initu = 1.0;
        initv = 1  - (s - .5) * 4;
      } else {
        initu = 1  - (s - .75) * 4;
        initv = 0;
      }
      lateralRoot.InitPoint(v->sp, initu, initv);
    }

    cell c;
    c->treeId = 1;
    c->id = _idCounter;
    c->timeStep = _time;
    T.addCell(c, vs);

    std::vector<Point3d> polygon;
    Point3d center;
    forall(const junction& j, T.S.neighbors(c))
    {
      polygon.push_back(j->sp.Pos());
      center += j->sp.Pos();
    }
    center /= polygon.size();
    lateralRoot.SetPoint(c->sp, c->sp, center);
    
    // store initial area for current cell
    c->center = center;
    c->initialArea = geometry::polygonArea(polygon);
    c->area = c->initialArea;
    _idCounter++;
  }
  
  //----------------------------------------------------------------
  
  void initLateralRoot()
  {
    // render eight cells at the beginning for which each pair shares an
    // area ratio of 2:1. The bigger cell will always be located at the left
    // and right boundary
    
    std::size_t lineageCounter = 1;
    
    // init the inner smaller cells
    for( std::size_t w = 0; w < 2; w++ )
      for( std::size_t h = 0; h < 2; h++ )
      {
        double u = 1./3. + w*1./6.;
        double v = 0. + h*1./2.;
        this->generateCell( std::make_pair( u, v ),
                            std::make_pair( 1./6., 1./2. ),
                            lineageCounter );
        
        lineageCounter++;
      }
    
    // then the outer bigger cells
    for( std::size_t w = 0; w < 2; w++ )
      for( std::size_t h = 0; h < 2; h++ )
      {
        double u = 0. + w*2./3.;
        double v = 0. + h*1./2.;
        this->generateCell( std::make_pair( u, v ),
                            std::make_pair( 1./3., 1./2. ),
                            lineageCounter );
        
        lineageCounter++;
      }
  }
  //----------------------------------------------------------------
  
  // generate a cell starting at the bottom left with coord u and v with different lengths
  cell generateCell( const std::pair<double, double> &start,
                     const std::pair<double, double> &length,
                     const std::size_t treeId )
  {
    // set of junctions for the cell
    std::vector<junction> vs;
        
    double u = start.first;
    double v = start.second;
    
    for( std::size_t w = 0; w < 4; w++ )
    {
      junction j;
      // perhaps not required because the size is correct?
      //T.W.insert(j);
      double u,v;
      
      switch(w)
      {
        case 0: u = start.first; v = start.second; break;
        case 1: u = start.first; v = start.second + length.second; break;
        case 2: u = start.first + length.first; v = start.second + length.second; break;
        case 3: u = start.first + length.first; v = start.second; break;
        default: u = start.first; v = start.second; break;
      }
      
      lateralRoot.InitPoint( j->sp, u, v );
      junctionAlreadyShared( j->sp, j );
      vs.push_back(j);
    }
    
    cell c;
    c->treeId = treeId;
    c->id = _idCounter;
    c->timeStep = _time;
    
    T.addCell( c, vs );
    
    std::vector<Point3d> polygon;
    Point3d center;
    forall(const junction& j, T.S.neighbors(c))
    {
      polygon.push_back(j->sp.Pos());
      center += j->sp.Pos();
    }
    center /= polygon.size();
    c->center = center;
    
    // store initial area for current cell
    c->initialArea = geometry::polygonArea(polygon);
    c->area = c->initialArea;
    lateralRoot.SetPoint(c->sp, c->sp, center);
    
    // afterwards increment the id counter
    _idCounter++;
  }
  
  //----------------------------------------------------------------
  
  void junctionAlreadyShared( SurfacePoint sp, junction &js )
  {
    forall( const cell& c, T.C )
    {
      forall(const junction& j, T.S.neighbors(c))
      {
        if( j->sp.Pos() == sp.Pos() )
        {
          js = j;
          return;
        }
      }
    }
  }

  //----------------------------------------------------------------
  
  void initExportFile( const std::string &filename )
  {
    std::ofstream out( filename.c_str(), std::ofstream::out );

    // rotation information
    out << "0 0 0\n";
    // center of root
    out << "0 -2 0\n";
    // dimension
    out << "3 2 0\n";
    // header
    out << "ObjectID X Y Z Timepoint Radius Precursors Color Lineage TrackID TrackColor TrackGroup\n";
  }
  
  //----------------------------------------------------------------
  
  void exportLineageInformation( const std::string &filename,
                                 const cell& c )
  {
    std::ofstream out( filename.c_str(), std::ofstream::out | std::ofstream::app );
    
    std::vector<Point3d> polygon;
    Point3d center;
    forall(const junction& j, T.S.neighbors(c))
    {
      polygon.push_back(j->sp.Pos());
      center += j->sp.Pos();
    }
    center /= polygon.size();
    c->center = center;
    c->area = geometry::polygonArea(polygon);
    c->timeStep = _time;
    
    out << c->id << " "
        << c->center.i() << " "
        << c->center.j() << " "
        << c->center.k() << " "
        << c->timeStep << " "
        << c->area << " "
        << "\"{"; 
    
    std::size_t counter = 1;
    for( std::set<std::size_t>::iterator setIter = c->precursors.begin();
         setIter != c->precursors.end(); setIter++, counter++ )
    {
      if( counter != c->precursors.size() )
        out << *setIter << ", ";
      else
        out << *setIter << "}\" ";
    }
    
    if( c->precursors.size() == 0 )
      out << "}\" ";
        
    out << "Color "
        << c->treeId << " "
        << "TrackId "
        << "TrackColor "
        << "0\n";
        
    out.close();
  }
  
  //----------------------------------------------------------------
  
  void setStatus()
  {
		double time = lateralRoot.GetTime();
    setStatusMessage(QString("# vertices: %1 - # cells: %2 - step: %3 - timeStep: %4").arg(T.W.size()).arg(T.C.size()).arg(time).arg(_time));
  }

  //----------------------------------------------------------------
  
  void updateFromOld(const cell& cl, const cell& cr, const cell& c,
                     const MyTissue::division_data& ddata, MyTissue&)
  {
		// inherit treeid information
		cl->treeId = cr->treeId = c->treeId;
    // and set new cell ids for the daughter cells
    cl->id = _idCounter;
    _idCounter++;
    cr->id = _idCounter;
    _idCounter++;
    // set the current time steps
    cl->timeStep = _time;
    cr->timeStep = _time;
    // insert the new initial areas
    // left cell
    std::vector<Point3d> polygon;
    forall(const junction& j, T.S.neighbors(cl))
      polygon.push_back(j->sp.Pos());
    
    cl->initialArea = geometry::polygonArea(polygon);
    cl->area = cl->initialArea;
    
    // right cell
    polygon.clear();
    forall(const junction& j, T.S.neighbors(cr))
      polygon.push_back(j->sp.Pos());
    
    cr->initialArea = geometry::polygonArea(polygon);
    cr->area = cr->initialArea;
    
    // update precursors
    cl->precursors.insert( c->id );
    cr->precursors.insert( c->id );
    
    for( std::set<std::size_t>::iterator setIter = c->precursors.begin();
         setIter != c->precursors.end(); setIter++ )
    {
      cl->precursors.insert( *setIter );
      cr->precursors.insert( *setIter );
    }
  }

  //----------------------------------------------------------------
  
  bool step_divisions()
  {
    // Find cells to be divided
    std::list<cell> to_divide;
    forall(const cell& c, T.C)
    {
      std::vector<Point3d> polygon;
      Point3d center;
      forall(const junction& j, T.S.neighbors(c))
      {
        polygon.push_back(j->sp.Pos());
        center += j->sp.Pos();
      }
      center /= polygon.size();
      lateralRoot.SetPoint(c->sp, c->sp, center);
      double a = geometry::polygonArea(polygon);
      if( useRatio )
      {
        // only apply the division based on ratio with at least eight cells
        if( c->id > 7 )
        {
          // divide cells if their area size has increased by a certain
          // percentage
          double initialArea = c->initialArea;
          initialArea += initialArea*divisionAreaRatio;
          if( a > initialArea )
            to_divide.push_back(c);
        }
        else
        {
          if(a > divisionArea)
            to_divide.push_back(c);
        }
      }
      else
      {
        if(a > divisionArea)
          to_divide.push_back(c);
      }
    }

    // Divide the cells
    forall(const cell& c, to_divide)
    {
      T.divideCell(c);
    }

    return !to_divide.empty();
  }

  //----------------------------------------------------------------
  
  void step_tracking()
  {
    if( _time <= 250 )
    {
      forall(const cell& c, T.C)
      {
        this->exportLineageInformation( _fileName, c );
      }
      _time++;
    }
  }
  
  //----------------------------------------------------------------
  
  void step_growth()
  {
	  lateralRoot.GrowStep(dt);
    forall(const junction& v, T.W)
		  lateralRoot.GetPos(v->sp);
  }

  //----------------------------------------------------------------
  
  void step()
  {
    for(int i = 0 ; i < stepPerView ; ++i)
    {
      this->step_tracking();
      this->step_divisions();
      this->step_growth();
    }
    this->setStatus();
  }

  //----------------------------------------------------------------
  
	void initDraw(Viewer* viewer)
	{
    viewer->setSceneBoundingBox(Vec(-10.0, -10.0, -1.0), Vec(10.0, 10.0, 1.0));
	}

	//----------------------------------------------------------------
	
  void preDraw()
  {
    util::Palette::Color bg = palette.getColor(bgColor);
    glClearColor(bg.r(), bg.g(), bg.b(), 1);
    T.preDraw();
  }

  //----------------------------------------------------------------
  
  void postDraw()
  {
    T.postDraw();
  }

  //----------------------------------------------------------------
  
  void draw(Viewer* viewer)
  {
    forall(const cell& c, T.C)
    {
      //T.drawCell(c, .5);
      T.drawCell(c, this->cellColor(c), this->cellColor(c)*0.6 );
    }
  }

  //----------------------------------------------------------------
  
	// color for inner cells
  Colorf cellColor(const cell& c)
  {
  	return palette.getColor(c->treeId);
  }

  //----------------------------------------------------------------
  
	// color for contour of cells
  Colorf contourColor(const cell& c)
  {
    return palette.getColor(T.contourColor);
  }

  //----------------------------------------------------------------
  
  // Methods needed by the tissue
  Point3d position(const cell& c) const { return c->sp.Pos(); }
  Point3d position(const junction& c) const { return c->sp.Pos(); }

  void setPosition(const cell& c, const Point3d& p) { lateralRoot.SetPoint(c->sp, c->sp, p); }
  void setPosition(const junction& j, const Point3d& p) { lateralRoot.SetPoint(j->sp, j->sp, p); }

  void setPositionHint(const junction&, const junction&, const junction&, double) {}

  Point3d normal(const junction& ) const { return Point3d(0,0,1); }
  Point3d normal(const cell& ) const { return Point3d(0,0,1); }

};

#include "model.moc"
DEFINE_MODEL(MyModel);
