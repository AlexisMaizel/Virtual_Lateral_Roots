[Main]
Dt: 0.002 // 0.002 // how fast the shape of the end stage is reached, default 0.1
Growth: .1  // not considered in code
InitialCellNumber: 2 // initial cell number designed for VLR: 1, 2, 6 or 8
RealDataName: none // valid entries are:
// none, 120830_raw, 121204_raw_2014, 121211_raw, 130508_raw, 130607_raw, Average2, Average4, Average6, Average8
SubDivisionLevelOfCells: 10 // this value minus 1 defines the number of additional vertices per cell wall and therefore the subdivision level of a cell
BezierGrowthSurface: true // if the surface type is bezier then it can be chosen to choose the generated bezier surface that includes the growth tensor information
SurfaceScale: 1 // scale factor for surface based on real data points
UseAutomaticContourPoints: false // use contour points that are generated automatically
InitialSituationType: 0 // forced division situation at the beginning; can be 0, 1, 2, or 3 -> 0: not hardwired at the beginning, 1: force the initial situation beginning with two founder cells that divide anticlinally in a 1/3:2/3 ratio, 2: force the initial situation beginning with two founder cells that first divide anticlinally in a 1/3:2/3 ratio and afterwards periclinally resulting in 6 cells, 3: force the initial situation beginning with two founder cells that first divide anticlinally in a 1/3:2/3 ratio and afterwards anticlinally again resulting in 6 cells
CenterOfMassAfterLOD: true // compute the center of mass after applying a level of detail
Loop: false
LODThreshold: 1 // threshold for used edge criterion in LOD
AvoidTrianglesThreshold: 0 // in [0, 100] percentage distance threshold in order to avoid triangle cells; example: 20 -> 20% of total length of cell wall is the minimum distance that should be guaranteed between division line end point and junction of cell wall
LoadLastModel: false // if true load the last created model for which this variable was set to false -> this is required to rerun models that are generated with randomized parameters
OnlyGrowthInHeight: true // the idealized bezier surface changes in width and height -> setting this parameter to true results in a bezier surface only increasing in height
RenderMovies: false // if set to true then create movies from the same model (even if randomized) visualized by cells, spheres, and layercurves, all colored by layers and additionally one movie where the cells are colored by founder cells

[View]
StepPerView: 1
BackgroundColor: 255
RenderCells: true
RenderJunctions: false
RenderSpheres: false
RenderLayerCurves: false
RenderControlPoints: false
RenderBezierSurface: false
RenderCellCenter: false
RenderPCLine: false // only when PerToGrowth as division type is used

[Division]
DivisionArea: 400 // 1100 ... 0.025 threshold size of cells before they divide
DivisionAreaRatio: 0.1 // 0.145 ... 0.45 threshold of division area ratio in percentage, for example 0.5 means that a cell divides if its initial area has grown by 50%
EqualAreaRatio: 1.
UseAreaRatio: true // only use area ratio for divisions
UseCombinedAreaRatio: true // use the area ratio and the area threshold to prevent the cells becoming smaller step by step
UseWallRatio: false
DivisionWallRatio: 0.45 // 0.45 divide a cell if a wall of the cell is longer than a certain percentage of the initial length
UseAlternativeDivisionType: true
DivisionType: Besson-Dumais // Decussation PerToGrowth Energy Besson-Dumais Random RandomEqualAreas; defines the type of division which then replaces the chosen type of division set by the variable DivisionAlgorithm below
ProbabilityOfDecussationDivision: 100 // probability for having a decussation division (has to be in [0, 100])
DivisionAngleThreshold: 45. // angle threshold to distinguish between anticlinal and periclinal division
CellColoringType: 1 // cell coloring type: 0 -> cells are colored based on founder cells/lineage trees; 1 -> cells are colored based on layer assignments after each periclinal division; 2 -> cells are colored based on type: interior or boundary
FirstDivisionsAreaRatio: 0.001 // first division ratio if forced initial situation is used
SecondDivisionsAreaRatio: 0.3 // second division ratio if forced initial situation is used
TimeDelay: 20 // time delay for next division after the six cell stage if forced initial situation is used

[Tissue]
DivisionAlgorithm: ShortWall // ClosestWall ShortWall ClosestMid
CellPinch: 0.0 // real data point: 5.1, default: 0.75 displacement of the position of the newly inserted vertexes 0.155
CellMaxPinch: 0.0 // 1.5 ... 0.1 ... 0.155 max displacement of the position of the newly inserted vertexes
CellWallMin: 0.0 // for SW, CW, CM
StrictCellWallMin: false // for SW, CW, CM
CellWallSample: 0.0005 // 0.5 // for SW; not sure but influences the angle of a division wall, default 0.005
CellSampleDx: 0.001 // for SW, default 0.01

[TissueView]
CellColorBegin: 1
CellColorEnd: 2
CellColorCenter: 0.2
CellContourColor: 0
CellWallCorner: .01
DrawInsideCells: true
DrawCellBorders: true
CellBlending: 1
CellCenterBlending: 1
CellWallWidth: .003

[Surface]	
Surfaces: 2   			// Number of surfaces
SurfTimeScale: .2		// Surface time scale
SurfMaxDist: .2			// Max distance from surface for closest point search
Surface0: vlr1.s		// Surfaces	
Surface1: vlr2.s
SurfaceScale0: 1		// Scale factors
SurfaceScale1: 1
SurfaceTime0: 1			// Time until at next stage
SurfaceTime1: 6
