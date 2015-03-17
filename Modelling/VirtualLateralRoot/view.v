[Main]
Dt: 0.002 // 0.002 // how fast the shape of the end stage is reached, default 0.1
Growth: .1  // not considered in code
InitialCellNumber: 2 // initial cell number designed for VLR: 1, 2, 6 or 8
InitialCellsOfRealData: none // valid entries are:
// none, 120830_raw, 121204_raw_2014, 121211_raw, 130508_raw, 130607_raw, Average2, Average4, Average6, Average8
SubDivisionLevelOfCells: 1 // this value minus 1 defines the number of additional vertices per cell wall and therefore the subdivision level of a cell
ExportLineage: true // export lineage information of cells
ExportDivisionProperties: true // export division information of cells
SurfaceType: 0 // type of surface: 0 -> bezier surface, 1 -> surface based on triangulation of real data
BezierGrowthSurface: true // if the surface type is bezier then it can be chosen to choose the generated bezier surface that includes the growth tensor information
SurfaceScale: 1 // scale factor for surface based on real data points
UseAutomaticContourPoints: false // use contour points that are generated automatically
ForceInitialSituation: true // force the initial situation beginning with two founder cells that first divide anticlinally in a 1/3:2/3 ratio and afterwards periclinally resulting in 6 cells
SmootherCells: false // TODO: if this is true then there are two additional cell positions between two subsequent time steps that are interpolated linearly

[View]
StepPerView: 1
BackgroundColor: 255

[Division]
DivisionArea: 1100 // 1100 ... 0.025 threshold size of cells before they divide
DivisionAreaRatio: 0.62 // 0.145 ... 0.45 threshold of division area ratio in percentage, for example 0.5 means that a cell divides if its initial area has grown by 50%
UseAreaRatio: true // only use area ratio for divisions
UseCombinedAreaRatio: true // use the area ratio and the area threshold to prevent the cells becoming smaller step by step
UseWallRatio: false
DivisionWallRatio: 0.45 // 0.45 divide a cell if a wall of the cell is longer than a certain percentage of the initial length
DivisionType: Decussation // Decussation PerToGrowth ; defines the type of division which then replaces the chosen type of division set by the variable DivisionAlgorithm below
ProbabilityOfDecussationDivision: 0.9 // probability for having a decussation division (has to be in [0, 1])
DivisionAngleThreshold: 45. // angle threshold to distinguish between anticlinal and periclinal division
CellColoringType: 0 // cell coloring type: 0 -> cells are colored based on founder cells/lineage trees; 1 -> cells are colored based on layer assignments after each periclinal division; 2 -> cells are colored based on type: interior or boundary
FirstDivisionsAreaRatio: 0.001 // first division ratio if forced initial situation is used
SecondDivisionsAreaRatio: 0.3 // second division ratio if forced initial situation is used
TimeDelay: 20 // time delay for next division after the six cell stage if forced initial situation is used

[Tissue]
DivisionAlgorithm: ShortWall // ClosestWall ShortWall ClosestMid
CellPinch: 0.75 // real data point: 5.1, default: 0.75 displacement of the position of the newly inserted vertexes 0.155
CellMaxPinch: 0.1 // 1.5 ... 0.1 ... 0.155 max displacement of the position of the newly inserted vertexes
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
