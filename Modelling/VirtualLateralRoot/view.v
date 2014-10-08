[Main]
Dt: 0.0005 // how fast the shape of the end stage is reached, default 0.1
Growth: .1  // not considered in code
InitialCellNumber: 1 // initial cell number designed for VLR: 1, 2 or 8
InitialCellsOfRealData: 121211_raw // valid entries are:
// none, 120830_raw, 121204_raw_2014, 121211_raw, 130508_raw, 130607_raw
SubDivisionLevelOfCells: 1 // this value minus 1 defines the number of additional vertices per cell wall and therefore the subdivision level of a cell
ExportLineage: false // export lineage information of cells
ExportDivisionProperties: true // export division information of cells
SurfaceType: 1 // type of surface: 0 -> bezier surface, 1 -> surface based on triangulation of real data
SurfaceScale: 1 // scale factor for surface based on real data points

[View]
StepPerView: 1
BackgroundColor: 255

[Division]
DivisionArea: 1000 // 0.075 threshold size of cells before they divide
DivisionAreaRatio: 0.10 // 0.45 threshold of division area ratio in percentage, for example 0.5 means that a cell divides if its initial area has grown by 50%
UseAreaRatio: true // only use area ratio for divisions
UseCombinedAreaRatio: true // use the area ratio and the area threshold to prevent the cells becoming smaller step by step
UseWallRatio: false
DivisionWallRatio: 0.45 // 0.45 divide a cell if a wall of the cell is longer than a certain percentage of the initial length
UseDecussationDivision: false
ProbabilityOfDecussationDivision: 0.9 // probability for having a decussation division (has to be in [0, 1])
DivisionAngleThreshold: 45. // angle threshold to distinguish between anticlinal and periclinal division
CellColoringType: 0 // cell coloring type: 0 -> cells are colored based on founder cells/lineage trees; 1 -> cells are colored based on layer assignments after each periclinal division; 2 -> cells are colored based on type: interior or boundary

[Tissue]
DivisionAlgorithm: ShortWall // ClosestWall ShortWall ClosestMid
CellPinch: 30.1 // real data point: 10.75, default: 0.75 displacement of the position of the newly inserted vertexes
CellMaxPinch: 3.1 // 0.1 max displacement of the position of the newly inserted vertexes
CellMaxArea: 1 // not considered in code
CellWallMin: 0.1 // for SW, CW, CM
StrictCellWallMin: false // for SW, CW, CM
CellWallSample: 0.5 // for SW; not sure but influences the ratio of how a division wall is generated, default 0.005
OrderCellDivision: true  // not considered in code
CellSampleDx: 0.01//0.000000001 // for SW, default 0.01

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
