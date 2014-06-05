[Main]
Dt: .1 // how fast the shape of the end stage is reached, default 0.01
Growth: .1  // not considered in code
CellInitWalls: 4 // number of cell walls at the beginning
InitialConstellation: 0 // initial cell constellation: 0 -> starting with one cell, 1 -> starting with eight cells positioned more similar to lateral root
ExportLineage: true // export lineage information of cells

[View]
StepPerView: 1
BackgroundColor: 255

[Division]
DivisionArea: 0.06 // 0.075 threshold size of cells before they divide
UseAreaRatio: true
DivisionAreaRatio: 0.7 // 0.45 threshold of division area ratio in percentage, for example 0.5 means that a cell divides if its initial area has grown by 50%
UseWallRatio: true
DivisionWallRatio: 0.45 // 0.45 divide a cell if a wall of the cell is longer than a certain percentage of the initial length
UseDecussationDivision: true
ProbabilityOfDecussationDivision: 0.9 // probability for having a decussation division (has to be in [0, 1])
DivisionAngleThreshold: 45. // angle threshold to distinguish between anticlinal or periclinal division

[Tissue]
DivisionAlgorithm: ShortWall // ClosestWall ShortWall ClosestMid
CellPinch: 1 // 2 displacement of the position of the newly inserted vertexes
CellMaxPinch: 0.1 // 0.05 max displacement of the position of the newly inserted vertexes
CellMaxArea: 1 // not considered in code
CellWallMin: 0.0 // for SW, CW, CM
StrictCellWallMin: false // for SW, CW, CM
CellWallSample: 0.005 // for SW; not sure but influences the ratio of how a division wall is generated, default .005
OrderCellDivision: true  // not considered in code
CellSampleDx: 0.01 // for SW 

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
SurfaceTime1: 5
