----------------------------------util.ns.CreateApproxSpace----------------------------------------------
--
--   Lua - Script to compute the cylinder problem
--
--	This script sets up a problem for the Navier-Stokes discretization
--	and solves the cylinder problem
--
--   Author: Josef Dubsky, Andreas Vogel
--
--------------------------------------------------------------------------------

PluginRequired("NavierStokes")

ug_load_script("ug_util.lua")
ug_load_script("util/domain_disc_util.lua")
ug_load_script("navier_stokes_util.lua")
ug_load_script("util/conv_rates_static.lua")

dim 		= util.GetParamNumber("-dim", 2, "world dimension")
numRefs 	= util.GetParamNumber("-numRefs", 2, "number of grid refinements")
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0, "number of prerefinements (parallel)")
bConvRates  = util.HasParamOption("-convRate", "compute convergence rates")
bBenchmarkRates = util.HasParamOption("-benchRate", "compute benchmark rates")

bStokes 	= util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed")
bNoLaplace 	= util.HasParamOption("-nolaplace", "If defined, only laplace term used")
bExactJac 	= util.HasParamOption("-exactjac", "If defined, exact jacobian used")
bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "full", "Upwind type")
stab        = util.GetParam("-stab", "flow", "Stabilization type")
diffLength  = util.GetParam("-difflength", "cor", "Diffusion length type")

porder = 1
vorder = 1
discType = "fe"

local Viscosity	= 1e-3
local Um = 0.3
if dim == 3 then Um = 0.45 end
local H = 0.41
local L = 0.1
local Umean2 = math.pow(2/3*Um, 2)

local ref = {}
ref.CD = 5.57953523384
ref.CL = 0.010618948146
ref.DeltaP = 0.11752016697

if 	dim == 2 then 
	gridName = util.GetParam("-grid", "grids/cylinder.ugx")
	--gridName = util.GetParam("-grid", "grids/box.ugx")
	--gridName = util.GetParam("-grid", "grids/double-arrow-small.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_tri.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_box_tri_fine.ugx")
	--gridName = util.GetParam("-grid", "grids/cylinder_rotate_box_tri_fine.ugx")
elseif dim == 3 then
	gridName = util.GetParam("-grid", "grids/cylinder3d.ugx")
--	gridName = util.GetParam("-grid", "grids/cylinder3d_fine.ugx")
else print("Choosen Dimension not supported. Exiting."); exit(); end


-- Lets write some info about the choosen parameter
print(" Choosen Parater:")
print("    dim              = " .. dim)
print("    numTotalRefs     = " .. numRefs)
print("    numPreRefs       = " .. numPreRefs)
print("    grid             = " .. gridName)
print("    porder           = " .. porder)
print("    vorder           = " .. vorder)
print("    type             = " .. discType)
print("    only stokes      = " .. tostring(bStokes))
print("    no laplace       = " .. tostring(bNoLaplace))
print("    exact jacobian   = " .. tostring(bExactJac))
print("    peclet blend     = " .. tostring(bPecletBlend))
print("    upwind           = " .. upwind)
print("    stab             = " .. stab)
print("    diffLength       = " .. diffLength)

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

function CreateDomain()

	InitUG(dim, AlgebraType("CPU", 1))
	local dom = Domain()
	LoadDomain(dom, gridName)
	
	
	-- NOTE: Projector creation in script-code is deprecated. Instead one should
	--		 add projectors to individual subsets directly in ProMesh.
	if     dim == 2 then 
		ProjectVerticesToSphere(dom, {0.2, 0.2}, 0.05, 0.001)
		falloffProjector = SphereProjector(MakeVec(0.2, 0.2, 0), 0.1, 0.15)
	elseif dim == 3 then 
		falloffProjector = CylinderProjector(MakeVec(0.5, 0.2, 0.0), MakeVec(0, 0, 1), 0.04, 0.1)
	end

	local projHandler = ProjectionHandler(dom:subset_handler())
	dom:set_refinement_projector(projHandler)

	projHandler:set_projector("Inner", falloffProjector)
	projHandler:set_projector("CylinderWall", falloffProjector)
	if dim == 3 then
		projHandler:set_projector("BackWall", falloffProjector)
		projHandler:set_projector("FrontWall", falloffProjector)
	end

	-- Create a refiner instance. This is a factory method
	-- which automatically creates a parallel refiner if required.
	local refiner =  GlobalDomainRefiner(dom)
	
	write("Pre-Refining("..numPreRefs.."): ")
	for i=1,numPreRefs do write(i .. " ");	refiner:refine(); end
	write("done. Distributing...")
	if util.DistributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, wFct) == false then
		print("Error while Distributing Grid. Aborting.")
		exit();
	end
	write(" done. Post-Refining("..(numRefs-numPreRefs).."): ")
	for i=numPreRefs+1,numRefs do refiner:refine(); write(i-numPreRefs .. " "); end
	write("done.\n")
	
	--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "grid_p"..ProcRank()..".ugx", 0.5)
	
	return dom
end

function CreateApproxSpace(dom, discType, vorder, porder)

	local approxSpace = util.ns.CreateApproxSpace(dom, discType, vorder, porder)
	
	-- print statistic on the distributed dofs
	--approxSpace:init_levels()
	approxSpace:init_top_surface()
	approxSpace:print_statistic()
	--approxSpace:print_local_dof_statistic(2)
	
	return approxSpace
end

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------
globalNSDisc = nil
function CreateDomainDisc(approxSpace, discType, vorder, porder)

	local FctCmp = approxSpace:names()
	NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, discType)
	NavierStokesDisc:set_exact_jacobian(bExactJac)
	NavierStokesDisc:set_stokes(bStokes)
	NavierStokesDisc:set_laplace( not(bNoLaplace) )
	NavierStokesDisc:set_kinematic_viscosity( Viscosity );
	globalNSDisc = NavierStokesDisc
				
	local porder = approxSpace:lfeid(dim):order()
	local vorder = approxSpace:lfeid(0):order()
	
	--upwind if available
	if discType == "fv1" or discType == "fvcr" then
		NavierStokesDisc:set_upwind(upwind)
		NavierStokesDisc:set_peclet_blend(bPecletBlend)
	end
	
	-- fv1 must be stablilized
	if discType == "fv1" then
		NavierStokesDisc:set_stabilization(stab, diffLength)
		NavierStokesDisc:set_pac_upwind(true)
	end
	
	-- fe must be stabilized for (Pk, Pk) space
	if discType == "fe" and porder == vorder then
		NavierStokesDisc:set_stabilization(3)
	end
	if discType == "fe" then
		NavierStokesDisc:set_quad_order(math.pow(vorder, dim)+2)
	end
	if discType == "fv" then
		NavierStokesDisc:set_quad_order(math.pow(vorder, dim)+2)
	end
	
	-- setup Outlet
	--OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
	--OutletDisc:add("Outlet")
	
	-- setup Inlet
	function inletVel2d(x, y, t)
		return 4 * Um * y * (H-y) / (H*H), 0.0
	end
	function inletVel3d(x, y, z, t)
		return 16 * Um * y * z * (H-y) * (H-z) / (H*H*H*H), 0.0, 0.0
	end
	InletDisc = NavierStokesInflow(NavierStokesDisc)
	InletDisc:add("inletVel"..dim.."d", "Inlet, Outlet")
	
	--setup Walles
	WallDisc = NavierStokesWall(NavierStokesDisc)
	if dim == 2 then
		WallDisc:add("UpperWall,LowerWall,CylinderWall")
	elseif dim == 3 then
		WallDisc:add("UpperWall,LowerWall,CylinderWall,FrontWall,BackWall")	
	end
	
	-- Finally we create the discretization object which combines all the
	-- separate discretizations into one domain discretization.
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(NavierStokesDisc)
	domainDisc:add(InletDisc)
	domainDisc:add(WallDisc)
	--domainDisc:add(OutletDisc)
	
	return domainDisc
end

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------

function CreateSolver(approxSpace, discType, p)

	local base = LU()
	
	local smoother = nil
	if discType == "fvcr" or discType == "fecr" then 
		smoother = ComponentGaussSeidel(0.1, {"p"}, {1,2}, {1})
	elseif discType == "fv1" then 
		smoother = ILU()
		smoother:set_damp(0.7)
	elseif discType == "fe" and porder == vorder then
		smoother = ILU()
		smoother:set_damp(0.7)
	else
		smoother = ComponentGaussSeidel(0.1, {"p"}, {0}, {1})
	end
	
	local smooth = util.smooth.parseParams()
	smoother = util.smooth.create(smooth)

	local numPreSmooth, numPostSmooth, baseLev, cycle, bRAP = util.gmg.parseParams()
	local gmg = util.gmg.create(approxSpace, smoother, numPreSmooth, numPostSmooth,
							 cycle, base, baseLev, bRAP)
	gmg:add_prolongation_post_process(AverageComponent("p"))
	-- transfer = StdTransfer()
	-- transfer:enable_p1_lagrange_optimization(false)
	-- gmg:set_transfer(transfer)

	local sol = util.solver.parseParams()
	local solver = util.solver.create(sol, gmg)
	if bStokes then
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-99, true))
	else 
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-2, true))	
	end
		
	local convCheck = ConvCheck(500, 1e-11, 1e-99, true)
	
	local newtonSolver = NewtonSolver()
	newtonSolver:set_linear_solver(solver)
	newtonSolver:set_convergence_check(convCheck)
	newtonSolver:set_line_search(StandardLineSearch(10, 1.0, 0.9, true, true))
	--newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))
	
	return newtonSolver
end

function ComputeNonLinearSolution(u, domainDisc, solver)

	util.rates.static.StdComputeNonLinearSolution(u, domainDisc, solver)
	AdjustMeanValue(u, "p")
end

--------------------------------------------------------------------------------
-- Run Problem
--------------------------------------------------------------------------------

local p = vorder
local dom = CreateDomain()
local approxSpace = CreateApproxSpace(dom, discType, vorder, porder)
local domainDisc = CreateDomainDisc(approxSpace, discType, p)
local solver = CreateSolver(approxSpace, discType, p)
--solver:set_debug(GridFunctionDebugWriter(approxSpace))
		
print(solver:config_string())

local u = GridFunction(approxSpace)
u:set(0)

--	ComputeNonLinearSolution(u, CreateDomainDisc(approxSpace, "fe", p), solver)
ComputeNonLinearSolution(u, domainDisc, solver)

local FctCmp = approxSpace:names()
local VelCmp = {}

for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end
vtkWriter = VTKOutput()
vtkWriter:select(VelCmp, "velocity")
vtkWriter:select("p", "pressure")
vtkWriter:print("navier_stokes_"..dim.."d", u)

if dim == 2 then
	local DL = DragLift(u, "u,v,p", "CylinderWall", "Inner", Viscosity, 1.0, p+3)
	local C_D = 2*DL[1]/(Umean2*L)
	local C_L = 2*DL[2]/(Umean2*L)

	local PEval = GlobalGridFunctionNumberData(u, "p")
	local Delta_P = PEval:evaluate_global({0.15, 0.2}) - PEval:evaluate_global( {0.25, 0.2} )

	print("p1: "..PEval:evaluate_global({0.15, 0.2}))
	print("p2: "..PEval:evaluate_global({0.25, 0.2}))

	print("C_D - ref.CD: "..string.format("%.3e", C_D - ref.CD))
	print("C_L - ref.CL: "..string.format("%.3e", C_L - ref.CL))
	print("Delta_P - ref.DeltaP: "..string.format("%.3e", Delta_P - ref.DeltaP))
end	
