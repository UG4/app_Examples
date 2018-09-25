----------------------------------------------------------
--
--   Lua - Script to perform the Elder-Problem on an adaptive grid
--
--   Author: Andreas Vogel, Sebastian Reiter
--
----------------------------------------------------------

PrintBuildConfiguration()

ug_load_script("../scripts/ug_util.lua")
ug_load_script("util/load_balancing_util.lua")
ug_load_script("util/solver_util.lua")

dim = util.GetParamNumber("-dim", 2)

-- default parameters
if dim == 2 then
	gridName = "grids/elder_quads_8x2.ugx"
	zOutputTransform = 10
	defRefTol = 0.1
	defCoarsenTol = 0.01
elseif dim == 3 then
	gridName = "grids/elder_hex_8x8x2.ugx"
	zOutputTransform = 300
	defRefTol = 2.5
	defCoarsenTol = 0.25
else print("Dimension "..dim.." not supported");	exit(); end


numTimeSteps 	= util.GetParamNumber("-numTimeSteps", 100)
maxLvl			= util.GetParamNumber("-maxLvl",    5)
-- timestep in seconds: 3153600 sec = 0.1 year
--dt 				= util.GetParamNumber("-dt", 3.1536e6)
dt 				= util.GetParamNumber("-dt", 2.e6)
minStepSize     = util.GetParamNumber("-minStepSize", dt/40);
stepReducFac  	= util.GetParamNumber("-stepReducFac", 0.25);
verbose		  	= util.HasParamOption("-verbose");
refTol			= util.GetParamNumber("-refTol", defRefTol, "refinement tolerance for error-indicator. -1 -> no refs")
coarsenTol		= util.GetParamNumber("-coarsenTol", defCoarsenTol, "coarsening tolerance for error-indicator. -1 -> no coarsening")
stepsBetweenRedist = util.GetParamNumber("-stepsBetweenRedist", 1, "Number of timesteps before a new redistribution is attempted - only for parallel use. Default is 1")
numDampSteps	= util.GetParamNumber("-numDampSteps", 3, "Number of timesteps in which dt is damped. Damping is reduced with each time step. Default is 3")
maxDamping		= util.GetParamNumber("-maxDamping", 0.1, "The maximal damping with which dt is damped during damp-steps. Damping is reduced with each time step. Default is 0.1")
discardSolution = util.HasParamOption("-discardSolution", "solution won't be saved to a file.")
writeGridLayout = util.HasParamOption("-writeGridLayout", "adapted grid layouts will be written to a file.")
stepsBetweenWrite = util.GetParamNumber("-stepsBetweenWrite", 1, "Number of timesteps until the next solution or grid-layout is written to a file. Default is 1")
stepsBetweenRedist = math.floor(stepsBetweenRedist)
if(stepsBetweenRedist < 1) then stepsBetweenRedist = 1; end

balancer.partitioner = "dynBisection"
balancer.ParseParameters()
balancer.PrintParameters()
print("")

util.CheckAndPrintHelp()

--------------------------------------------------------------------------------
--  Setup Domain and Approximation Space
--------------------------------------------------------------------------------
InitUG(dim, AlgebraType("CPU", 1))


-- Create a domain for the specified grid file
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateDomain(gridName, 0, neededSubsets)


-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:add_fct("p", "Lagrange", 1)
--approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

--------------------------------------------------------------------------------
--  Setup FV Element Discretization
--------------------------------------------------------------------------------

ug_load_script("elder_user_data.lua")

-- User Data
Porosity = 0.1
MolecularDiffusion = 3.565e-6
Permeability = 4.845e-13
Viscosity = 1e-3;

Gravity = ConstUserVector(0.0)
Gravity:set_entry(dim-1, -9.81)

function DensityFct(c) return 1000 + 200 * c end
function DDensityFct_c(c) return 200 end

-- create dirichlet boundary for concentration
dirichletBND = DirichletBoundary()
dirichletBND:add("ConcentrationDirichletBnd", "c", "Boundary")
dirichletBND:add("PressureDirichletBnd", "p", "Boundary")

-- molecular Diffusion
Diffusion = ScaleAddLinkerMatrix()
Diffusion:add(Porosity, MolecularDiffusion)

-- Density
Density = LuaUserFunctionNumber("DensityFct", 1);
Density:set_deriv(0, "DDensityFct_c");

	-- Darcy Velocity
DarcyVelocity = DarcyVelocityLinker(); 
DarcyVelocity:set_permeability(Permeability)
DarcyVelocity:set_viscosity(Viscosity)
DarcyVelocity:set_density(Density)
DarcyVelocity:set_gravity(Gravity)
print("Darcy Velocity created.")

-- set the product Density * Porosity (only Porosity for Bussinesq)
rhophi = Porosity

FlowEq = ConvectionDiffusion("p", "Inner", "fv1")
FlowEq:set_mass(rhophi)
FlowEq:set_mass_scale(0.0)
FlowEq:set_flux(DarcyVelocity)
print("Flow Equation created.")

TransportEq = ConvectionDiffusion("c", "Inner", "fv1")
--TransportEq:set_upwind(FullUpwind())
TransportEq:set_mass_scale(rhophi)
TransportEq:set_velocity(DarcyVelocity)
TransportEq:set_diffusion(Diffusion)
print("Transport Equation created.")

Density:set_input(0, TransportEq:value())
DarcyVelocity:set_pressure_gradient(FlowEq:gradient())

constraints = OneSideP1Constraints()
--constraints = SymP1Constraints() -- only for fe

-- add Element Discretization to discretization
domainDisc = DomainDiscretization(approxSpace)

domainDisc:add(TransportEq)
domainDisc:add(FlowEq)
domainDisc:add(dirichletBND)
domainDisc:add(constraints)

-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()

--------------------------------------------------------------------------------
--  Solver
--------------------------------------------------------------------------------

solverDesc = {
	type = "newton",

	convCheck = {
		type		= "standard",
		iterations 	= 30,		-- maximum number of iterations
		absolute	= 5e-6,		-- absolut value of defect to be reached; usually 1e-7 - 1e-9
		reduction	= 1e-10,	-- reduction factor of defect to be reached; usually 1e-6 - 1e-8
		verbose		= true			-- print convergence rates if true
	},

	linSolver = 
	{
		type = "bicgstab",
		precond = {
	        type            = "gmg",        -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
	        adaptive        = false,
	        smoother		= "ilu",
	        cycle           = "V",          -- gmg-cycle ["V", "F", "W"]
	        preSmooth       = 3,            -- number presmoothing steps
	        postSmooth      = 3,            -- number postsmoothing steps
	        rap             = false,        -- comutes RAP-product instead of assembling if true
			rim				= false,			-- smooth on surface rim
			emulateFullRefined	= false,	-- emulate full grid (works with rap=true only)
	        baseLevel       = 0,               -- gmg - baselevel
	        gatheredBaseSolverIfAmbiguous = true,
	        baseSolver = "lu",
	        approxSpace = approxSpace
		},

		convCheck = {
	        type            = "standard",
	        iterations      = 100,          -- number of iterations
	        absolute        = 1e-8,        -- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be larger than in newton section)
	        reduction       = 1e-3,         -- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be larger than in newton section)
	        verbose         = true,         -- print convergence rates if true
		}
	}
}

newtonSolver = util.solver.CreateSolver(solverDesc)
newtonSolver:init(op)
linSolver = solverDesc.linSolver.instance

--------------------------------------------------------------------------------
--  Find Start Grid
--------------------------------------------------------------------------------

u = GridFunction(approxSpace)
refiner = HangingNodeDomainRefiner(dom)
loadBalancer = balancer.CreateLoadBalancer(dom)

write("Creating initial grid levels:")
for i = 1,maxLvl do
	--write(" " .. i)
	print("Pre-Adaption " .. i)
	
	-- 1. Intepolate start value
	Interpolate("ConcentrationStart", u, "c", 0.0)
	
	-- 2. estimate error and mark
	--MarkForAdaption_GradientIndicator(refiner, u, "c", 1e-8, 0.5, 0.2, maxLvl-1);
	MarkForAdaption_AbsoluteGradientIndicator(refiner, u, "c", refTol, -1, 0, maxLvl)
	
	-- 3. refine
	refiner:refine() 
	refiner:clear_marks()
	
	if(loadBalancer ~= nil) then
		loadBalancer:rebalance()
		loadBalancer:create_quality_record("initial-redist-"..i..":")
	end
end
--write(" done\n")
print ("Pre-Adaption done")
				
--SaveDebugHierarchy(dom, "initial")
approxSpace:print_statistic()
print(dom:domain_info():to_string())

--------------------------------------------------------------------------------
--  TimeLoop
--------------------------------------------------------------------------------

time = 0.0
step = 0

-- set initial value
print("Interpolation start values")
Interpolate("PressureStart", u, "p", time)
Interpolate("ConcentrationStart", u, "c", time)

-- filename
filename = "Elder"

-- write start solution
if(discardSolution == false) then
	print("Writing start values")
	out = VTKOutput()
	out:select_nodal("c", "c")
	out:select_nodal("p", "p")
	out:select_element(DarcyVelocity, "DarcyVelocity")
	out:print(filename, u, step, time)
end

-- create new grid function for old value
uOld = GridFunction(approxSpace)
uOld:assign(u)

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

local linSolStatFile = io.open("LinSolverStats.txt", "w+")
linSolStatFile:write("# "  .. "Timestep" .. "Num Step" .. " \t " .. "Defect" .. " \t " .. "Rate" .. " \t " .. "Avg rate" .. " \n")		
	if (not linSolStatFile) then
		write("Gnuplot Error: cannot open output file: '")
		write("LinSolverStats.txt" .. " '\n");
		return 1
	end
io.close(linSolStatFile)

numTotalNewtonSteps = 0

for step = 1, numTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGINS at " .. time .. " ++++++")
	
	if(writeGridLayout == true) and (step % stepsBetweenWrite == 0) then
		SaveParallelGridLayout(dom:grid(),
				"layout-step-"..step.."-p"..ProcRank()..".ugx", zOutputTransform)
	end
	
	-- choose time step
	do_dt = dt
	
	-- start with smaller timestep for first steps
	if step <= numDampSteps then
	--	perform linear interpolation between maximal damped time step and normal time step
		local ia = (step - 1) / numDampSteps
		do_dt = do_dt * ia + (1 - ia) * do_dt * maxDamping;
	end

	print("Size of timestep dt: " .. do_dt)
	
	bSuccess = false;
	
	while bSuccess == false do
		-- setup time Disc for old solutions and timestep
		timeDisc:prepare_step(solTimeSeries, do_dt)
		
		-- prepare newton solver
		if newtonSolver:prepare(u) == false then 
			print ("Newton solver failed at step "..step.."."); exit();
		end 
		
		-- apply newton solver
		if newtonSolver:apply(u) == false then 
			do_dt = do_dt * stepReducFac;
			print("\n++++++ Newton solver failed. Trying decreased stepsize " .. do_dt);
			if(do_dt < minStepSize) then
				print("++++++ Time Step to small. Cannot solve problem.");
				exit();
			end
		else
			bSuccess = true; 
		end
		
		numTotalNewtonSteps = numTotalNewtonSteps + newtonSolver:num_newton_steps()
	end
	
	-- linear solver statistics
	linSolStat = {}
	linSolStat.step = linSolver:convergence_check():step()
	linSolStat.defect = linSolver:convergence_check():defect()
	linSolStat.rate = linSolver:convergence_check():reduction()
	linSolStat.avg_rate = linSolver:convergence_check():avg_rate()
	
	print(" ##### LIN SOLVER STATS AT STEP " .. step .. " #####") 
	print(" ## steps: " .. linSolStat.step)
	print(" ## defect: " .. linSolStat.defect)
	print(" ## rate: " .. linSolStat.rate)
	print(" ## avg_red: " .. linSolStat.avg_rate)
	print(" ##############################################")
	
	local linSolStatFile = io.open("LinSolverStats.txt", "a")
	linSolStatFile:write(" " .. step .. " \t " .. linSolStat.step .. " \t " .. linSolStat.defect .. " \t " .. linSolStat.rate .. " \t " .. linSolStat.avg_rate .. " \n")		
	io.close(linSolStatFile)
	
	-- update new time
	time = solTimeSeries:time(0) + do_dt

 --	2: PERFORM GRID ADAPTION
	local curNumVrts = ParallelSum(dom:grid():num_vertices())
	print("  Current number of vertices (global): " .. curNumVrts)
	print("  Coarsening")
	local adptionRunning = true
	local i = 1
	while (adptionRunning == true) do
		MarkForAdaption_AbsoluteGradientIndicator(refiner, u, "c", -1, coarsenTol, 0, maxLvl)
	
		local nodesBefore = ParallelSum(dom:grid():num_vertices())
		refiner:coarsen() 
		local nodesAfter = ParallelSum(dom:grid():num_vertices())
		refiner:clear_marks()
				
		if(nodesAfter == nodesBefore) then
			adptionRunning = false
			print("  New number of vertices left (global): " .. nodesAfter)
		else
			write (" " .. i .. " ")
			i = i + 1
		end
	end
	
	print("  Refining")
	adptionRunning = true
	i = 1
	while (adptionRunning == true) do
		MarkForAdaption_AbsoluteGradientIndicator(refiner, u, "c", refTol, -1, 0, maxLvl)
		
		local nodesBefore = ParallelSum(dom:grid():num_vertices())
		refiner:refine() 
		local nodesAfter = ParallelSum(dom:grid():num_vertices())
		refiner:clear_marks()
		
		if(nodesAfter == nodesBefore) then
			adptionRunning = false
			print("  New number of vertices left (global): " .. nodesAfter)
		else
			write (" " .. i .. " ")
			i = i + 1
		end	
	end

	if((step ~= numTimeSteps) and (loadBalancer ~= nil) and (step % stepsBetweenRedist == 0)) then
		loadBalancer:rebalance()
		loadBalancer:create_quality_record("step-"..step..":")
	end
	
	-- plot solution
	if(discardSolution == false) and (step % stepsBetweenWrite == 0) then
		out:print(filename, u, step, time)
	end
	
	-- get oldest solution
	oldestSol = solTimeSeries:oldest()

	-- copy values into oldest solution (we reuse the memory here)
	VecScaleAssign(oldestSol, 1.0, u)
	
	-- push oldest solutions with new values to front, oldest sol pointer is poped from end
	solTimeSeries:push_discard_oldest(oldestSol, time)

	print("dof statistics")
	approxSpace:print_statistic()
	print("domain info:")
	print(dom:domain_info():to_string())

	print("++++++ TIMESTEP " .. step .. "  END ++++++");
end

-- end timeseries, produce gathering file
if(discardSolution == false) then
	out:write_time_pvd(filename, u)
end

if loadBalancer ~= nil then
	print("Distribution quality statistics:")
	loadBalancer:print_quality_records()
end


print("Solver info:")
print("  num newton steps:      " .. numTotalNewtonSteps)
print("  num lin-solver calls:  " .. newtonSolver:total_linsolver_calls())
print("  num lin-solver steps:  " .. newtonSolver:total_linsolver_steps())
print("")

