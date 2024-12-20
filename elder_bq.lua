----------------------------------------------------------
--
--  Lua-script for simulation of the Elder problem (with boussinesq approximation)
--
--  Author: Ekaterina Vasilyeva
--
--  Originally implemented for AMFS 2021
--
--  Based on Examples/elder_adapt.lua
----------------------------------------------------------

ug_load_script ("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")

------------------------------------------------------------------------------------------
-- Geometry and time interval
------------------------------------------------------------------------------------------

dim			= 2
gridName	= "grids/elder_quads_8x2_bnd.ugx"

numRefs 	= 5 -- number of refinements of the initial grid
numPreRefs 	= 1 -- this grid level is considered as the most coarse one for the numerics
endTime		= 10 * 365 * 24 * 60 * 60 -- [s] = 10 years
dt			= 10 * 24 * 60 * 60 -- [s] = 10 days

upwind		= "partial" -- upwind type for the transport equation: "no", "full" or "partial"

vtk_file_name = "Elder_Bq_GL" .. numRefs -- VTK output file name base

------------------------------------------------------------------------------------------
-- Geological parameters, initial and boundary conditions
------------------------------------------------------------------------------------------

Porosity = 0.1
MolecularDiffusion = 3.565e-6
Permeability = 4.845e-13
Viscosity = 1e-3

function PressureStart(x, y, t, si)
	return 9810 * (150 - y)
end

------------------------------------------------------------------------------------------
-- Dirichlet boundary conditions for the salt on the top boundary
------------------------------------------------------------------------------------------


function ConcentrationDirichletBnd(x, y, t)
		if x > 150 and x < 450 then
			return true, 1.0
		end
	
	    return false, 0.0
end

------------------------------------------------------------------------------------------
-- Initialize ug4; use the block algebra for the pairs (c, p) (2 components per block)
------------------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 2));

--------------------------------------------------------------------------------
--  Setup Domain and Approximation Space
--------------------------------------------------------------------------------

-- create the domain, load the grid and refine it
dom = util.CreateDomain(gridName, numPreRefs, {"Inner", "LeftRightBnd", "BottomBnd", "TopBnd","TopCorners"})
balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

-- set up approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:add_fct("p", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
--approxSpace:print_layout_statistic() -- makes sense only for parallel runs
approxSpace:print_statistic()

-- Order the DoFs:
OrderLex (approxSpace, "y")

--------------------------------------------------------------------------------
--  Setup the spatial FV Discretization
--------------------------------------------------------------------------------

Gravity = ConstUserVector(0.0)
Gravity:set_entry(1, -9.81)

-- create dirichlet boundary for concentration
dirichletBND = DirichletBoundary()
dirichletBND:add("ConcentrationDirichletBnd", "c", "TopBnd")
dirichletBND:add(0, "c", "BottomBnd")
-- TopCorners
dirichletBND:add(0, "p", "TopCorners")

-- density
function DensityFct(c) return 1000 + 200 * c end
function DDensityFct_c(c) return 200 end
Density = LuaUserFunctionNumber("DensityFct", 1);
Density:set_deriv(0, "DDensityFct_c");

-- Darcy Velocity (as a linker)
DarcyVelocity = DarcyVelocityLinker(); 
DarcyVelocity:set_permeability(Permeability)
DarcyVelocity:set_viscosity(Viscosity)
DarcyVelocity:set_density(Density)
DarcyVelocity:set_gravity(Gravity)

-- flow equation
FlowEq = ConvectionDiffusionFV1("p", "Inner")
FlowEq:set_mass_scale(0.0)
FlowEq:set_flux(DarcyVelocity)

-- transport equation
TransportEq = ConvectionDiffusionFV1("c", "Inner")
TransportEq:set_upwind(UpwindFV1(upwind))
TransportEq:set_mass_scale(Porosity)
TransportEq:set_velocity(DarcyVelocity)
TransportEq:set_diffusion(Porosity * MolecularDiffusion)

-- associate the density with the unknown of the transport equation (i.e. c)
Density:set_input(0, TransportEq:value())
-- associate the Darcy velocity with the unknown of the flow equation (i.e. p)
DarcyVelocity:set_pressure_gradient(FlowEq:gradient())

-- set up the global discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(TransportEq)
domainDisc:add(FlowEq)
domainDisc:add(dirichletBND)

--------------------------------------------------------------------------------
--  Set up the solver
--------------------------------------------------------------------------------

util.solver.defaults.approxSpace = approxSpace
solverDesc = {
	type = "newton",
	
	lineSearch = {
		type			= "standard",-- ["standard", "none"]
		maxSteps		= 8,		-- maximum number of line search steps
		lambdaStart		= 1,		-- start value for scaling parameter
		lambdaReduce	= 0.5,		-- reduction factor for scaling parameter
		acceptBest 		= true,		-- check for best solution if true
		checkAll		= false		-- check all maxSteps steps if true 
	},
	
	convCheck = {
		type			= "standard",
		iterations 		= 30,		-- maximum number of iterations
		absolute		= 5e-6,		-- absolut value of defect to be reached; usually 1e-7 - 1e-9
		reduction		= 1e-10,	-- reduction factor of defect to be reached; usually 1e-6 - 1e-8
		verbose			= true		-- print convergence rates if true
	},

	linSolver = 
	{
		type = "bicgstab",
		precond = {
	        type		= "gmg",	-- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
	        smoother	= "ilu",	-- pre- and postsmoother, only for gmg ["ilu", "ilut", "jac", "gs", "sgs"]
	        cycle		= "V",		-- gmg-cycle ["V", "F", "W"]
	        preSmooth	= 2,		-- number presmoothing steps
	        postSmooth	= 2,		-- number postsmoothing steps
	        rap			= false,	-- comutes RAP-product instead of assembling if true
	        baseLevel	= numPreRefs,-- gmg - coarsest level
	        baseSolver	= "lu",
		},

		convCheck = {
	        type		= "standard",
	        iterations	= 100,		-- number of iterations
	        absolute	= 1e-8,		-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (may not be larger than in newton section)
	        reduction	= 1e-5,		-- reduction factor of defect to be reached
	        verbose		= true,		-- print convergence rates if true
		}
	}
}

solver = util.solver.CreateSolver(solverDesc)

-----------------------------------------------------------------------------------------
-- Ploter for the results
------------------------------------------------------------------------------------------

out = VTKOutput()
out:clear_selection()
out:select_nodal("c", "c")
out:select_nodal("p", "p")
out:select(DarcyVelocity, "q")

------------------------------------------------------------------------------------------
-- Set up and apply the time stepping scheme
------------------------------------------------------------------------------------------

-- Set up the time discretization and the time stepping scheme
timeDisc = ThetaTimeStep(domainDisc, 1.0) -- Implicit Euler: 1.0
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(solver)
timeIntegrator:attach_observer(VTKOutputObserver(vtk_file_name, out))

timeIntegrator:set_time_step(dt)

------------------------------------------------------------------------------------------
-- Prepare the initial guess for the pressure
------------------------------------------------------------------------------------------

-- grid function for the solution
u = GridFunction(approxSpace)

-- Set the initial condition
Interpolate("PressureStart", u, "p")
Interpolate(0, u, "c")

-- Fix the mass fraction and solve the linear problem for the pressure

fixer = DirichletBoundary()
domainDisc:add(fixer)
fixer:invert_subset_selection()
fixer:add("c", "")

solver:init(AssembledOperator(domainDisc))
solver:prepare(u)

-- apply the solver for the stationary pressure problem
if not solver:apply(u) then
	print("===> THE PREPARATION PHASE FAILED! <===")
	exit()
end

domainDisc:remove (fixer)

------------------------------------------------------------------------------------------
-- Compute the time steps
------------------------------------------------------------------------------------------

timeIntegrator:apply(u, endTime, u, 0.0)

------------------------------------------------------------------------------------------
-- End of File
------------------------------------------------------------------------------------------
