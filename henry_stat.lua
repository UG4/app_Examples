----------------------------------------------------------
--
--  Lua-script for simulation of the Henry problem (with boussinesq approximation)
--
--  Author: Dmitry Logashenko
--
----------------------------------------------------------
PrintBuildConfiguration()
ug_load_script("ug_util.lua") 
ug_load_script ("util/load_balancing_util.lua")

dim = 2 -- dimensionality of the problem
InitUG (dim, AlgebraType ("CPU", 2)) -- use the block algebra (2x2 blocks for the point-block discr.)

------------------------------------------------------------------------------------------
-- Geological parameters, initial and boundary conditions
------------------------------------------------------------------------------------------

phi = 0.35 -- porosity
D = 18.8571e-6 * phi -- constant diffusion/dispersion
K = 1.019368e-9 -- permeability
mu = 1e-3

recharge = 6.6e-5 -- assumed to be a constant here (or 3.3e-5)

boussinesq = true -- if to use the Boussinesq approximation

------------------------------------------------------------------------------------------
-- Geometry and parameters of the numerical methods
------------------------------------------------------------------------------------------

gridName	= "grids/henry_quad_2x1.ugx" -- grid of quadrilaterals
--gridName	= "grids/henry_tri_4.ugx" -- alternative grid (of triangles)

numRefs 	= 7 -- finest grid level
numPreRefs 	= 3 -- coarse grid level for GMG

upwind		= "full" -- upwind type for the transport equation: "no", "full" or "partial"

vtk_file_name = "HenryStat" -- VTK output file name base
if boussinesq then vtk_file_name = vtk_file_name .. "-Boussinesq" end
vtk_file_name = vtk_file_name .. "-GL" .. numRefs

-- density
function rho_fct (c) return 1000 + 25 * c end
function Drho_fct_c(c) return 25 end
rho = LuaUserFunctionNumber (rho_fct, 1);
rho:set_deriv (0, Drho_fct_c);

-- gravity
g_val = -9.81
g = ConstUserVector (0)
g:set_entry (1, g_val)

------------------------------------------------------------------------------------------
-- Dirichlet boundary conditions for the pressure and the salt mass fraction
------------------------------------------------------------------------------------------

function hydroStatPressure (x, y) -- hydrostatic pressure used at the sea side
	return g_val * rho_fct (1) * y
end

------------------------------------------------------------------------------------------
-- Initial guesses for the pressure and the mass fraction
------------------------------------------------------------------------------------------

-- initMassFrac = function (x, y) if 2 * (x - 2) > y then return 1 end return 0 end

initPressure = hydroStatPressure

--------------------------------------------------------------------------------
--  Domain and Approximation Space
--------------------------------------------------------------------------------

-- load and refine the grid
dom = util.CreateDomain (gridName, numPreRefs, {"Medium", "Top", "Bottom", "Sea", "Inflow"})
balancer.RefineAndRebalanceDomain (dom, numRefs - numPreRefs)

print ("Domain info:")
print (dom:domain_info():to_string ())

-- create the approximation space
approxSpace = ApproximationSpace (dom); util.solver.defaults.approxSpace = approxSpace
approxSpace:add_fct ("c", "Lagrange", 1)
approxSpace:add_fct ("p", "Lagrange", 1)
approxSpace:init_levels ()
approxSpace:init_top_surface ()

print("---> DoF statistics:")
approxSpace:print_statistic()

-- Order the DoFs:
OrderLex (approxSpace, "y")

--------------------------------------------------------------------------------
--  Spatial Discretization (vertex-centered FV)
--------------------------------------------------------------------------------

domainDisc = DomainDiscretization(approxSpace); assembledOp = AssembledOperator (domainDisc)

-- Dirichlet BC
DirichletBC = DirichletBoundary (); domainDisc:add (DirichletBC)
-- Sea side
DirichletBC:add (1, "c", "Sea")
DirichletBC:add (hydroStatPressure, "p", "Sea")
-- Inflow
DirichletBC:add (0, "c", "Inflow")

-- Neumann BC
NeumannBC = NeumannBoundaryFV1 ("p"); domainDisc:add (NeumannBC)
NeumannBC:add (- recharge * rho_fct (0), "Inflow", "Medium")

-- Darcy Velocity (as a linker)
DarcyVel = DarcyVelocityLinker ()
DarcyVel:set_permeability (K)
DarcyVel:set_viscosity (mu)
DarcyVel:set_density (rho)
DarcyVel:set_gravity (g)

-- The density as a factor in the fluxes
if boussinesq then
	rho_factor = rho_fct (0) -- the density is considered to be constant in this factor
else
	rho_factor = rho
end

-- Continuity equation
ContinuityEq = ConvectionDiffusionFV1 ("p", "Medium"); domainDisc:add (ContinuityEq)
ContinuityEq:set_flux (rho_factor * DarcyVel)

-- Transport equation
TransportEq = ConvectionDiffusionFV1 ("c", "Medium"); domainDisc:add (TransportEq)
TransportEq:set_upwind (UpwindFV1 (upwind))
TransportEq:set_velocity (rho_factor * DarcyVel)
TransportEq:set_diffusion (rho_factor * ToUserMatrix (D))

-- connect the density with the mass fraction (as the unknown of the transport eq.)
rho:set_input (0, TransportEq:value ())
-- connect the velocity with the pressure gradient (as the gradient of the unknown of the continuity eq.)
DarcyVel:set_pressure_gradient (ContinuityEq:gradient ())

--------------------------------------------------------------------------------
--  The non-linear solver
--------------------------------------------------------------------------------
solver = util.solver.CreateSolver (
{
	type = "newton",
	
	lineSearch = {
		type			= "standard",-- ["standard", "none"]
		maxSteps		= 20,		-- maximum number of line search steps
		lambdaReduce	= 0.5,		-- reduction factor for scaling parameter
		acceptBest 		= true,		-- check for best solution if true
		checkAll		= false		-- check all maxSteps steps if true 
	},
	
	convCheck = {
		type			= "standard",
		iterations 		= 30,		-- maximum number of iterations
		absolute		= 1e-9,	-- absolut value of defect to be reached
		reduction		= 1e-6,		-- reduction factor of defect to be reached
		verbose			= true		-- print convergence rates if true
	},

	linSolver = 
	{
		type = "bicgstab",
		precond = {
	        type		= "gmg",	-- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
	        smoother	= {			-- pre- and postsmoother, only for gmg
	        	type = "ilu",		-- type of the smoothers ["ilu", "ilut", "jac", "gs", "sgs"]
	        	beta = 0.75,		-- parameters of the smoothers [beta for ilu only]
	        	overlap = true		-- parallelization parameters
	        },
	        cycle		= "V",		-- gmg-cycle ["V", "F", "W"]
	        preSmooth	= 4,		-- number presmoothing steps
	        postSmooth	= 4,		-- number postsmoothing steps
	        rap			= true,	-- comutes RAP-product instead of assembling if true
	        baseLevel	= numPreRefs,-- gmg - coarsest level
	        baseSolver	= "lu",
		},

		convCheck = {
	        type		= "standard",
	        iterations	= 100,		-- number of iterations
	        absolute	= 0.5e-9,	-- absolut value of defact to be reached
	        reduction	= 1e-4,		-- reduction factor of defect to be reached
	        verbose		= true,		-- print convergence rates if true
		}
	}
})

------------------------------------------------------------------------------------------
-- Prepare the initial guess
------------------------------------------------------------------------------------------

-- grid function for the solution
u = GridFunction (approxSpace)
u:set (0)
-- Interpolate (initMassFrac, u, "c")
Interpolate (hydroStatPressure, u, "p")

-- Fix the mass fraction and solve the linear problem for the pressure

fixer = DirichletBoundary (); domainDisc:add (fixer)
fixer:invert_subset_selection ()
fixer:add ("c", "")

solver:init (assembledOp)
solver:prepare (u)

-- apply the solver
print ("===> Computation of the initial guess for the pressure <===")
if not solver:apply (u) then
	print ("===> THE PREPARATION PHASE FAILED! <===")
	exit ()
end

domainDisc:remove (fixer)

------------------------------------------------------------------------------------------
-- Solve the stationary problem
------------------------------------------------------------------------------------------

-- initialize the solver
solver:init (assembledOp)
solver:prepare (u)

-- apply the solver
print ("===> Computation of the solution of the problem <===")
if not solver:apply (u) then
	print ("===> SOLVER FAILED! <===")
	exit ()
end

------------------------------------------------------------------------------------------
-- Output the results
------------------------------------------------------------------------------------------

out = VTKOutput (); out:clear_selection ()
out:select_nodal ("c", "c")
out:select_nodal ("p", "p")
out:select (DarcyVel, "q")
out:print (vtk_file_name, u)

------------------------------------------------------------------------------------------
-- End of File
------------------------------------------------------------------------------------------
