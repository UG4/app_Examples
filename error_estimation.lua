--------------------------------------------------------------
-- Solving a Poisson problem on parallel adaptive multigrid --
-- hierarchies using residual-based error estimation.       --
--                                                          --
-- author: Markus Breit                                     --
-- date:   2019-04-02                                       --
--------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- dimension of the problem
dim = util.GetParamNumber("-dim", 2)
if dim ~= 2 and dim ~= 3 then
	print("Dimension " .. dim .. " not supported. Only 2 and 3 are. Aborting.")
	exit()
end
InitUG(dim, AlgebraType("CPU", 1))


-------------------------------------
-- process command line parameters --
-------------------------------------
-- grid file
gridName = util.GetParam("-grid", "Examples/grids/laplace_sample_grid_"..dim.."d.ugx")

-- adaption parameters
maxLvl = util.GetParamNumber("-maxLvl", 15)
refTol = util.GetParamNumber("-refTol", 5e-2, "refinement tolerance")


-------------------------------------------
-- create domain and approximation space --
-------------------------------------------
-- load and check domain
dom = Domain()
LoadDomain(dom, gridName)
if util.CheckSubsets(dom, {"Inner", "Boundary"}) == false then
	print("Subset check failed. Aborting...")
	exit()
end

-- approx space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_surfaces()
approxSpace:init_top_surface()
approxSpace:print_statistic()

-- prepare refiner
refiner = HangingNodeDomainRefiner(dom)


--------------------
-- discretizatiom --
--------------------
-- define source and boundary terms
function source2d(x, y, t)
	local r = math.sqrt(x*x + y*y)
	if r < 1e-9 then return 0.0 end
	return -200*(1.0 + math.exp(-200*r+160) - 200*r + 200*math.exp(-200*r+160)*r)
				* math.exp(-200*r+160) * math.pow(1.0 + math.exp(-200*r+160),-3) / r
end
function dirichletBnd2d(x, y, t)
	local r = math.sqrt(x*x + y*y)
	return true, 1.0 / (1.0 + math.exp(-200*(r-0.8)))
end

function source3d(x, y, z, t)
	local r = math.sqrt(x*x + y*y + z*z)
	if r < 1e-9 then return 0.0 end
	return -400*(1.0 + math.exp(-200*r+160) - 100*r + 100*math.exp(-200*r+160)*r)
				* math.exp(-200*r+160) * math.pow(1.0 + math.exp(-200*r+160),-3) / r
end
function dirichletBnd3d(x, y, z, t)
	local r = math.sqrt(x*x + y*y + z*z)
	return true, 1.0 / (1.0 + math.exp(-200*(r-0.8)))
end


-- set up element discretization and constraints
elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
elemDisc:set_diffusion(1.0)
elemDisc:set_source("source"..dim.."d")

dirichletBND = DirichletBoundary()
dirichletBND:add("dirichletBnd"..dim.."d", "c", "Boundary")

hangingNodeConstraint = SymP1Constraints()


-- setup error estimators
ee = SideAndElemErrEstData(2, 2, "Inner")
eeMult = MultipleSideAndElemErrEstData(approxSpace)
eeMult:add(ee, "c")

elemDisc:set_error_estimator(ee)
dirichletBND:set_error_estimator(eeMult)


-- put together global domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)
domainDisc:add(hangingNodeConstraint)

-- prepare linear operator
linOp = AssembledLinearOperator()
linOp:set_discretization(domainDisc)


--------------------
-- prepare solver --
--------------------
-- GMG setup
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(LU())
gmg:set_gathered_base_solver_if_ambiguous(true)
gmg:set_smoother(GaussSeidel())
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)
gmg:set_smooth_on_surface_rim(true)
gmg:set_cycle_type(1)  -- V-cycle
gmg:set_rap(true)
	
-- linear solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(ConvCheck(100, 1e-12, 1e-12))


------------------
-- solve system --
------------------
-- in parallel environments, use a load balancer to distribute the grid
balancer.partitioner = "parmetis"
loadBalancer = balancer.CreateLoadBalancer(dom)


-- prepare grid functions for solution and rhs
u = GridFunction(approxSpace)
u:set(0.0)
b = GridFunction(approxSpace)

-- prepare adaptive remeshing
refStrat = StdRefinementMarking(refTol, maxLvl)
coarseStrat = StdCoarseningMarking(refTol, 8.0)


-- prepare output of error indicators
approxSpace_vtk = ApproximationSpace(dom)
approxSpace_vtk:add_fct("eta_squared", "piecewise-constant")
u_vtk = GridFunction(approxSpace_vtk)

out_error = VTKOutput()
out_error:clear_selection()
out_error:select_all(false)
out_error:select_element("eta_squared", "error")


-- start adaptation loop
breakAtEndOfLoop = false
for i = 1, maxLvl do
	print(" #######  START adaption " .. i .."  #######")
	
	-- rebalance
	if loadBalancer ~= nil then
		loadBalancer:rebalance()
		loadBalancer:create_quality_record("step "..i..":")
	end
	
	-- init operator
	linOp:init_op_and_rhs(b)
	
	-- set dirichlet values in solution
	linOp:set_dirichlet_values(u)

	-- init solver for linear operator
	if not linSolver:init(linOp) then
		print("Solver initialization failed. Aborting")
		exit()
	end
	
	-- apply solver
	if not linSolver:apply(u,b) then
		print("Solver application failed. Aborting")
		exit()
	end

	-- export solution to vtk	
	WriteGridFunctionToVTK(u, "solution_"..i)
	
	-- estimate error and mark
	domainDisc:calc_error(u, u_vtk)
	WriteGridFunctionToVTK(u_vtk, "error_indicators_"..i)
	if i < maxLvl then
		-- coarsen if possible
		domainDisc:mark_with_strategy(refiner, coarseStrat)
		numElemBeforeCoarsening = dom:domain_info():num_elements()
		numCoarsen = refiner:num_marked_elements()
		if numCoarsen >= numElemBeforeCoarsening/5 then
			print("Coarsening.")
			refiner:coarsen()
			domainDisc:calc_error(u)
		end
		refiner:clear_marks()
		
		-- refine if required
		domainDisc:mark_with_strategy(refiner, refStrat)
		if refiner:num_marked_elements() > 0 then
			print("Refining.")
			refiner:refine()
		else
			print("Solution has reached desired tolerance.")
			breakAtEndOfLoop = true
		end
		refiner:clear_marks()
		
		-- print approximation space statistics
		approxSpace:print_statistic()
	end
	
	print(" #######  END adaption " .. i .."  #######")

	if breakAtEndOfLoop then
		break
	end
end

--[[
-- save final grid hierarchy
offset = 3.0
if dim == 2 then offset = 0.4 end
SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(),
	"grid_hierarchy_transformed_p" .. ProcRank() .. ".ugx", offset)
--]]

-- print distribution quality statistics
if loadBalancer ~= nil then
	loadBalancer:print_quality_records()
end

