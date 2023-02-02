-- Copyright (c) 2010-2017:  G-CSC, Goethe University Frankfurt
-- Authors: Andreas Vogel, Sebastian Reiter
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.


-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")

-- Parse parameters and print help
gridName	= util.GetParam("-grid", "grids/cooler.ugx",
							"filename of underlying grid")
numRefs		= util.GetParamNumber("-numRefs", 2, "number of refinements")

steadyState	= util.HasParamOption("-steadyState", "If specified, the steady state of the problem is computed. Else a time-dependent problem is computed.")

endTime 	= util.GetParamNumber("-endTime", 0.4, "simulated time frame in seconds")
dt			= util.GetParamNumber("-dt", 0.02, "time step size")

util.CheckAndPrintHelp("Cooler");


-- initialize ug with the world dimension 3 and an algebra system with scalar coefficients
InitUG(3, AlgebraType("CPU", 1));


-- Load a domain without initial refinements.
requiredSubsets = {"cooler", "air", "staticAir", "cpu", "inflow", "outflow"}
dom = util.CreateDomain(gridName, 0, requiredSubsets)

-- Refine the domain (redistribution is handled internally for parallel runs)
print("refining...")
-- This balancing setup makes sense for structured grids with uniform refinement
balancerDesc = {
	partitioner = {
		name = "staticBisection",
		clusteredSiblings = false
	},

	hierarchy = {
		name 						= "noRedists",
		minElemsPerProcPerLevel		= redistElemThreshold,
		maxRedistProcs				= redistProcs,
	},
}
util.refinement.CreateRegularHierarchy(dom, numRefs, true, balancerDesc)




-- set up approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("t", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("approximation space:")
approxSpace:print_statistic()


-- set up discretization
coolerDisc = ConvectionDiffusion("t", "cooler", "fv1")
coolerDisc:set_diffusion(10)

airDisc = ConvectionDiffusion("t", "air", "fv1")
airDisc:set_diffusion(0.1)
airDisc:set_velocity({0, 10, 0})
airDisc:set_upwind(FullUpwind())

staticAirDisc = ConvectionDiffusion("t", "staticAir", "fv1")
staticAirDisc:set_diffusion(0.1)

flowBnd = DirichletBoundary()
flowBnd:add(20, "t", "inflow")
flowBnd:add(20, "t", "outflow")

cpuBnd = DirichletBoundary()
cpuBnd:add(80, "t", "cpu")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(coolerDisc)
domainDisc:add(airDisc)
domainDisc:add(staticAirDisc)
domainDisc:add(flowBnd)
domainDisc:add(cpuBnd)


-- set up solver (using 'util/solver_util.lua')
solverDesc = {
	type = "bicgstab",
	precond = {
		type		= "gmg",
		approxSpace	= approxSpace,
		smoother	= "ilu",
		baseSolver	= "lu"
	}
}

solver = util.solver.CreateSolver(solverDesc)


print("\nsolving...")
u = GridFunction(approxSpace)
u:set(20.0)


if steadyState then
	local A = AssembledLinearOperator(domainDisc)
	local b = GridFunction(approxSpace)
	domainDisc:adjust_solution(u)
	domainDisc:assemble_linear(A, b)

	solver:init(A, u)
	solver:apply(u, b)

	solFileName = "sol_cooler"
	print("writing solution to '" .. solFileName .. "'...")
	WriteGridFunctionToVTK(u, solFileName)
	SaveVectorForConnectionViewer(u, solFileName .. ".vec")
else
	local startTime = 0
	util.SolveLinearTimeProblem(u, domainDisc, solver, VTKOutput(), "sol_cooler",
								"ImplEuler", 1, startTime, endTime, dt); 
end

print("done")

