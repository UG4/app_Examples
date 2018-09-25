-- Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
-- Authors: Raphael Prohl, Sebastian Reiter
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

-- run this script like this:
-- 'ugshell -ex Examples/solmech-simple-deformation.lua'
-- or use a different geometry:
-- 'ugshell -ex Examples/solmech-simple-deformation.lua -grid grids/solmech_cylinder.ugx'

PluginRequired("SmallStrainMechanics")

-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")
ug_load_script("util/profiler_util.lua")

-- Parse parameters and print help
params = {}
params.gridName	= util.GetParam("-grid", "grids/springboard.ugx", "filename of underlying grid")
params.numRefs	= util.GetParamNumber("-numRefs", 2, "number of refinements")
params.order	= util.GetParamNumber("-order", 1, "order of the function space")
params.smoother	= util.GetParam("-smoother", "gs", "Smoother for the multigrid method. Options are 'gs' and 'ilu'.")


util.CheckAndPrintHelp("Solid Mechanics: Simple Deformation");


-- initialize ug with the world dimension 2 and scalar matrix coefficients
InitUG(3, AlgebraType("CPU", 1));


-- Load a domain without initial refinements.
requiredSubsets = {"Flex", "Force", "Fixed"}
dom = util.CreateDomain(params.gridName, 0, requiredSubsets)


-- Refine the domain (redistribution is handled internally for parallel runs)
print("refining...")
util.refinement.CreateRegularHierarchy(dom, params.numRefs, true)


-- set up approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("ux", "Lagrange", params.order)          
approxSpace:add_fct("uy", "Lagrange", params.order)          
approxSpace:add_fct("uz", "Lagrange", params.order)

approxSpace:init_levels()
approxSpace:init_top_surface()

print("discretization:")
approxSpace:print_statistic()


-- set up discretization
-- coefficients are made up to get a reasonable displacement of the springboard
matLaw = HookeLaw()
matLaw:set_hooke_elasticity_tensor(4000000,3200000)

elemDisc = SmallStrainMechanics("ux, uy, uz", "Flex")
elemDisc:set_material_law(matLaw)
elemDisc:set_mass_scale(0)
elemDisc:set_quad_order(2)

print(elemDisc:config_string())

dirichletBND = DirichletBoundary()
-- dirichletBND:add( 0, "ux", "Force")
-- dirichletBND:add( 0, "uy", "Force")
-- dirichletBND:add( -0.5, "uz", "Force")

dirichletBND:add( 0, "ux", "Fixed")
dirichletBND:add( 0, "uy", "Fixed")
dirichletBND:add( 0, "uz", "Fixed")

neumannBND = NeumannBoundaryFE("uz")
neumannBND:add(600, "Force", "Flex") -- force a human excerts on the springboard while standing


domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)
domainDisc:add(neumannBND)

A = AssembledLinearOperator(domainDisc)
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(0.0)
domainDisc:adjust_solution(u)
domainDisc:assemble_linear(A, b)


--solve
print("\nsolving...")

-- set up solver (using 'util/solver_util.lua')
smoothers = {
	gs = {type = "gs", consistentInterfaces = true},
	ilu = {type = "ilu", overlap = true},
}

solverDesc = {
	type = "bicgstab",
	precond = {
		type		= "gmg",
		approxSpace	= approxSpace,
		smoother	= smoothers[params.smoother],
		baseSolver	= "lu"
	}
}

solver = util.solver.CreateSolver(solverDesc)
solver:init(A, u)
solver:apply(u, b)


--write results
solFileName = "sol_simple_deformation"
print("writing solution to '" .. solFileName .. "'...")

vtkWriter = VTKOutput()
vtkWriter:select_all(false)	-- don't write ux, uy, uz separately
vtkWriter:select_nodal("ux,uy,uz", "displacement")	-- write a displacement vector (ux, uy, uz)
vtkWriter:print(solFileName, u, 0, 0, false)

SaveVectorForConnectionViewer(u, solFileName .. ".vec")


if GetProfilerAvailable() == true then
	WriteProfileData("solmech-simple-deformation.pdxml")
	print("")
	print("--- PROFILES --- ")
	util.PrintProfile_TotalTime("main                           ")
	util.PrintProfile_TotalTime("  assemble_linear              ")
	util.PrintProfile_TotalTime("  LS_InitPrecond               ")
	util.PrintProfile_TotalTime("  LS_ApplyReturnDefect         ")
end

print("done")