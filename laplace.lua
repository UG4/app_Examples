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
dim			= util.GetParamNumber("-dim", 3, "Dimension of the problem", {2,3})

gridName	= util.GetParam("-grid", "grids/laplace_sphere_" .. dim .. "d.ugx",
							"filename of underlying grid")
numRefs		= util.GetParamNumber("-numRefs", 4, "number of refinements")


util.CheckAndPrintHelp("Laplace-Problem");


-- initialize ug with the world dimension 2 and scalar matrix coefficients
InitUG(dim, AlgebraType("CPU", 1));


-- Load a domain without initial refinements.
requiredSubsets = {"Inner", "bndNegative", "bndPositive"}
dom = util.CreateDomain(gridName, 0, requiredSubsets)

-- Refine the domain (redistribution is handled internally for parallel runs)
print("refining...")
util.refinement.CreateRegularHierarchy(dom, numRefs, true)

-- set up approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("approximation space:")
approxSpace:print_statistic()


-- set up discretization
-- Please have a look at this page for more information on the
-- ConvectionDiffusion discretization object:
-- http://ug4.github.io/docs/plugins/classug_1_1_convection_diffusion_plugin_1_1_convection_diffusion_base.html#details

elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
elemDisc:set_diffusion(1.0)
elemDisc:set_source(0)

dirichletBND = DirichletBoundary()
dirichletBND:add(-1, "c", "bndNegative")
dirichletBND:add(1, "c", "bndPositive")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)


-- set up solver (using 'util/solver_util.lua')
solverDesc = {
	type = "bicgstab",
	precond = {
		type		= "gmg",
		approxSpace	= approxSpace,
		smoother	= "jac",
		baseSolver	= "lu"
	}
}

solver = util.solver.CreateSolver(solverDesc)


print("\nsolving...")
A = AssembledLinearOperator(domainDisc)
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(0.0)
domainDisc:adjust_solution(u)
domainDisc:assemble_linear(A, b)

solver:init(A, u)
solver:apply(u, b)


solFileName = "sol_laplace_" .. dim .. "d"
print("writing solution to '" .. solFileName .. "'...")
WriteGridFunctionToVTK(u, solFileName)
SaveVectorForConnectionViewer(u, solFileName .. ".vec")

print("done")
